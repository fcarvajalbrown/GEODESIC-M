use crate::device::GpuContext;
use geodesic_core::BackendError;
use wgpu::util::DeviceExt;

pub struct NonbondedInput<'a> {
    pub pos_x: &'a [f64],
    pub pos_y: &'a [f64],
    pub pos_z: &'a [f64],
    pub sigma: &'a [f64],
    pub epsilon: &'a [f64],
    pub offsets: &'a [u32],
    pub neighbors: &'a [u32],
    pub r_cutoff: f64,
    pub r_switch: f64,
    pub box_size: [f64; 3],
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct ParamsUniform {
    n_atoms: u32,
    r_cutoff: f32,
    r_switch: f32,
    _pad0: f32,
    box_size: [f32; 3],
    _pad1: f32,
}

pub struct NonbondedKernel {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

fn storage_entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

impl NonbondedKernel {
    pub fn new(ctx: &GpuContext) -> Result<Self, BackendError> {
        let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("nonbonded"),
            source: wgpu::ShaderSource::Wgsl(include_str!("nonbonded.wgsl").into()),
        });
        let entries = vec![
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            storage_entry(1, true),
            storage_entry(2, true),
            storage_entry(3, true),
            storage_entry(4, true),
            storage_entry(5, true),
            storage_entry(6, false),
        ];
        let layout = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("nonbonded"),
            entries: &entries,
        });
        let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("nonbonded"),
            bind_group_layouts: &[&layout],
            push_constant_ranges: &[],
        });
        let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("nonbonded"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "main",
            compilation_options: Default::default(),
            cache: None,
        });
        Ok(Self { pipeline, layout })
    }

    pub fn evaluate(&self, ctx: &GpuContext, input: &NonbondedInput) -> (Vec<[f32; 3]>, f32) {
        let n = input.pos_x.len();
        let positions: Vec<[f32; 4]> = (0..n)
            .map(|i| [input.pos_x[i] as f32, input.pos_y[i] as f32, input.pos_z[i] as f32, 0.0])
            .collect();
        let sigma: Vec<f32> = input.sigma.iter().map(|&x| x as f32).collect();
        let epsilon: Vec<f32> = input.epsilon.iter().map(|&x| x as f32).collect();
        let params = ParamsUniform {
            n_atoms: n as u32,
            r_cutoff: input.r_cutoff as f32,
            r_switch: input.r_switch as f32,
            _pad0: 0.0,
            box_size: [input.box_size[0] as f32, input.box_size[1] as f32, input.box_size[2] as f32],
            _pad1: 0.0,
        };

        let dev = &ctx.device;
        let mk_storage = |data: &[u8], label: &str| {
            dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: data,
                usage: wgpu::BufferUsages::STORAGE,
            })
        };
        let params_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let pos_buf = mk_storage(bytemuck::cast_slice(&positions), "positions");
        let sig_buf = mk_storage(bytemuck::cast_slice(&sigma), "sigma");
        let eps_buf = mk_storage(bytemuck::cast_slice(&epsilon), "epsilon");
        let off_buf = mk_storage(bytemuck::cast_slice(input.offsets), "offsets");
        let nbr_buf = mk_storage(bytemuck::cast_slice(input.neighbors), "neighbors");

        let out_len = (n * 4 * 4) as u64;
        let out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out_force"),
            size: out_len,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let read_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: out_len,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind = dev.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("nonbonded"),
            layout: &self.layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: pos_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: sig_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: eps_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: off_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: nbr_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: out_buf.as_entire_binding() },
            ],
        });

        let mut enc = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: None });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind, &[]);
            let groups = (n as u32).div_ceil(64);
            pass.dispatch_workgroups(groups.max(1), 1, 1);
        }
        enc.copy_buffer_to_buffer(&out_buf, 0, &read_buf, 0, out_len);
        ctx.queue.submit(Some(enc.finish()));

        let slice = read_buf.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        dev.poll(wgpu::Maintain::Wait);
        let data = slice.get_mapped_range();
        let raw: &[[f32; 4]] = bytemuck::cast_slice(&data);
        let mut forces = Vec::with_capacity(n);
        let mut energy = 0.0f32;
        for v in raw.iter() {
            forces.push([v[0], v[1], v[2]]);
            energy += v[3];
        }
        drop(data);
        read_buf.unmap();
        (forces, 0.5 * energy)
    }
}

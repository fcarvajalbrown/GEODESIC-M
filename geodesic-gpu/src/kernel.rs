use crate::device::GpuContext;
use geodesic_core::{AtomData, BackendError, SimParams};
use wgpu::util::DeviceExt;

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
    n_atoms: usize,
    params_buf: wgpu::Buffer,
    pos_buf: wgpu::Buffer,
    sig_buf: wgpu::Buffer,
    eps_buf: wgpu::Buffer,
    off_buf: wgpu::Buffer,
    nbr_buf: wgpu::Buffer,
    nbr_capacity: usize,
    out_buf: wgpu::Buffer,
    read_buf: wgpu::Buffer,
    bind: wgpu::BindGroup,
    pos_scratch: Vec<[f32; 4]>,
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

#[allow(clippy::too_many_arguments)]
fn make_bind(
    device: &wgpu::Device,
    layout: &wgpu::BindGroupLayout,
    params_buf: &wgpu::Buffer,
    pos_buf: &wgpu::Buffer,
    sig_buf: &wgpu::Buffer,
    eps_buf: &wgpu::Buffer,
    off_buf: &wgpu::Buffer,
    nbr_buf: &wgpu::Buffer,
    out_buf: &wgpu::Buffer,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("nonbonded"),
        layout,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: pos_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: sig_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: eps_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: off_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: nbr_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: out_buf.as_entire_binding() },
        ],
    })
}

impl NonbondedKernel {
    pub fn new(ctx: &GpuContext, atoms: &AtomData, params: &SimParams) -> Result<Self, BackendError> {
        let dev = &ctx.device;
        let shader = dev.create_shader_module(wgpu::ShaderModuleDescriptor {
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
        let layout = dev.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("nonbonded"),
            entries: &entries,
        });
        let pipeline_layout = dev.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("nonbonded"),
            bind_group_layouts: &[&layout],
            push_constant_ranges: &[],
        });
        let pipeline = dev.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("nonbonded"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "main",
            compilation_options: Default::default(),
            cache: None,
        });

        let n = atoms.mass.len();
        let params_uniform = ParamsUniform {
            n_atoms: n as u32,
            r_cutoff: params.r_cutoff as f32,
            r_switch: params.r_switch as f32,
            _pad0: 0.0,
            box_size: [params.box_size[0] as f32, params.box_size[1] as f32, params.box_size[2] as f32],
            _pad1: 0.0,
        };
        let sigma: Vec<f32> = atoms.sigma.iter().map(|&x| x as f32).collect();
        let epsilon: Vec<f32> = atoms.epsilon.iter().map(|&x| x as f32).collect();

        let params_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("params"),
            contents: bytemuck::bytes_of(&params_uniform),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let sig_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("sigma"),
            contents: bytemuck::cast_slice(&sigma),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let eps_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("epsilon"),
            contents: bytemuck::cast_slice(&epsilon),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let pos_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("positions"),
            size: (n * 16) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let off_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("offsets"),
            size: ((n + 1) * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let nbr_capacity = 1usize;
        let nbr_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("neighbors"),
            size: (nbr_capacity * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out_force"),
            size: (n * 16) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let read_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: (n * 16) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let bind = make_bind(dev, &layout, &params_buf, &pos_buf, &sig_buf, &eps_buf, &off_buf, &nbr_buf, &out_buf);

        Ok(Self {
            pipeline,
            layout,
            n_atoms: n,
            params_buf,
            pos_buf,
            sig_buf,
            eps_buf,
            off_buf,
            nbr_buf,
            nbr_capacity,
            out_buf,
            read_buf,
            bind,
            pos_scratch: vec![[0.0f32; 4]; n],
        })
    }

    pub fn upload_neighbors(&mut self, ctx: &GpuContext, offsets: &[u32], neighbors: &[u32]) {
        ctx.queue.write_buffer(&self.off_buf, 0, bytemuck::cast_slice(offsets));
        if neighbors.len() > self.nbr_capacity {
            self.nbr_capacity = neighbors.len();
            self.nbr_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("neighbors"),
                size: (self.nbr_capacity * 4) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.bind = make_bind(
                &ctx.device,
                &self.layout,
                &self.params_buf,
                &self.pos_buf,
                &self.sig_buf,
                &self.eps_buf,
                &self.off_buf,
                &self.nbr_buf,
                &self.out_buf,
            );
        }
        if !neighbors.is_empty() {
            ctx.queue.write_buffer(&self.nbr_buf, 0, bytemuck::cast_slice(neighbors));
        }
    }

    pub fn evaluate(
        &mut self,
        ctx: &GpuContext,
        pos_x: &[f64],
        pos_y: &[f64],
        pos_z: &[f64],
    ) -> (Vec<[f32; 3]>, f32) {
        let n = self.n_atoms;
        for i in 0..n {
            self.pos_scratch[i] = [pos_x[i] as f32, pos_y[i] as f32, pos_z[i] as f32, 0.0];
        }
        ctx.queue.write_buffer(&self.pos_buf, 0, bytemuck::cast_slice(&self.pos_scratch));

        let dev = &ctx.device;
        let mut enc = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: None });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &self.bind, &[]);
            let groups = (n as u32).div_ceil(64);
            pass.dispatch_workgroups(groups.max(1), 1, 1);
        }
        let out_len = (n * 16) as u64;
        enc.copy_buffer_to_buffer(&self.out_buf, 0, &self.read_buf, 0, out_len);
        ctx.queue.submit(Some(enc.finish()));

        let slice = self.read_buf.slice(..);
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
        self.read_buf.unmap();
        (forces, 0.5 * energy)
    }
}

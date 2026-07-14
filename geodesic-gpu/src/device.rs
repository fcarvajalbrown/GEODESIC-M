use geodesic_core::BackendError;

pub struct GpuContext {
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
}

pub fn try_new() -> Result<GpuContext, BackendError> {
    pollster::block_on(async {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::DX12 | wgpu::Backends::VULKAN,
            ..Default::default()
        });
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await
            .ok_or(BackendError::NoAdapter)?;
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("geodesic-gpu"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::Performance,
                },
                None,
            )
            .await
            .map_err(|_| BackendError::DeviceLost)?;
        Ok(GpuContext { device, queue })
    })
}

pub fn context_or_skip() -> Option<GpuContext> {
    match try_new() {
        Ok(ctx) => Some(ctx),
        Err(BackendError::NoAdapter) => {
            eprintln!("skipping GPU test: no adapter (DX12/Vulkan) available");
            None
        }
        Err(e) => {
            eprintln!("skipping GPU test: {e}");
            None
        }
    }
}

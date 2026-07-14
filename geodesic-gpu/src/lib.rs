//! geodesic-gpu: wgpu compute backend for the non-bonded LJ force loop (M2).

pub mod device;
pub mod gpu_backend;
pub mod kernel;
pub mod neighbor_csr;

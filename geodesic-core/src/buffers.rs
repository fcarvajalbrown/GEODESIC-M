/// One per thread, length N, allocated once at startup and zeroed at the
/// start of each step. The master thread reduces all buffers into
/// SimState::force_{x,y,z} in thread-index order (deterministic).
#[derive(Debug, Clone)]
pub struct ForceBuffer {
    pub fx: Vec<f64>, // length N
    pub fy: Vec<f64>,
    pub fz: Vec<f64>,
}

/// Snapshot written to the ring buffer every frame_interval steps.
/// Positions are downcast to f32 here — the simulation itself stays f64.
#[derive(Debug)]
pub struct TrajectoryFrame {
    pub step: u64,
    pub time_ps: f64,

    // f32 positions for GPU upload to the renderer — half the bandwidth of f64
    pub pos_x: Vec<f32>,
    pub pos_y: Vec<f32>,
    pub pos_z: Vec<f32>,

    // Per-atom PSL flexibility score (computed post-hoc; 0.0 until PSL runs)
    pub flexibility: Vec<f32>,

    pub potential_energy: f64,
    pub kinetic_energy: f64,
}

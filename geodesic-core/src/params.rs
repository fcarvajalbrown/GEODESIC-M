/// Immutable run configuration, created once from config.toml and shared
/// via Arc<SimParams> across threads.
#[derive(Debug, Clone)]
pub struct SimParams {
    pub n_atoms: usize,
    pub n_steps: u64,
    pub dt: f64,           // timestep (ps); typically 0.004 ps = 4 fs with Geodesic BAOAB
    pub box_size: [f64; 3], // simulation box (Å); cubic assumed in v1

    pub r_cutoff: f64,
    pub r_skin: f64,
    pub r_switch: f64,

    pub max_constr_iter: u32, // I_max for Geodesic A-step
    pub constr_tol: f64,      // convergence threshold for |λ_i|

    pub frame_interval: u32, // steps between trajectory snapshots
    pub n_threads: usize,    // T; recorded in log for reproducibility

    // E; defines the Jacobi metric g^J_ij = 2(E - V) m_i δ_ij at every step
    pub total_energy: f64,
}

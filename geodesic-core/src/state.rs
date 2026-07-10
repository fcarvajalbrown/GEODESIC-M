/// Mutable integration state, SoA layout, passed through every BAB step.
pub struct SimState {
    // Positions (Å)
    pub pos_x: Vec<f64>,
    pub pos_y: Vec<f64>,
    pub pos_z: Vec<f64>,

    // Velocities (Å/ps)
    pub vel_x: Vec<f64>,
    pub vel_y: Vec<f64>,
    pub vel_z: Vec<f64>,

    // Net forces (kcal/mol·Å) — overwritten each step
    pub force_x: Vec<f64>,
    pub force_y: Vec<f64>,
    pub force_z: Vec<f64>,

    pub potential_energy: f64,
    pub kinetic_energy: f64,

    pub step: u64,
}

impl SimState {
    /// All position/velocity/force arrays zero-initialized to length `n_atoms`.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            pos_x: vec![0.0; n_atoms],
            pos_y: vec![0.0; n_atoms],
            pos_z: vec![0.0; n_atoms],
            vel_x: vec![0.0; n_atoms],
            vel_y: vec![0.0; n_atoms],
            vel_z: vec![0.0; n_atoms],
            force_x: vec![0.0; n_atoms],
            force_y: vec![0.0; n_atoms],
            force_z: vec![0.0; n_atoms],
            potential_energy: 0.0,
            kinetic_energy: 0.0,
            step: 0,
        }
    }
}

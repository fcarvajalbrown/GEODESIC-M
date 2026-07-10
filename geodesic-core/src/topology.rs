/// Force field connectivity, SoA throughout — the bonded force loop reads
/// each array contiguously in sequence.
pub struct BondedTopology {
    // Bond stretching — one entry per bond
    pub bond_i: Vec<u32>,
    pub bond_j: Vec<u32>,
    pub bond_k: Vec<f64>,  // force constant (kcal/mol·Å²)
    pub bond_r0: Vec<f64>, // equilibrium length (Å)

    // Angle bending — one entry per angle
    pub angle_i: Vec<u32>,
    pub angle_j: Vec<u32>, // central atom
    pub angle_k: Vec<u32>,
    pub angle_kth: Vec<f64>, // force constant (kcal/mol·rad²)
    pub angle_th0: Vec<f64>, // equilibrium angle (rad)

    // Dihedral torsion — one entry per dihedral
    pub dihed_i: Vec<u32>,
    pub dihed_j: Vec<u32>,
    pub dihed_k: Vec<u32>,
    pub dihed_l: Vec<u32>,
    pub dihed_kphi: Vec<f64>,  // barrier height (kcal/mol)
    pub dihed_n: Vec<u32>,     // multiplicity
    pub dihed_delta: Vec<f64>, // phase (rad)

    // Holonomic constraints for the Geodesic BAOAB A-step
    // (typically bond lengths involving hydrogen)
    pub constr_i: Vec<u32>,
    pub constr_j: Vec<u32>,
    pub constr_dsq: Vec<f64>, // target |r_i - r_j|² (Å²)

    // Non-bonded exclusions: 1-2, 1-3, and 1-4 pairs, filtered out of the
    // Verlet list at build time so bonded neighbors don't also get LJ forces
    pub excl_i: Vec<u32>,
    pub excl_j: Vec<u32>,
}

/// Verlet pair list, rebuilt when any atom displaces more than
/// (r_skin - r_cutoff) / 2 since the last build.
pub struct NeighborList {
    // Flat list of all pairs (i, j) with i < j within r_skin
    pub pair_i: Vec<u32>,
    pub pair_j: Vec<u32>,

    // Positions at last rebuild — used for the displacement check
    pub ref_x: Vec<f64>,
    pub ref_y: Vec<f64>,
    pub ref_z: Vec<f64>,

    pub r_cutoff: f64, // force goes to zero beyond this
    pub r_skin: f64,   // list radius; r_skin > r_cutoff
    pub r_switch: f64, // switching function onset
}

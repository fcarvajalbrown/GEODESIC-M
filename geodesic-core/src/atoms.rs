/// Static per-atom properties, read-only after initialization.
pub struct AtomData {
    // LJ parameters — SoA, used in non-bonded inner loop
    pub epsilon: Vec<f64>, // kcal/mol
    pub sigma: Vec<f64>,   // Å
    pub mass: Vec<f64>,    // amu
    pub charge: Vec<f64>,  // elementary charge (reserved; not used in v1)

    // Metadata — AoS, used only for I/O and GUI coloring
    pub meta: Vec<AtomMeta>,
}

pub struct AtomMeta {
    pub element: Element,
    pub residue_id: u32,
    pub atom_name: [u8; 4], // PDB atom name field
    pub chain_id: u8,
}

/// Extended as the prmtop parser (Phase 2) encounters new element types.
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Element {
    H,
    C,
    N,
    O,
    S,
    Unknown,
}

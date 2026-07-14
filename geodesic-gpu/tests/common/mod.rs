use geodesic_core::{AtomData, BondedTopology, NeighborList, SimParams, SimState};
use geodesic_engine::force::nonbonded;
use geodesic_engine::neighbor;

pub const FIX_DIR: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../geodesic-engine/tests/fixtures");

pub fn load(name: &str) -> (SimState, AtomData, BondedTopology) {
    let prmtop = std::fs::read_to_string(format!("{FIX_DIR}/{name}.prmtop")).unwrap();
    let (atoms, topology) = geodesic_io::prmtop::parse(&prmtop).unwrap();
    let inpcrd = std::fs::read_to_string(format!("{FIX_DIR}/{name}.inpcrd")).unwrap();
    let state = geodesic_io::inpcrd::parse(&inpcrd, atoms.mass.len(), false).unwrap();
    (state, atoms, topology)
}

pub fn params(n_atoms: usize) -> SimParams {
    SimParams {
        n_atoms,
        n_steps: 0,
        dt: 0.004,
        box_size: [100.0, 100.0, 100.0],
        periodic: true,
        r_cutoff: 12.0,
        r_skin: 14.0,
        r_switch: 10.0,
        max_constr_iter: 100,
        constr_tol: 1e-10,
        frame_interval: 1,
        n_threads: 1,
        total_energy: 0.0,
    }
}

pub fn clone_positions(src: &SimState) -> SimState {
    let mut s = SimState::new(src.pos_x.len());
    s.pos_x = src.pos_x.clone();
    s.pos_y = src.pos_y.clone();
    s.pos_z = src.pos_z.clone();
    s
}

/// Wrap a fresh state via the neighbor build, then compute the CPU non-bonded
/// reference forces on that same wrapped state.
pub fn cpu_nonbonded_reference(
    state: &SimState,
    atoms: &AtomData,
    topology: &BondedTopology,
    p: &SimParams,
) -> (SimState, NeighborList, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = state.pos_x.len();
    let mut s = clone_positions(state);
    let list = neighbor::build(&mut s, p, topology);
    let (mut fx, mut fy, mut fz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    nonbonded::compute_pair_forces(
        &s, atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, p.box_size,
        &mut fx, &mut fy, &mut fz,
    );
    (s, list, fx, fy, fz)
}

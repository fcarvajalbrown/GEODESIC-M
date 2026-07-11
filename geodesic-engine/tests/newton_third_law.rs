use geodesic_core::{AtomData, BondedTopology, SimParams, SimState};
use geodesic_engine::force::{bonded, nonbonded};
use geodesic_engine::neighbor;

/// SAD.md §13.3: for every pair the force engine evaluates, F_ij = -F_ji, so
/// the net force over the whole system must vanish. Runs on all fixture
/// systems.
fn load(name: &str) -> (SimState, AtomData, BondedTopology) {
    let prmtop_text = std::fs::read_to_string(format!("tests/fixtures/{name}.prmtop")).unwrap();
    let (atoms, topology) = geodesic_io::prmtop::parse(&prmtop_text).unwrap();
    let inpcrd_text = std::fs::read_to_string(format!("tests/fixtures/{name}.inpcrd")).unwrap();
    let state = geodesic_io::inpcrd::parse(&inpcrd_text, atoms.mass.len(), false).unwrap();
    (state, atoms, topology)
}

fn params(n_atoms: usize) -> SimParams {
    SimParams {
        n_atoms,
        n_steps: 0,
        dt: 0.004,
        box_size: [100.0, 100.0, 100.0],
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

fn check_net_force_vanishes(name: &str) {
    let (state, atoms, topology) = load(name);
    let n = state.pos_x.len();
    let p = params(n);

    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut fz = vec![0.0; n];

    bonded::compute_bond_forces(&state, &topology, &mut fx, &mut fy, &mut fz);
    bonded::compute_angle_forces(&state, &topology, &mut fx, &mut fy, &mut fz);
    bonded::compute_dihedral_forces(&state, &topology, &mut fx, &mut fy, &mut fz);

    let mut wrap_state = SimState::new(n);
    wrap_state.pos_x = state.pos_x.clone();
    wrap_state.pos_y = state.pos_y.clone();
    wrap_state.pos_z = state.pos_z.clone();
    let list = neighbor::build(&mut wrap_state, &p, &topology);
    nonbonded::compute_pair_forces(
        &state, &atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, p.box_size,
        &mut fx, &mut fy, &mut fz,
    );

    let n_eps = n as f64 * 1e-10;
    let (sx, sy, sz) = (fx.iter().sum::<f64>(), fy.iter().sum::<f64>(), fz.iter().sum::<f64>());
    assert!(sx.abs() < n_eps, "{name}: net force x = {sx}");
    assert!(sy.abs() < n_eps, "{name}: net force y = {sy}");
    assert!(sz.abs() < n_eps, "{name}: net force z = {sz}");
}

#[test]
fn lj_pair_net_force_vanishes() {
    check_net_force_vanishes("lj_pair");
}

#[test]
fn harmonic_dimer_net_force_vanishes() {
    check_net_force_vanishes("harmonic_dimer");
}

#[test]
fn water_box_4_net_force_vanishes() {
    check_net_force_vanishes("water_box_4");
}

#[test]
fn ala_dipeptide_net_force_vanishes() {
    check_net_force_vanishes("ala_dipeptide");
}

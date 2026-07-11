use geodesic_core::{AtomData, BondedTopology, SimParams, SimState};
use geodesic_engine::force::{bonded, nonbonded};
use geodesic_engine::neighbor;

/// SAD.md §13.1/§13.2: the finite-difference gradient check run against the
/// committed fixture files (as opposed to the hand-built minimal systems in
/// bonded_gradient.rs/nonbonded_gradient.rs, which isolate one force term at
/// a time). This proves the fixtures parse into a physically consistent
/// system and exercises bond+angle+LJ together, not just each in isolation.
fn load(name: &str) -> (SimState, AtomData, BondedTopology) {
    let prmtop_text =
        std::fs::read_to_string(format!("tests/fixtures/{name}.prmtop")).unwrap();
    let (atoms, topology) = geodesic_io::prmtop::parse(&prmtop_text).unwrap();
    let inpcrd_text =
        std::fs::read_to_string(format!("tests/fixtures/{name}.inpcrd")).unwrap();
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

/// Total potential (bond + angle + dihedral + LJ) and its analytic force,
/// with no constraint promotion — this checks the force-field terms
/// themselves, independent of the constraint solver (SAD.md §13.2 vs §13.6
/// are deliberately separate tests).
fn total_potential_and_force(
    state: &SimState,
    atoms: &AtomData,
    topology: &BondedTopology,
    p: &SimParams,
) -> (f64, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = state.pos_x.len();
    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut fz = vec![0.0; n];
    let mut v = 0.0;

    v += bonded::compute_bond_forces(state, topology, &mut fx, &mut fy, &mut fz);
    v += bonded::compute_angle_forces(state, topology, &mut fx, &mut fy, &mut fz);
    v += bonded::compute_dihedral_forces(state, topology, &mut fx, &mut fy, &mut fz);

    let mut wrap_state = SimState::new(n);
    wrap_state.pos_x = state.pos_x.clone();
    wrap_state.pos_y = state.pos_y.clone();
    wrap_state.pos_z = state.pos_z.clone();
    let list = neighbor::build(&mut wrap_state, p, topology);
    v += nonbonded::compute_pair_forces(
        state, atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, p.box_size,
        &mut fx, &mut fy, &mut fz,
    );

    (v, fx, fy, fz)
}

fn check_gradient_matches_fixture(name: &str) {
    let (state, atoms, topology) = load(name);
    let p = params(state.pos_x.len());
    let n = state.pos_x.len();
    let eps = 1e-6;

    let (_, fx, fy, fz) = total_potential_and_force(&state, &atoms, &topology, &p);

    for i in 0..n {
        for (axis, analytic) in [(0, fx[i]), (1, fy[i]), (2, fz[i])] {
            let mut plus = SimState::new(n);
            plus.pos_x = state.pos_x.clone();
            plus.pos_y = state.pos_y.clone();
            plus.pos_z = state.pos_z.clone();
            let mut minus = SimState::new(n);
            minus.pos_x = state.pos_x.clone();
            minus.pos_y = state.pos_y.clone();
            minus.pos_z = state.pos_z.clone();

            match axis {
                0 => {
                    plus.pos_x[i] += eps;
                    minus.pos_x[i] -= eps;
                }
                1 => {
                    plus.pos_y[i] += eps;
                    minus.pos_y[i] -= eps;
                }
                _ => {
                    plus.pos_z[i] += eps;
                    minus.pos_z[i] -= eps;
                }
            }

            let (v_plus, ..) = total_potential_and_force(&plus, &atoms, &topology, &p);
            let (v_minus, ..) = total_potential_and_force(&minus, &atoms, &topology, &p);
            let numeric = -(v_plus - v_minus) / (2.0 * eps);
            let rel_err = (analytic - numeric).abs() / (numeric.abs() + 1.0);
            assert!(
                rel_err < 1e-4,
                "{name}: atom {i} axis {axis}: analytic={analytic}, numeric={numeric}, rel_err={rel_err}"
            );
        }
    }
}

#[test]
fn lj_pair_gradient_matches_finite_difference() {
    check_gradient_matches_fixture("lj_pair");
}

#[test]
fn harmonic_dimer_gradient_matches_finite_difference() {
    check_gradient_matches_fixture("harmonic_dimer");
}

#[test]
fn water_box_4_gradient_matches_finite_difference() {
    check_gradient_matches_fixture("water_box_4");
}

#[test]
fn ala_dipeptide_gradient_matches_finite_difference() {
    check_gradient_matches_fixture("ala_dipeptide");
}

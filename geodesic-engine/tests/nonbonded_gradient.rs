use geodesic_core::{AtomData, AtomMeta, Element, NeighborList, SimState};
use geodesic_engine::force::nonbonded::compute_pair_forces;

fn two_atom_system(r: f64) -> (SimState, AtomData) {
    let mut state = SimState::new(2);
    state.pos_x = vec![0.0, r];

    let meta = AtomMeta {
        element: Element::Unknown,
        residue_id: 0,
        atom_name: *b"AT  ",
        chain_id: 0,
    };
    let atoms = AtomData {
        epsilon: vec![1.0, 1.0],
        sigma: vec![2.0, 2.0],
        mass: vec![1.0, 1.0],
        charge: vec![0.0, 0.0],
        meta: vec![meta, meta],
    };
    (state, atoms)
}

fn neighbor_list(r_cutoff: f64, r_switch: f64) -> NeighborList {
    NeighborList {
        pair_i: vec![0],
        pair_j: vec![1],
        ref_x: vec![0.0, 0.0],
        ref_y: vec![0.0, 0.0],
        ref_z: vec![0.0, 0.0],
        r_cutoff,
        r_skin: r_cutoff + 1.0,
        r_switch,
    }
}

fn potential_at(r: f64, r_cutoff: f64, r_switch: f64) -> f64 {
    let (state, atoms) = two_atom_system(r);
    let list = neighbor_list(r_cutoff, r_switch);
    let box_size = [1000.0, 1000.0, 1000.0];
    let mut fx = vec![0.0; 2];
    let mut fy = vec![0.0; 2];
    let mut fz = vec![0.0; 2];
    compute_pair_forces(
        &state, &atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, box_size,
        &mut fx, &mut fy, &mut fz,
    )
}

fn analytic_force_x1(r: f64, r_cutoff: f64, r_switch: f64) -> f64 {
    let (state, atoms) = two_atom_system(r);
    let list = neighbor_list(r_cutoff, r_switch);
    let box_size = [1000.0, 1000.0, 1000.0];
    let mut fx = vec![0.0; 2];
    let mut fy = vec![0.0; 2];
    let mut fz = vec![0.0; 2];
    compute_pair_forces(
        &state, &atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, box_size,
        &mut fx, &mut fy, &mut fz,
    );
    fx[1]
}

fn check_gradient(r: f64, r_cutoff: f64, r_switch: f64) {
    let eps = 1e-6;
    let v_plus = potential_at(r + eps, r_cutoff, r_switch);
    let v_minus = potential_at(r - eps, r_cutoff, r_switch);
    let numeric = -(v_plus - v_minus) / (2.0 * eps);
    let analytic = analytic_force_x1(r, r_cutoff, r_switch);
    let rel_err = (analytic - numeric).abs() / (numeric.abs() + 1.0);
    assert!(
        rel_err < 1e-5,
        "r={r}: analytic={analytic}, numeric={numeric}, rel_err={rel_err}"
    );
}

#[test]
fn gradient_matches_finite_difference_below_switch() {
    // r well inside r_switch=10 -> unswitched LJ regime
    check_gradient(2.5, 12.0, 10.0);
    check_gradient(3.0, 12.0, 10.0);
}

#[test]
fn gradient_matches_finite_difference_in_switch_region() {
    check_gradient(10.5, 12.0, 10.0);
    check_gradient(11.5, 12.0, 10.0);
}

#[test]
fn switching_function_is_continuous_at_boundaries() {
    let r_cutoff = 12.0;
    let r_switch = 10.0;
    // value continuity: just inside vs just past r_switch, and near r_cutoff vs zero
    let just_below = potential_at(r_switch - 1e-7, r_cutoff, r_switch);
    let just_above = potential_at(r_switch + 1e-7, r_cutoff, r_switch);
    assert!((just_below - just_above).abs() < 1e-4);

    let near_cutoff = potential_at(r_cutoff - 1e-6, r_cutoff, r_switch);
    assert!(near_cutoff.abs() < 1e-3, "potential should vanish approaching r_cutoff, got {near_cutoff}");
}

#[test]
fn repulsive_force_pushes_atoms_apart() {
    // r < sigma -> strongly repulsive regime
    let f = analytic_force_x1(1.5, 12.0, 10.0);
    assert!(f > 0.0, "atom 1 (at larger x) should be pushed further away: f={f}");
}

#[test]
fn newtons_third_law_holds() {
    let (state, atoms) = two_atom_system(3.0);
    let list = neighbor_list(12.0, 10.0);
    let box_size = [1000.0, 1000.0, 1000.0];
    let mut fx = vec![0.0; 2];
    let mut fy = vec![0.0; 2];
    let mut fz = vec![0.0; 2];
    compute_pair_forces(
        &state, &atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, box_size,
        &mut fx, &mut fy, &mut fz,
    );
    assert!((fx[0] + fx[1]).abs() < 1e-12);
}

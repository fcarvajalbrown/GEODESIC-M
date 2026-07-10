use geodesic_core::{BondedTopology, SimParams, SimState};
use geodesic_engine::neighbor;

fn empty_bonded() -> BondedTopology {
    BondedTopology {
        bond_i: vec![],
        bond_j: vec![],
        bond_k: vec![],
        bond_r0: vec![],
        angle_i: vec![],
        angle_j: vec![],
        angle_k: vec![],
        angle_kth: vec![],
        angle_th0: vec![],
        dihed_i: vec![],
        dihed_j: vec![],
        dihed_k: vec![],
        dihed_l: vec![],
        dihed_kphi: vec![],
        dihed_n: vec![],
        dihed_delta: vec![],
        constr_i: vec![],
        constr_j: vec![],
        constr_dsq: vec![],
        excl_i: vec![],
        excl_j: vec![],
    }
}

fn params(box_size: [f64; 3]) -> SimParams {
    SimParams {
        n_atoms: 0,
        n_steps: 0,
        dt: 0.004,
        box_size,
        r_cutoff: 5.0,
        r_skin: 6.0,
        r_switch: 4.0,
        max_constr_iter: 100,
        constr_tol: 1e-6,
        frame_interval: 1,
        n_threads: 1,
        total_energy: 0.0,
    }
}

#[test]
fn finds_pairs_within_skin_and_excludes_far_pairs() {
    let mut state = SimState::new(3);
    state.pos_x = vec![0.0, 3.0, 20.0];
    let p = params([1000.0, 1000.0, 1000.0]);
    let list = neighbor::build(&mut state, &p, &empty_bonded());

    // (0,1) within r_skin=6, (0,2) and (1,2) far beyond it
    assert_eq!(list.pair_i, vec![0]);
    assert_eq!(list.pair_j, vec![1]);
}

#[test]
fn bonded_exclusions_are_filtered_out() {
    let mut state = SimState::new(2);
    state.pos_x = vec![0.0, 1.0];
    let p = params([1000.0, 1000.0, 1000.0]);
    let mut bonded = empty_bonded();
    bonded.excl_i = vec![0];
    bonded.excl_j = vec![1];

    let list = neighbor::build(&mut state, &p, &bonded);
    assert!(list.pair_i.is_empty(), "excluded pair should not appear in the neighbor list");
}

#[test]
fn minimum_image_finds_pair_across_periodic_boundary() {
    let mut state = SimState::new(2);
    // box=10; atoms at x=0.5 and x=9.5 are 1.0 apart through the boundary,
    // 9.0 apart the "raw" way -> only detected as neighbors via min-image
    state.pos_x = vec![0.5, 9.5];
    let p = params([10.0, 10.0, 10.0]);
    let list = neighbor::build(&mut state, &p, &empty_bonded());

    assert_eq!(list.pair_i, vec![0]);
    assert_eq!(list.pair_j, vec![1]);
}

#[test]
fn positions_are_wrapped_into_box() {
    let mut state = SimState::new(1);
    state.pos_x = vec![-1.0];
    state.pos_y = vec![11.0];
    let p = params([10.0, 10.0, 10.0]);
    let _list = neighbor::build(&mut state, &p, &empty_bonded());

    assert!((state.pos_x[0] - 9.0).abs() < 1e-12);
    assert!((state.pos_y[0] - 1.0).abs() < 1e-12);
}

#[test]
fn rebuild_triggered_past_half_skin_displacement() {
    let mut state = SimState::new(1);
    let p = params([1000.0, 1000.0, 1000.0]);
    let list = neighbor::build(&mut state, &p, &empty_bonded());

    // half-skin = (6.0 - 5.0) / 2 = 0.5
    state.pos_x[0] = 0.4;
    assert!(!neighbor::needs_rebuild(&state, &list));
    state.pos_x[0] = 0.6;
    assert!(neighbor::needs_rebuild(&state, &list));
}

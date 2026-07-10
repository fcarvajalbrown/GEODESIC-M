use geodesic_core::{AtomData, AtomMeta, BondedTopology, ConvergenceError, Element};
use geodesic_engine::constraint::solve;

fn empty_topology() -> BondedTopology {
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

/// Two atoms (H, C) held at an equilibrium separation of 1.0 A by a single
/// holonomic constraint — SAD.md §13.6's "harmonic dimer".
fn harmonic_dimer() -> (BondedTopology, AtomData) {
    let mut topo = empty_topology();
    topo.constr_i = vec![0];
    topo.constr_j = vec![1];
    topo.constr_dsq = vec![1.0]; // 1.0 A target separation, squared

    let atoms = AtomData {
        epsilon: vec![0.0, 0.0],
        sigma: vec![0.0, 0.0],
        mass: vec![1.008, 12.0],
        charge: vec![0.0, 0.0],
        meta: vec![
            AtomMeta { element: Element::H, residue_id: 0, atom_name: *b"H1  ", chain_id: 0 },
            AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C1  ", chain_id: 0 },
        ],
    };

    (topo, atoms)
}

#[test]
fn converges_within_max_iter_for_valid_input() {
    let (topo, atoms) = harmonic_dimer();
    let ref_x = vec![0.0, 1.0];
    let ref_y = vec![0.0, 0.0];
    let ref_z = vec![0.0, 0.0];
    // trial position after a drift step: bond stretched 5% off target
    let mut pos_x = vec![0.0, 1.05];
    let mut pos_y = vec![0.0, 0.0];
    let mut pos_z = vec![0.0, 0.0];

    let result = solve(
        &topo, &atoms, &ref_x, &ref_y, &ref_z, &mut pos_x, &mut pos_y, &mut pos_z, 50, 1e-10, 0,
    );

    assert!(result.is_ok(), "expected convergence, got {:?}", result);
    let iters = result.unwrap();
    assert!(iters <= 50);
    assert!(iters > 0);
}

#[test]
fn manifold_adherence_after_solve() {
    let (topo, atoms) = harmonic_dimer();
    let ref_x = vec![0.0, 1.0];
    let ref_y = vec![0.0, 0.0];
    let ref_z = vec![0.0, 0.0];
    let mut pos_x = vec![0.0, 1.05];
    let mut pos_y = vec![0.0, 0.0];
    let mut pos_z = vec![0.0, 0.0];

    solve(&topo, &atoms, &ref_x, &ref_y, &ref_z, &mut pos_x, &mut pos_y, &mut pos_z, 50, 1e-10, 0)
        .expect("solver should converge");

    let dx = pos_x[0] - pos_x[1];
    let dy = pos_y[0] - pos_y[1];
    let dz = pos_z[0] - pos_z[1];
    let dsq = dx * dx + dy * dy + dz * dz;
    let d0sq = topo.constr_dsq[0];

    assert!(
        (dsq - d0sq).abs() < 1e-8,
        "constrained bond length^2 {dsq} deviates from target {d0sq} by more than 1e-8"
    );
}

#[test]
fn max_iter_one_returns_convergence_error_not_wrong_result() {
    let (topo, atoms) = harmonic_dimer();
    let ref_x = vec![0.0, 1.0];
    let ref_y = vec![0.0, 0.0];
    let ref_z = vec![0.0, 0.0];
    let mut pos_x = vec![0.0, 1.05];
    let mut pos_y = vec![0.0, 0.0];
    let mut pos_z = vec![0.0, 0.0];

    let result = solve(
        &topo, &atoms, &ref_x, &ref_y, &ref_z, &mut pos_x, &mut pos_y, &mut pos_z, 1, 1e-10, 7,
    );

    match result {
        Err(ConvergenceError::ConstraintSolver { step, max_iter, atom_i, atom_j, .. }) => {
            assert_eq!(step, 7);
            assert_eq!(max_iter, 1);
            assert_eq!((atom_i, atom_j), (0, 1));
        }
        other => panic!("expected ConvergenceError::ConstraintSolver, got {:?}", other),
    }
}

#[test]
fn no_constraints_converges_immediately() {
    let topo = empty_topology();
    let atoms = AtomData {
        epsilon: vec![0.0],
        sigma: vec![0.0],
        mass: vec![1.0],
        charge: vec![0.0],
        meta: vec![AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C1  ", chain_id: 0 }],
    };
    let ref_pos = vec![0.0];
    let mut pos_x = vec![0.0];
    let mut pos_y = vec![0.0];
    let mut pos_z = vec![0.0];

    let result = solve(
        &topo, &atoms, &ref_pos, &ref_pos, &ref_pos, &mut pos_x, &mut pos_y, &mut pos_z, 50,
        1e-10, 0,
    );
    assert_eq!(result.unwrap(), 0);
}

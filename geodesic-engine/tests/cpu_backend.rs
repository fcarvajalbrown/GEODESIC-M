use geodesic_core::{AtomData, AtomMeta, BondedTopology, ComputeBackend, Element, SimParams, SimState};
use geodesic_engine::cpu_backend::CpuBackend;
use geodesic_engine::force::nonbonded;
use geodesic_engine::neighbor;

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

fn params(n_atoms: usize, n_threads: usize) -> SimParams {
    SimParams {
        n_atoms,
        n_steps: 0,
        dt: 0.004,
        box_size: [100.0, 100.0, 100.0],
        r_cutoff: 8.0,
        r_skin: 10.0,
        r_switch: 6.0,
        max_constr_iter: 100,
        constr_tol: 1e-10,
        frame_interval: 1,
        n_threads,
        total_energy: 0.0,
    }
}

/// 12 atoms on a small perturbed grid — close enough that many pairs fall
/// within r_cutoff=8.0, so a 4-thread strip decomposition actually spans
/// multiple non-trivial chunks.
fn lj_cluster() -> (SimState, AtomData) {
    let mut state = SimState::new(12);
    let coords = [
        (0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (3.0, 3.0, 0.0),
        (0.0, 0.0, 3.0), (3.0, 0.0, 3.0), (0.0, 3.0, 3.0), (3.0, 3.0, 3.0),
        (1.5, 1.5, 1.5), (4.5, 1.5, 1.5), (1.5, 4.5, 1.5), (1.5, 1.5, 4.5),
    ];
    for (n, &(x, y, z)) in coords.iter().enumerate() {
        state.pos_x[n] = x;
        state.pos_y[n] = y;
        state.pos_z[n] = z;
    }

    let meta = AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C1  ", chain_id: 0 };
    let atoms = AtomData {
        epsilon: vec![0.2; 12],
        sigma: vec![1.5; 12],
        mass: vec![12.0; 12],
        charge: vec![0.0; 12],
        meta: vec![meta; 12],
    };
    (state, atoms)
}

#[test]
fn compute_forces_matches_direct_nonbonded_call_single_thread() {
    let (mut state, _) = lj_cluster();
    let p = params(12, 1);

    let mut backend = CpuBackend::new(lj_cluster().1, empty_topology(), &p);
    backend.build_neighbor_list(&mut state, &p);
    let via_backend = backend.compute_forces(&state).clone();

    let (mut state_direct, atoms_direct) = lj_cluster();
    let list = neighbor::build(&mut state_direct, &p, &empty_topology());
    let mut fx = vec![0.0; 12];
    let mut fy = vec![0.0; 12];
    let mut fz = vec![0.0; 12];
    nonbonded::compute_pair_forces(
        &state_direct, &atoms_direct, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch,
        p.box_size, &mut fx, &mut fy, &mut fz,
    );

    for i in 0..12 {
        assert!((via_backend.fx[i] - fx[i]).abs() < 1e-12, "fx[{i}] mismatch");
        assert!((via_backend.fy[i] - fy[i]).abs() < 1e-12, "fy[{i}] mismatch");
        assert!((via_backend.fz[i] - fz[i]).abs() < 1e-12, "fz[{i}] mismatch");
    }
}

/// SAD.md §7.2's determinism claim is that T (thread count) is a
/// reproducibility *parameter* — the same T must reproduce bit-for-bit run
/// after run, since it fixes the static strip partition. It does NOT claim
/// results are bit-identical *across different* T: floating-point addition
/// isn't associative, so summing one atom's pair contributions via several
/// per-thread partial sums (4 threads) vs. one running sum (1 thread) can
/// legitimately differ in the last bit even though both used the exact same
/// static, non-work-stealing decomposition.
#[test]
fn nonbonded_reduction_is_repeatable_for_a_fixed_thread_count() {
    for n_threads in [1, 4] {
        let (mut state_a, atoms_a) = lj_cluster();
        let p = params(12, n_threads);
        let mut backend_a = CpuBackend::new(atoms_a, empty_topology(), &p);
        backend_a.build_neighbor_list(&mut state_a, &p);
        let run_a = backend_a.compute_forces(&state_a).clone();

        let (mut state_b, atoms_b) = lj_cluster();
        let mut backend_b = CpuBackend::new(atoms_b, empty_topology(), &p);
        backend_b.build_neighbor_list(&mut state_b, &p);
        let run_b = backend_b.compute_forces(&state_b).clone();

        assert_eq!(backend_a.n_threads(), n_threads);
        for i in 0..12 {
            assert_eq!(run_a.fx[i].to_bits(), run_b.fx[i].to_bits(), "T={n_threads}: fx[{i}] not repeatable");
            assert_eq!(run_a.fy[i].to_bits(), run_b.fy[i].to_bits(), "T={n_threads}: fy[{i}] not repeatable");
            assert_eq!(run_a.fz[i].to_bits(), run_b.fz[i].to_bits(), "T={n_threads}: fz[{i}] not repeatable");
        }
    }
}

/// Cross-check that different thread counts still agree to numerical
/// tolerance (not bit-for-bit — see the repeatability test above for why).
#[test]
fn nonbonded_reduction_agrees_across_thread_counts_within_tolerance() {
    let (mut state1, atoms1) = lj_cluster();
    let p1 = params(12, 1);
    let mut backend1 = CpuBackend::new(atoms1, empty_topology(), &p1);
    backend1.build_neighbor_list(&mut state1, &p1);
    let f1 = backend1.compute_forces(&state1).clone();

    let (mut state4, atoms4) = lj_cluster();
    let p4 = params(12, 4);
    let mut backend4 = CpuBackend::new(atoms4, empty_topology(), &p4);
    backend4.build_neighbor_list(&mut state4, &p4);
    let f4 = backend4.compute_forces(&state4).clone();

    for i in 0..12 {
        assert!((f1.fx[i] - f4.fx[i]).abs() < 1e-10, "fx[{i}] mismatch beyond floating-point tolerance");
        assert!((f1.fy[i] - f4.fy[i]).abs() < 1e-10, "fy[{i}] mismatch beyond floating-point tolerance");
        assert!((f1.fz[i] - f4.fz[i]).abs() < 1e-10, "fz[{i}] mismatch beyond floating-point tolerance");
    }
}

#[test]
fn reduce_forces_matches_last_compute_forces() {
    let (mut state, atoms) = lj_cluster();
    let topo = empty_topology();
    let p = params(12, 4);
    let mut backend = CpuBackend::new(atoms, topo, &p);
    backend.build_neighbor_list(&mut state, &p);
    let computed = backend.compute_forces(&state).clone();
    let reduced = backend.reduce_forces();

    for i in 0..12 {
        assert_eq!(computed.fx[i].to_bits(), reduced.fx[i].to_bits());
        assert_eq!(computed.fy[i].to_bits(), reduced.fy[i].to_bits());
        assert_eq!(computed.fz[i].to_bits(), reduced.fz[i].to_bits());
    }
}

#[test]
fn geodesic_drift_respects_bond_constraint() {
    let mut topo = empty_topology();
    topo.constr_i = vec![0];
    topo.constr_j = vec![1];
    topo.constr_dsq = vec![1.0];

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

    let mut state = SimState::new(2);
    state.pos_x = vec![0.0, 1.0];
    state.vel_x = vec![-0.2, 0.2];

    let p = params(2, 1);
    let mut backend = CpuBackend::new(atoms, topo, &p);
    backend.geodesic_drift(&mut state, 0.1).expect("should converge");

    let dx = state.pos_x[0] - state.pos_x[1];
    assert!((dx * dx - 1.0).abs() < 1e-9, "bond^2 drifted: {}", dx * dx);
}

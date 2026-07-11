use geodesic_core::{AtomData, AtomMeta, BondedTopology, Element, SimState};
use geodesic_engine::constraint::constrain_velocities;
use geodesic_engine::integrator::{drift_and_constrain, half_kick, FORCE_TO_ACCEL_ANG_PER_PS2};

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

fn harmonic_dimer_topology() -> (BondedTopology, AtomData) {
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
    (topo, atoms)
}

#[test]
fn half_kick_applies_dt_half_over_mass_times_force() {
    let mut state = SimState::new(2);
    state.force_x = vec![10.0, -4.0];
    state.force_y = vec![0.0, 2.0];
    state.force_z = vec![1.0, 0.0];
    let atoms = AtomData {
        epsilon: vec![0.0, 0.0],
        sigma: vec![0.0, 0.0],
        mass: vec![2.0, 4.0],
        charge: vec![0.0, 0.0],
        meta: vec![
            AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C1  ", chain_id: 0 },
            AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C2  ", chain_id: 0 },
        ],
    };

    half_kick(&mut state, &atoms, 0.5);

    let c = FORCE_TO_ACCEL_ANG_PER_PS2;
    assert!((state.vel_x[0] - 10.0 * 0.5 / 2.0 * c).abs() < 1e-11);
    assert!((state.vel_z[0] - 1.0 * 0.5 / 2.0 * c).abs() < 1e-11);
    assert!((state.vel_x[1] - (-4.0) * 0.5 / 4.0 * c).abs() < 1e-11);
    assert!((state.vel_y[1] - 2.0 * 0.5 / 4.0 * c).abs() < 1e-11);
}

#[test]
fn drift_without_constraints_is_plain_free_drift() {
    let topo = empty_topology();
    let atoms = AtomData {
        epsilon: vec![0.0],
        sigma: vec![0.0],
        mass: vec![1.0],
        charge: vec![0.0],
        meta: vec![AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C1  ", chain_id: 0 }],
    };
    let mut state = SimState::new(1);
    state.pos_x = vec![1.0];
    state.vel_x = vec![2.0];

    drift_and_constrain(&mut state, &topo, &atoms, 0.5, 10, 1e-10).unwrap();

    assert!((state.pos_x[0] - 2.0).abs() < 1e-14); // 1.0 + 0.5*2.0
    assert!((state.vel_x[0] - 2.0).abs() < 1e-14); // unchanged
}

#[test]
fn drift_and_constrain_keeps_bond_on_manifold_and_resyncs_velocity() {
    let (topo, atoms) = harmonic_dimer_topology();
    let mut state = SimState::new(2);
    state.pos_x = vec![0.0, 1.0];
    // velocity that would stretch the bond if left unconstrained
    state.vel_x = vec![-0.2, 0.2];

    drift_and_constrain(&mut state, &topo, &atoms, 0.1, 50, 1e-12).unwrap();

    let dx = state.pos_x[0] - state.pos_x[1];
    assert!((dx * dx - 1.0).abs() < 1e-9, "bond^2 drifted: {}", dx * dx);

    // velocity must equal the actual displacement / dt, not the pre-drift value
    let expected_vx0 = (state.pos_x[0] - 0.0) / 0.1;
    assert!((state.vel_x[0] - expected_vx0).abs() < 1e-9);
}

/// Free rigid rotor: two masses on a rigid rod, zero external force,
/// translating and rotating about the center of mass. The BAB + RATTLE
/// velocity projection loop must conserve kinetic energy and bond length
/// over many steps — parameters and expected tolerances cross-checked
/// against an independent Python reference implementation of the same
/// algorithm (2000 steps, dt=0.004: KE drift ~2e-12, bond error ~2e-16).
#[test]
fn free_rigid_rotor_conserves_energy_and_bond_length() {
    let (topo, atoms) = harmonic_dimer_topology();
    let mut state = SimState::new(2);
    state.pos_x = vec![0.0, 1.0];

    let v_com = [0.1, 0.05, 0.0];
    let omega = [0.0, 0.3, 0.0];
    state.vel_x = vec![v_com[0] + omega[0], v_com[0] - omega[0]];
    state.vel_y = vec![v_com[1] + omega[1], v_com[1] - omega[1]];
    state.vel_z = vec![v_com[2] + omega[2], v_com[2] - omega[2]];
    // no forces: state.force_{x,y,z} stay zero, half_kick is then a no-op

    let ke0: f64 = (0..2)
        .map(|i| 0.5 * atoms.mass[i] * (state.vel_x[i].powi(2) + state.vel_y[i].powi(2) + state.vel_z[i].powi(2)))
        .sum();

    let dt = 0.004;
    for _ in 0..2000 {
        half_kick(&mut state, &atoms, dt / 2.0);
        drift_and_constrain(&mut state, &topo, &atoms, dt, 50, 1e-12).unwrap();
        half_kick(&mut state, &atoms, dt / 2.0);
        constrain_velocities(
            &topo, &atoms, &state.pos_x, &state.pos_y, &state.pos_z, &mut state.vel_x,
            &mut state.vel_y, &mut state.vel_z, 50, 1e-12, state.step,
        )
        .unwrap();
        state.step += 1;
    }

    let ke1: f64 = (0..2)
        .map(|i| 0.5 * atoms.mass[i] * (state.vel_x[i].powi(2) + state.vel_y[i].powi(2) + state.vel_z[i].powi(2)))
        .sum();
    let dx = state.pos_x[0] - state.pos_x[1];
    let dy = state.pos_y[0] - state.pos_y[1];
    let dz = state.pos_z[0] - state.pos_z[1];
    let bond_sq = dx * dx + dy * dy + dz * dz;

    assert!((ke1 - ke0).abs() < 1e-8, "kinetic energy drifted by {}", ke1 - ke0);
    assert!((bond_sq - 1.0).abs() < 1e-8, "bond^2 drifted to {}", bond_sq);
}

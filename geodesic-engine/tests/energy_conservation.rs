use geodesic_core::AtomData;
use geodesic_engine::constraint::{constrain_velocities, promote_hydrogen_bonds};
use geodesic_engine::integrator::{drift_and_constrain, half_kick};

/// SAD.md §13.4: tests the Geodesic BAB integrator itself (not force
/// correctness — that's §13.2/fixture_gradient_check.rs), using the
/// harmonic_dimer fixture. The H-C bond is promoted to a rigid constraint
/// (touches hydrogen), so this exercises the full load -> promote ->
/// BAB+RATTLE pipeline end to end, not just the algorithm on a hand-built
/// system (see integrator.rs's free_rigid_rotor_conserves_energy_and_bond_length
/// for that narrower check).
#[test]
fn harmonic_dimer_conserves_energy_over_100k_steps() {
    let prmtop_text = std::fs::read_to_string("tests/fixtures/harmonic_dimer.prmtop").unwrap();
    let (atoms, mut topology) = geodesic_io::prmtop::parse(&prmtop_text).unwrap();
    let inpcrd_text = std::fs::read_to_string("tests/fixtures/harmonic_dimer.inpcrd").unwrap();
    let mut state = geodesic_io::inpcrd::parse(&inpcrd_text, atoms.mass.len(), false).unwrap();

    promote_hydrogen_bonds(&mut topology, &atoms);
    assert_eq!(topology.constr_i.len(), 1, "H-C bond should have been promoted to a constraint");
    assert_eq!(topology.bond_i.len(), 0, "no harmonic bonds should remain");

    // translation + rotation about the center of mass
    state.vel_x = vec![-0.15, 0.15];
    state.vel_y = vec![0.08, -0.02];
    state.vel_z = vec![0.0, 0.0];
    // no bonded/LJ forces on this fixture once the bond is constrained and
    // epsilon=0, so state.force_{x,y,z} legitimately stay zero throughout

    // A hand-picked initial velocity generally has a component along the
    // bond direction, which a perfectly rigid constraint cannot support --
    // real MD setups always project the initial velocity onto the
    // constraint's tangent space before measuring E(0), exactly like every
    // subsequent step does after its own force kick. Skipping this here
    // isn't a smaller version of the same physics; it measures E(0) against
    // an inconsistent initial condition and reports the one-time projection
    // as spurious "drift" (confirmed independently in Python: without this,
    // step 0 alone accounts for the entire ~27.5% energy change, constant
    // for all 100k steps after; with it, drift is ~1e-9 throughout).
    constrain_velocities(
        &topology, &atoms, &state.pos_x, &state.pos_y, &state.pos_z, &mut state.vel_x,
        &mut state.vel_y, &mut state.vel_z, 100, 1e-10, 0,
    )
    .unwrap();

    let kinetic_energy = |atoms: &AtomData, vx: &[f64], vy: &[f64], vz: &[f64]| -> f64 {
        (0..2)
            .map(|i| 0.5 * atoms.mass[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]))
            .sum()
    };

    let e0 = kinetic_energy(&atoms, &state.vel_x, &state.vel_y, &state.vel_z);
    assert!(e0 > 0.0, "test needs nonzero initial kinetic energy to be meaningful");

    let dt = 0.004;
    for _ in 0..100_000 {
        half_kick(&mut state, &atoms, dt / 2.0);
        drift_and_constrain(&mut state, &topology, &atoms, dt, 100, 1e-10).unwrap();
        half_kick(&mut state, &atoms, dt / 2.0);
        constrain_velocities(
            &topology, &atoms, &state.pos_x, &state.pos_y, &state.pos_z, &mut state.vel_x,
            &mut state.vel_y, &mut state.vel_z, 100, 1e-10, state.step,
        )
        .unwrap();
        state.step += 1;
    }

    let e1 = kinetic_energy(&atoms, &state.vel_x, &state.vel_y, &state.vel_z);
    let rel_drift = (e1 - e0).abs() / e0.abs();
    assert!(rel_drift < 1e-4, "relative energy drift {rel_drift} exceeds SAD.md §13.4's 1e-4 tolerance");

    let dx = state.pos_x[0] - state.pos_x[1];
    let dy = state.pos_y[0] - state.pos_y[1];
    let dz = state.pos_z[0] - state.pos_z[1];
    let bond_sq = dx * dx + dy * dy + dz * dz;
    assert!((bond_sq - topology.constr_dsq[0]).abs() < 1e-6, "bond drifted off the manifold");
}

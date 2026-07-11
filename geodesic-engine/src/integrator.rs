use crate::constraint;
use geodesic_core::{AtomData, BondedTopology, ConvergenceError, SimState};

/// Converts F/m from kcal/(mol·Å·amu) to Å/ps² so a kick with dt in ps is
/// dimensionally correct: 1 kcal/(mol·Å·amu) = 20.455² Å/ps², where 20.455 is
/// the AKMA time constant the DCD and inpcrd velocity units also use. Without
/// it the dynamics run on the wrong timescale and kinetic and potential energy
/// are not comparable (V stays kcal/mol while KE would be 418.4x too small).
pub const FORCE_TO_ACCEL_ANG_PER_PS2: f64 = 20.455 * 20.455;

/// B: velocity half-kick, v <- v + dt_half * M^-1 * F (SAD.md §2.3). Uses
/// whatever forces are currently in `state.force_{x,y,z}` — the caller is
/// responsible for having evaluated them at the matching positions first.
pub fn half_kick(state: &mut SimState, atoms: &AtomData, dt_half: f64) {
    for i in 0..state.pos_x.len() {
        let inv_m = 1.0 / atoms.mass[i];
        state.vel_x[i] += dt_half * state.force_x[i] * inv_m * FORCE_TO_ACCEL_ANG_PER_PS2;
        state.vel_y[i] += dt_half * state.force_y[i] * inv_m * FORCE_TO_ACCEL_ANG_PER_PS2;
        state.vel_z[i] += dt_half * state.force_z[i] * inv_m * FORCE_TO_ACCEL_ANG_PER_PS2;
    }
}

/// A: geodesic drift on the constraint manifold C (SAD.md §2.3). Takes the
/// unconstrained step r <- r + dt*v, projects back onto C via
/// `constraint::solve`, then resyncs velocity to the actual (constrained)
/// displacement, v <- (r_new - r_old) / dt. For an unconstrained atom this
/// resync is a no-op (r_new - r_old = dt*v exactly); for a constrained atom
/// it is exactly the SHAKE position correction divided by dt — the same
/// velocity update RATTLE's position stage applies implicitly (Andersen
/// 1983). It does not, by itself, remove any along-bond velocity a later
/// force kick introduces — see `constraint::constrain_velocities` for that.
pub fn drift_and_constrain(
    state: &mut SimState,
    topology: &BondedTopology,
    atoms: &AtomData,
    dt: f64,
    max_iter: u32,
    tol: f64,
) -> Result<u32, ConvergenceError> {
    let ref_x = state.pos_x.clone();
    let ref_y = state.pos_y.clone();
    let ref_z = state.pos_z.clone();

    for i in 0..state.pos_x.len() {
        state.pos_x[i] += dt * state.vel_x[i];
        state.pos_y[i] += dt * state.vel_y[i];
        state.pos_z[i] += dt * state.vel_z[i];
    }

    let iters = constraint::solve(
        topology,
        atoms,
        &ref_x,
        &ref_y,
        &ref_z,
        &mut state.pos_x,
        &mut state.pos_y,
        &mut state.pos_z,
        max_iter,
        tol,
        state.step,
    )?;

    for i in 0..state.pos_x.len() {
        state.vel_x[i] = (state.pos_x[i] - ref_x[i]) / dt;
        state.vel_y[i] = (state.pos_y[i] - ref_y[i]) / dt;
        state.vel_z[i] = (state.pos_z[i] - ref_z[i]) / dt;
    }

    Ok(iters)
}

use geodesic_core::{AtomData, BondedTopology, ConvergenceError, Element};

/// Moves bonds touching a hydrogen atom out of the harmonic bond list and
/// into `constr_i`/`constr_j`/`constr_dsq` as rigid holonomic constraints
/// for the Geodesic BAOAB A-step (SAD.md §2.3, §7.2). Mirrors AMBER's SHAKE
/// convention (`ntc=2`): a prmtop stores the full bond list regardless, and
/// whether X-H bonds are rigidified is a runtime policy decision
/// (`config.toml`'s `constraints.constrain_hydrogens`), not something the
/// file format itself encodes.
pub fn promote_hydrogen_bonds(topology: &mut BondedTopology, atoms: &AtomData) {
    let n = topology.bond_i.len();
    let mut keep_i = Vec::with_capacity(n);
    let mut keep_j = Vec::with_capacity(n);
    let mut keep_k = Vec::with_capacity(n);
    let mut keep_r0 = Vec::with_capacity(n);

    for idx in 0..n {
        let i = topology.bond_i[idx] as usize;
        let j = topology.bond_j[idx] as usize;
        let r0 = topology.bond_r0[idx];
        let touches_hydrogen =
            atoms.meta[i].element == Element::H || atoms.meta[j].element == Element::H;

        if touches_hydrogen {
            topology.constr_i.push(topology.bond_i[idx]);
            topology.constr_j.push(topology.bond_j[idx]);
            topology.constr_dsq.push(r0 * r0);
        } else {
            keep_i.push(topology.bond_i[idx]);
            keep_j.push(topology.bond_j[idx]);
            keep_k.push(topology.bond_k[idx]);
            keep_r0.push(r0);
        }
    }

    topology.bond_i = keep_i;
    topology.bond_j = keep_j;
    topology.bond_k = keep_k;
    topology.bond_r0 = keep_r0;
}

/// Iterative Lagrangian solver projecting trial positions back onto the
/// constraint manifold C = { |r_i - r_j|^2 = d0^2 for every constrained
/// pair } (SAD.md §2.2, §2.3). `ref_{x,y,z}` are the on-manifold positions
/// before the drift that produced `pos_{x,y,z}`; the reference bond vector
/// fixes each correction's direction (standard SHAKE).
///
/// Corrections within one iteration are computed from positions fixed at
/// the start of that iteration and summed per atom before being applied —
/// Jacobi relaxation, not in-place Gauss-Seidel — so the result is
/// independent of constraint order or how constraints are chunked across
/// threads (SAD.md §7.2: dispatch via `rayon::par_iter` with a deterministic
/// reduction must not change the answer).
///
/// Convergence is measured on max|lambda_i| across all constraints
/// (SAD.md §7.2, §10.2's `constr_tol`). Returns the iteration count on
/// success. After `max_iter` iterations without convergence, returns
/// `ConvergenceError` rather than silently returning an unconverged result
/// (SAD.md §12.2) — the position buffers are left in their last, still
/// partially-corrected state, which callers must treat as invalid.
#[allow(clippy::too_many_arguments)]
pub fn solve(
    topology: &BondedTopology,
    atoms: &AtomData,
    ref_x: &[f64],
    ref_y: &[f64],
    ref_z: &[f64],
    pos_x: &mut [f64],
    pos_y: &mut [f64],
    pos_z: &mut [f64],
    max_iter: u32,
    tol: f64,
    step: u64,
) -> Result<u32, ConvergenceError> {
    let n_constr = topology.constr_i.len();
    if n_constr == 0 {
        return Ok(0);
    }

    let n_atoms = pos_x.len();
    let mut delta_x = vec![0.0; n_atoms];
    let mut delta_y = vec![0.0; n_atoms];
    let mut delta_z = vec![0.0; n_atoms];

    let mut last_max_lambda = 0.0_f64;
    let mut last_worst = (0usize, 0usize, 0usize);

    for iter in 0..max_iter {
        delta_x.fill(0.0);
        delta_y.fill(0.0);
        delta_z.fill(0.0);

        let mut max_lambda = 0.0_f64;
        let mut worst = (0usize, 0usize, 0usize);

        for n in 0..n_constr {
            let i = topology.constr_i[n] as usize;
            let j = topology.constr_j[n] as usize;
            let dsq_target = topology.constr_dsq[n];

            let dx = pos_x[i] - pos_x[j];
            let dy = pos_y[i] - pos_y[j];
            let dz = pos_z[i] - pos_z[j];
            let diff = dx * dx + dy * dy + dz * dz - dsq_target;

            let rx = ref_x[i] - ref_x[j];
            let ry = ref_y[i] - ref_y[j];
            let rz = ref_z[i] - ref_z[j];

            let inv_mi = 1.0 / atoms.mass[i];
            let inv_mj = 1.0 / atoms.mass[j];
            let dot_new_old = dx * rx + dy * ry + dz * rz;
            let denom = 2.0 * (inv_mi + inv_mj) * dot_new_old;
            let lambda = if denom == 0.0 { 0.0 } else { diff / denom };

            delta_x[i] -= lambda * inv_mi * rx;
            delta_y[i] -= lambda * inv_mi * ry;
            delta_z[i] -= lambda * inv_mi * rz;
            delta_x[j] += lambda * inv_mj * rx;
            delta_y[j] += lambda * inv_mj * ry;
            delta_z[j] += lambda * inv_mj * rz;

            let lambda_abs = lambda.abs();
            if lambda_abs > max_lambda {
                max_lambda = lambda_abs;
                worst = (n, i, j);
            }
        }

        for a in 0..n_atoms {
            pos_x[a] += delta_x[a];
            pos_y[a] += delta_y[a];
            pos_z[a] += delta_z[a];
        }

        if max_lambda < tol {
            return Ok(iter + 1);
        }

        last_max_lambda = max_lambda;
        last_worst = worst;
    }

    Err(ConvergenceError::ConstraintSolver {
        step,
        constraint_idx: last_worst.0,
        atom_i: last_worst.1,
        atom_j: last_worst.2,
        residual: last_max_lambda,
        max_iter,
    })
}

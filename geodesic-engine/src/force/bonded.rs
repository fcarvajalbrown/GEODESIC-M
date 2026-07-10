use geodesic_core::{BondedTopology, SimState};

/// V_b = k(r - r0)^2 (SAD.md §2.1.1). Returns total bond potential energy,
/// accumulates into the force slices.
pub fn compute_bond_forces(
    state: &SimState,
    bonded: &BondedTopology,
    force_x: &mut [f64],
    force_y: &mut [f64],
    force_z: &mut [f64],
) -> f64 {
    let mut potential = 0.0;
    for n in 0..bonded.bond_i.len() {
        let i = bonded.bond_i[n] as usize;
        let j = bonded.bond_j[n] as usize;
        let k = bonded.bond_k[n];
        let r0 = bonded.bond_r0[n];

        let dx = state.pos_x[j] - state.pos_x[i];
        let dy = state.pos_y[j] - state.pos_y[i];
        let dz = state.pos_z[j] - state.pos_z[i];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r == 0.0 {
            continue;
        }

        let dv_dr = 2.0 * k * (r - r0);
        potential += k * (r - r0) * (r - r0);

        let (ux, uy, uz) = (dx / r, dy / r, dz / r);
        // dV/dr > 0 (stretched) pulls i toward j
        force_x[i] += dv_dr * ux;
        force_y[i] += dv_dr * uy;
        force_z[i] += dv_dr * uz;
        force_x[j] -= dv_dr * ux;
        force_y[j] -= dv_dr * uy;
        force_z[j] -= dv_dr * uz;
    }
    potential
}

/// V_theta = k(theta - theta0)^2, theta at the central atom j of i-j-k
/// (SAD.md §2.1.1).
pub fn compute_angle_forces(
    state: &SimState,
    bonded: &BondedTopology,
    force_x: &mut [f64],
    force_y: &mut [f64],
    force_z: &mut [f64],
) -> f64 {
    let mut potential = 0.0;
    for n in 0..bonded.angle_i.len() {
        let i = bonded.angle_i[n] as usize;
        let j = bonded.angle_j[n] as usize;
        let k_idx = bonded.angle_k[n] as usize;
        let k_const = bonded.angle_kth[n];
        let theta0 = bonded.angle_th0[n];

        let b1 = [
            state.pos_x[i] - state.pos_x[j],
            state.pos_y[i] - state.pos_y[j],
            state.pos_z[i] - state.pos_z[j],
        ];
        let b2 = [
            state.pos_x[k_idx] - state.pos_x[j],
            state.pos_y[k_idx] - state.pos_y[j],
            state.pos_z[k_idx] - state.pos_z[j],
        ];
        let r1 = (b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]).sqrt();
        let r2 = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
        if r1 == 0.0 || r2 == 0.0 {
            continue;
        }

        let dot = b1[0] * b2[0] + b1[1] * b2[1] + b1[2] * b2[2];
        let cos_theta = (dot / (r1 * r2)).clamp(-1.0, 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-8);
        let theta = cos_theta.acos();

        let dv_dtheta = 2.0 * k_const * (theta - theta0);
        potential += k_const * (theta - theta0) * (theta - theta0);

        let a = dv_dtheta / sin_theta;
        let mut f_i = [0.0; 3];
        let mut f_k = [0.0; 3];
        for c in 0..3 {
            f_i[c] = a * (b2[c] / (r1 * r2) - cos_theta * b1[c] / (r1 * r1));
            f_k[c] = a * (b1[c] / (r1 * r2) - cos_theta * b2[c] / (r2 * r2));
        }

        force_x[i] += f_i[0];
        force_y[i] += f_i[1];
        force_z[i] += f_i[2];
        force_x[k_idx] += f_k[0];
        force_y[k_idx] += f_k[1];
        force_z[k_idx] += f_k[2];
        force_x[j] -= f_i[0] + f_k[0];
        force_y[j] -= f_i[1] + f_k[1];
        force_z[j] -= f_i[2] + f_k[2];
    }
    potential
}

/// V_phi = k[1 + cos(n*phi - delta)] over i-j-k-l (SAD.md §2.1.1).
///
/// f_i and f_l come from the Blondel & Karplus (1996) plane-normal formula.
/// f_j and f_k are recovered from f_i/f_l via p = b1·b2/|b2|^2 and
/// q = b3·b2/|b2|^2: f_j = -(1+p)*f_i + q*f_l, f_k = p*f_i - (1+q)*f_l.
/// Verified by full symbolic chain-rule differentiation of
/// phi = atan2((m×n)·b2/|b2|, m·n) against three independent geometries
/// (machine-epsilon agreement) — see memory.md for the derivation.
pub fn compute_dihedral_forces(
    state: &SimState,
    bonded: &BondedTopology,
    force_x: &mut [f64],
    force_y: &mut [f64],
    force_z: &mut [f64],
) -> f64 {
    let mut potential = 0.0;
    for n in 0..bonded.dihed_i.len() {
        let i = bonded.dihed_i[n] as usize;
        let j = bonded.dihed_j[n] as usize;
        let k = bonded.dihed_k[n] as usize;
        let l = bonded.dihed_l[n] as usize;
        let k_phi = bonded.dihed_kphi[n];
        let mult = bonded.dihed_n[n] as f64;
        let delta = bonded.dihed_delta[n];

        let b1 = sub(pos(state, j), pos(state, i));
        let b2 = sub(pos(state, k), pos(state, j));
        let b3 = sub(pos(state, l), pos(state, k));

        let m = cross(b1, b2);
        let nvec = cross(b2, b3);
        let m_len2 = dot(m, m);
        let n_len2 = dot(nvec, nvec);
        let b2_len = dot(b2, b2).sqrt();
        if m_len2 == 0.0 || n_len2 == 0.0 || b2_len == 0.0 {
            continue;
        }

        let x = dot(m, nvec);
        let y = dot(cross(m, nvec), b2) / b2_len;
        let phi = y.atan2(x);

        potential += k_phi * (1.0 + (mult * phi - delta).cos());
        let dv_dphi = k_phi * mult * (mult * phi - delta).sin();

        // Blondel & Karplus (1996): forces via the two plane normals.
        let f_i = scale(m, -dv_dphi * b2_len / m_len2);
        let f_l = scale(nvec, dv_dphi * b2_len / n_len2);

        let p = dot(b1, b2) / (b2_len * b2_len);
        let q = dot(b3, b2) / (b2_len * b2_len);
        let f_j = add(scale(f_i, -(p + 1.0)), scale(f_l, q));
        let f_k = sub(scale(f_i, p), scale(f_l, q + 1.0));

        add_force(force_x, force_y, force_z, i, f_i);
        add_force(force_x, force_y, force_z, j, f_j);
        add_force(force_x, force_y, force_z, k, f_k);
        add_force(force_x, force_y, force_z, l, f_l);
    }
    potential
}

type Vec3 = [f64; 3];

fn pos(state: &SimState, i: usize) -> Vec3 {
    [state.pos_x[i], state.pos_y[i], state.pos_z[i]]
}
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn add_force(fx: &mut [f64], fy: &mut [f64], fz: &mut [f64], idx: usize, f: Vec3) {
    fx[idx] += f[0];
    fy[idx] += f[1];
    fz[idx] += f[2];
}

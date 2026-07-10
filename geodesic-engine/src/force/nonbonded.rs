use geodesic_core::{AtomData, NeighborList, SimState};

/// Lennard-Jones with a quintic switching function over [r_switch, r_cutoff]
/// (SAD.md §2.1.2). S(u) = 1 - 10u^3 + 15u^4 - 6u^5, u = (r - r_switch) /
/// (r_cutoff - r_switch) — zero value AND first/second derivative at both
/// ends, so the switched force has no discontinuity at either boundary.
/// Lorentz-Berthelot combining rules for cross-pair sigma/epsilon.
///
/// Accumulates into the given output slices and returns total potential
/// energy. Operates serially over whatever pairs it's given — parallel
/// partitioning across threads is cpu_backend.rs's job (SAD.md §7.2), not
/// this module's.
pub fn compute_pair_forces(
    state: &SimState,
    atoms: &AtomData,
    neighbor_list: &NeighborList,
    box_size: [f64; 3],
    force_x: &mut [f64],
    force_y: &mut [f64],
    force_z: &mut [f64],
) -> f64 {
    let r_cutoff = neighbor_list.r_cutoff;
    let r_switch = neighbor_list.r_switch;
    let mut potential = 0.0;

    for (&i, &j) in neighbor_list.pair_i.iter().zip(neighbor_list.pair_j.iter()) {
        let (i, j) = (i as usize, j as usize);

        let dx = min_image(state.pos_x[j] - state.pos_x[i], box_size[0]);
        let dy = min_image(state.pos_y[j] - state.pos_y[i], box_size[1]);
        let dz = min_image(state.pos_z[j] - state.pos_z[i], box_size[2]);
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 > r_cutoff * r_cutoff || r2 == 0.0 {
            continue;
        }
        let r = r2.sqrt();

        let sigma = 0.5 * (atoms.sigma[i] + atoms.sigma[j]);
        let epsilon = (atoms.epsilon[i] * atoms.epsilon[j]).sqrt();
        if epsilon == 0.0 {
            continue;
        }

        let sr6 = (sigma / r).powi(6);
        let sr12 = sr6 * sr6;
        let v_lj = 4.0 * epsilon * (sr12 - sr6);
        let f_lj = 24.0 * epsilon / r * (2.0 * sr12 - sr6); // -dV/dr, radial magnitude

        let (v, f_radial) = if r <= r_switch {
            (v_lj, f_lj)
        } else {
            let u = (r - r_switch) / (r_cutoff - r_switch);
            let u2 = u * u;
            let s = 1.0 - 10.0 * u2 * u + 15.0 * u2 * u2 - 6.0 * u2 * u2 * u;
            let ds_dr = -30.0 * u2 * (1.0 - u) * (1.0 - u) / (r_cutoff - r_switch);
            // F = -d(V*S)/dr = F_lj*S - V*dS/dr
            (v_lj * s, f_lj * s - v_lj * ds_dr)
        };

        potential += v;

        let fx = f_radial * dx / r;
        let fy = f_radial * dy / r;
        let fz = f_radial * dz / r;
        force_x[i] -= fx;
        force_y[i] -= fy;
        force_z[i] -= fz;
        force_x[j] += fx;
        force_y[j] += fy;
        force_z[j] += fz;
    }

    potential
}

fn min_image(d: f64, box_len: f64) -> f64 {
    d - box_len * (d / box_len).round()
}

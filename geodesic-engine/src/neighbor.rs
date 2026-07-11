use geodesic_core::{BondedTopology, NeighborList, SimParams, SimState};
use std::collections::HashSet;

/// Wraps positions into [0, box_size) (SAD.md §2.4), then builds the
/// Verlet pair list using minimum-image distances, excluding 1-2/1-3/1-4
/// bonded pairs (SAD.md §2.5). Serial by construction — list rebuilds are
/// infrequent, and a plain double loop is unambiguously deterministic.
pub fn build(state: &mut SimState, params: &SimParams, bonded: &BondedTopology) -> NeighborList {
    // Non-periodic systems (GBSA, vacuum) have no box to wrap into; wrapping
    // per-atom would split a molecule straddling a boundary, and the bonded
    // terms use raw (non-min-image) differences (SAD.md §2.4 applies only when
    // PBC is on). Leave coordinates whole in that case.
    if params.periodic {
        wrap_into_box(state, params.box_size);
    }

    let n = state.pos_x.len();
    let excluded: HashSet<(u32, u32)> = bonded
        .excl_i
        .iter()
        .zip(bonded.excl_j.iter())
        .map(|(&i, &j)| if i < j { (i, j) } else { (j, i) })
        .collect();

    let r_skin_sq = params.r_skin * params.r_skin;
    let mut pair_i = Vec::new();
    let mut pair_j = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i as u32, j as u32)) {
                continue;
            }
            let dx = min_image(state.pos_x[j] - state.pos_x[i], params.box_size[0]);
            let dy = min_image(state.pos_y[j] - state.pos_y[i], params.box_size[1]);
            let dz = min_image(state.pos_z[j] - state.pos_z[i], params.box_size[2]);
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 <= r_skin_sq {
                pair_i.push(i as u32);
                pair_j.push(j as u32);
            }
        }
    }

    NeighborList {
        pair_i,
        pair_j,
        ref_x: state.pos_x.clone(),
        ref_y: state.pos_y.clone(),
        ref_z: state.pos_z.clone(),
        r_cutoff: params.r_cutoff,
        r_skin: params.r_skin,
        r_switch: params.r_switch,
    }
}

/// SAD.md §8.4: rebuild when any atom displaces more than (r_skin -
/// r_cutoff) / 2 since the last build; squared displacement avoids a
/// square root on the hot path.
pub fn needs_rebuild(state: &SimState, list: &NeighborList) -> bool {
    let half_skin = (list.r_skin - list.r_cutoff) / 2.0;
    let threshold_sq = half_skin * half_skin;
    for i in 0..state.pos_x.len() {
        let dx = state.pos_x[i] - list.ref_x[i];
        let dy = state.pos_y[i] - list.ref_y[i];
        let dz = state.pos_z[i] - list.ref_z[i];
        if dx * dx + dy * dy + dz * dz > threshold_sq {
            return true;
        }
    }
    false
}

fn min_image(d: f64, box_len: f64) -> f64 {
    d - box_len * (d / box_len).round()
}

fn wrap_into_box(state: &mut SimState, box_size: [f64; 3]) {
    for i in 0..state.pos_x.len() {
        state.pos_x[i] = wrap_coord(state.pos_x[i], box_size[0]);
        state.pos_y[i] = wrap_coord(state.pos_y[i], box_size[1]);
        state.pos_z[i] = wrap_coord(state.pos_z[i], box_size[2]);
    }
}

fn wrap_coord(x: f64, box_len: f64) -> f64 {
    x - box_len * (x / box_len).floor()
}

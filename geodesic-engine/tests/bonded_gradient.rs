use geodesic_core::{BondedTopology, SimState};
use geodesic_engine::force::bonded::{compute_angle_forces, compute_bond_forces, compute_dihedral_forces};

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

// Generic finite-difference check: perturb each coordinate of each atom,
// compare -dV/dx (numeric) against the analytic force component.
fn check_gradient<F>(state: &SimState, compute: F)
where
    F: Fn(&SimState) -> (f64, Vec<f64>, Vec<f64>, Vec<f64>),
{
    let n = state.pos_x.len();
    let (_, fx, fy, fz) = compute(state);
    let eps = 1e-6;

    for i in 0..n {
        for (axis, analytic) in [(0, fx[i]), (1, fy[i]), (2, fz[i])] {
            let mut plus = clone_state(state);
            let mut minus = clone_state(state);
            match axis {
                0 => {
                    plus.pos_x[i] += eps;
                    minus.pos_x[i] -= eps;
                }
                1 => {
                    plus.pos_y[i] += eps;
                    minus.pos_y[i] -= eps;
                }
                _ => {
                    plus.pos_z[i] += eps;
                    minus.pos_z[i] -= eps;
                }
            }
            let (v_plus, ..) = compute(&plus);
            let (v_minus, ..) = compute(&minus);
            let numeric = -(v_plus - v_minus) / (2.0 * eps);
            let rel_err = (analytic - numeric).abs() / (numeric.abs() + 1.0);
            assert!(
                rel_err < 1e-4,
                "atom {i} axis {axis}: analytic={analytic}, numeric={numeric}, rel_err={rel_err}"
            );
        }
    }
}

fn clone_state(s: &SimState) -> SimState {
    let mut c = SimState::new(s.pos_x.len());
    c.pos_x = s.pos_x.clone();
    c.pos_y = s.pos_y.clone();
    c.pos_z = s.pos_z.clone();
    c
}

fn zero_forces(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    (vec![0.0; n], vec![0.0; n], vec![0.0; n])
}

#[test]
fn bond_gradient_matches_finite_difference() {
    let mut state = SimState::new(2);
    state.pos_x = vec![0.0, 1.7];
    let mut topo = empty_topology();
    topo.bond_i = vec![0];
    topo.bond_j = vec![1];
    topo.bond_k = vec![300.0];
    topo.bond_r0 = vec![1.5];

    check_gradient(&state, |s| {
        let (mut fx, mut fy, mut fz) = zero_forces(2);
        let v = compute_bond_forces(s, &topo, &mut fx, &mut fy, &mut fz);
        (v, fx, fy, fz)
    });
}

#[test]
fn angle_gradient_matches_finite_difference() {
    let mut state = SimState::new(3);
    // bent geometry, not collinear
    state.pos_x = vec![1.0, 0.0, 0.3];
    state.pos_y = vec![0.0, 0.0, 1.0];
    let mut topo = empty_topology();
    topo.angle_i = vec![0];
    topo.angle_j = vec![1];
    topo.angle_k = vec![2];
    topo.angle_kth = vec![50.0];
    topo.angle_th0 = vec![1.9106]; // ~109.5 degrees in radians

    check_gradient(&state, |s| {
        let (mut fx, mut fy, mut fz) = zero_forces(3);
        let v = compute_angle_forces(s, &topo, &mut fx, &mut fy, &mut fz);
        (v, fx, fy, fz)
    });
}

#[test]
fn dihedral_gradient_matches_finite_difference() {
    let mut state = SimState::new(4);
    // non-degenerate dihedral geometry (not collinear/coplanar)
    state.pos_x = vec![0.0, 0.0, 1.0, 1.5];
    state.pos_y = vec![0.0, 1.0, 1.0, 2.0];
    state.pos_z = vec![0.0, 0.0, 0.0, 0.7];
    let mut topo = empty_topology();
    topo.dihed_i = vec![0];
    topo.dihed_j = vec![1];
    topo.dihed_k = vec![2];
    topo.dihed_l = vec![3];
    topo.dihed_kphi = vec![2.0];
    topo.dihed_n = vec![2];
    topo.dihed_delta = vec![0.0];

    check_gradient(&state, |s| {
        let (mut fx, mut fy, mut fz) = zero_forces(4);
        let v = compute_dihedral_forces(s, &topo, &mut fx, &mut fy, &mut fz);
        (v, fx, fy, fz)
    });
}

#[test]
fn dihedral_gradient_matches_finite_difference_second_geometry() {
    // different, less symmetric geometry — catches a formula that only
    // happens to work for one special-case layout
    let mut state = SimState::new(4);
    state.pos_x = vec![0.3, -0.4, 0.9, 2.1];
    state.pos_y = vec![-0.6, 0.5, 1.3, 0.8];
    state.pos_z = vec![1.1, 0.2, -0.3, 1.4];
    let mut topo = empty_topology();
    topo.dihed_i = vec![0];
    topo.dihed_j = vec![1];
    topo.dihed_k = vec![2];
    topo.dihed_l = vec![3];
    topo.dihed_kphi = vec![3.5];
    topo.dihed_n = vec![3];
    topo.dihed_delta = vec![0.7];

    check_gradient(&state, |s| {
        let (mut fx, mut fy, mut fz) = zero_forces(4);
        let v = compute_dihedral_forces(s, &topo, &mut fx, &mut fy, &mut fz);
        (v, fx, fy, fz)
    });
}

#[test]
fn dihedral_total_force_sums_to_zero() {
    let mut state = SimState::new(4);
    state.pos_x = vec![0.0, 0.0, 1.0, 1.5];
    state.pos_y = vec![0.0, 1.0, 1.0, 2.0];
    state.pos_z = vec![0.0, 0.0, 0.0, 0.7];
    let mut topo = empty_topology();
    topo.dihed_i = vec![0];
    topo.dihed_j = vec![1];
    topo.dihed_k = vec![2];
    topo.dihed_l = vec![3];
    topo.dihed_kphi = vec![2.0];
    topo.dihed_n = vec![2];
    topo.dihed_delta = vec![0.0];

    let (mut fx, mut fy, mut fz) = zero_forces(4);
    compute_dihedral_forces(&state, &topo, &mut fx, &mut fy, &mut fz);
    assert!(fx.iter().sum::<f64>().abs() < 1e-10);
    assert!(fy.iter().sum::<f64>().abs() < 1e-10);
    assert!(fz.iter().sum::<f64>().abs() < 1e-10);
}

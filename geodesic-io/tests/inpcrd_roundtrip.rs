#[test]
fn positions_only() {
    let text = std::fs::read_to_string("tests/fixtures/lj_pair.inpcrd").unwrap();
    let state = geodesic_io::inpcrd::parse(&text, 2, false).unwrap();

    assert_eq!(state.pos_x, vec![0.0, 3.0]);
    assert_eq!(state.pos_y, vec![0.0, 0.0]);
    assert_eq!(state.pos_z, vec![0.0, 0.0]);
    assert_eq!(state.vel_x, vec![0.0, 0.0]);
}

#[test]
fn positions_and_velocities() {
    let text = std::fs::read_to_string("tests/fixtures/lj_pair_vel.inpcrd").unwrap();
    let state = geodesic_io::inpcrd::parse(&text, 2, false).unwrap();

    assert_eq!(state.pos_x, vec![0.0, 3.0]);
    // stored 1.0 * 20.455 Å per (1/20.455 ps) -> Å/ps
    assert!((state.vel_x[0] - 20.455).abs() < 1e-9);
    assert_eq!(state.vel_x[1], 0.0);
    assert_eq!(state.vel_y, vec![0.0, 0.0]);
}

#[test]
fn atom_count_mismatch_is_rejected() {
    let text = std::fs::read_to_string("tests/fixtures/lj_pair.inpcrd").unwrap();
    let err = geodesic_io::inpcrd::parse(&text, 3, false).unwrap_err();
    assert!(matches!(err, geodesic_core::ConfigError::PhysicallyInvalid { .. }));
}

use geodesic_core::Element;

#[test]
fn parse_positions_from_pdb() {
    let text = std::fs::read_to_string("tests/fixtures/lj_pair.pdb").unwrap();
    let (state, meta) = geodesic_io::pdb::parse_positions(&text).unwrap();

    assert_eq!(state.pos_x, vec![0.0, 3.0]);
    assert_eq!(state.pos_y, vec![0.0, 0.0]);
    assert_eq!(state.pos_z, vec![0.0, 0.0]);

    assert_eq!(meta.len(), 2);
    assert_eq!(meta[0].atom_name, *b"AR1 ");
    assert_eq!(meta[0].chain_id, b'A');
    assert_eq!(meta[0].residue_id, 1);
    // "AR" isn't in the Element enum yet -> Unknown, same as prmtop's inference
    assert_eq!(meta[0].element, Element::Unknown);
}

#[test]
fn write_and_reparse_snapshot() {
    let path = std::env::temp_dir().join("geodesic_pdb_snapshot_test.pdb");
    let text = std::fs::read_to_string("tests/fixtures/lj_pair.pdb").unwrap();
    let (state, meta) = geodesic_io::pdb::parse_positions(&text).unwrap();

    geodesic_io::pdb::write_snapshot(
        &path,
        &state.pos_x,
        &state.pos_y,
        &state.pos_z,
        &meta,
        Some(&[0.5, 1.5]),
    )
    .unwrap();

    let written = std::fs::read_to_string(&path).unwrap();
    let (state2, meta2) = geodesic_io::pdb::parse_positions(&written).unwrap();
    assert_eq!(state2.pos_x, state.pos_x);
    assert_eq!(state2.pos_z, state.pos_z);
    assert_eq!(meta2.len(), meta.len());

    let bfactor_line = written.lines().next().unwrap();
    assert!(bfactor_line.contains("0.50"));

    std::fs::remove_file(&path).ok();
}

use geodesic_core::Element;

#[test]
fn lj_pair_atom_data() {
    let text = std::fs::read_to_string("tests/fixtures/lj_pair.prmtop").unwrap();
    let (atoms, bonded) = geodesic_io::prmtop::parse(&text).unwrap();

    assert_eq!(atoms.mass.len(), 2);
    assert_eq!(atoms.mass[0], 1.0);
    assert_eq!(atoms.mass[1], 2.0);
    assert_eq!(atoms.charge, vec![0.0, 0.0]);

    // A=16384, B=256 -> sigma=(A/B)^(1/6)=2.0, epsilon=B^2/(4A)=1.0
    assert!((atoms.sigma[0] - 2.0).abs() < 1e-8);
    assert!((atoms.epsilon[0] - 1.0).abs() < 1e-8);
    assert!((atoms.sigma[1] - 2.0).abs() < 1e-8);
    assert!((atoms.epsilon[1] - 1.0).abs() < 1e-8);

    assert_eq!(atoms.meta[0].atom_name, *b"AR1 ");
    assert_eq!(atoms.meta[1].atom_name, *b"AR2 ");
    assert_eq!(atoms.meta[0].element, Element::Unknown);
    assert_eq!(atoms.meta[0].residue_id, 0);
    assert_eq!(atoms.meta[1].residue_id, 0);

    assert_eq!(bonded.bond_i.len(), 0);
    assert_eq!(bonded.angle_i.len(), 0);
    assert_eq!(bonded.dihed_i.len(), 0);
    assert_eq!(bonded.excl_i.len(), 0);
}

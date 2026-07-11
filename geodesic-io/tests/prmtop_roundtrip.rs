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

/// Real AMBER prmtop files store CHARGE pre-scaled by 18.2223 (so Coulomb
/// energy is q_i*q_j/r directly in kcal/mol) -- AtomData::charge must come
/// out in plain elementary-charge units, not the AMBER-internal scaling.
#[test]
fn charge_is_unscaled_from_amber_internal_units() {
    let text = std::fs::read_to_string("tests/fixtures/lj_pair.prmtop")
        .unwrap()
        .replace(
            "%FLAG CHARGE\n%FORMAT(5E16.8)\n  0.00000000E+00  0.00000000E+00",
            "%FLAG CHARGE\n%FORMAT(5E16.8)\n  7.28892000E+00 -7.28892000E+00",
        );
    let (atoms, _) = geodesic_io::prmtop::parse(&text).unwrap();

    // 7.28892 / 18.2223 = 0.4 elementary charge
    assert!((atoms.charge[0] - 0.4).abs() < 1e-8, "got {}", atoms.charge[0]);
    assert!((atoms.charge[1] + 0.4).abs() < 1e-8, "got {}", atoms.charge[1]);
}

/// SAD.md §13.1/§13.8: ala_dipeptide.prmtop is the standard small-molecule
/// benchmark. This is a real AmberTools-generated file (tleap, ff96,
/// `sequence { ACE ALA NME }`), sourced from
/// choderalab/YankTools/testsystems/data/alanine-dipeptide-gbsa -- not
/// hand-typed, so what's checked here is structural fact about the real
/// molecule (residue composition, atom count, bond count, net charge),
/// not fabricated force-field numbers.
#[test]
fn ala_dipeptide_atom_data() {
    let text = std::fs::read_to_string("tests/fixtures/ala_dipeptide.prmtop").unwrap();
    let (atoms, bonded) = geodesic_io::prmtop::parse(&text).unwrap();

    assert_eq!(atoms.mass.len(), 22, "ACE-ALA-NME has 22 atoms");
    assert_eq!(bonded.bond_i.len(), 21);

    let net_charge: f64 = atoms.charge.iter().sum();
    assert!(net_charge.abs() < 1e-6, "capped dipeptide should be net-neutral, got {net_charge}");

    // residue boundaries: ACE=1..7, ALA=7..17, NME=17..23 (1-based, from
    // RESIDUE_POINTER 1,7,17) -> atom 0 is in ACE, atom 21 is in NME
    assert_eq!(atoms.meta[0].residue_id, 0);
    assert_eq!(atoms.meta[21].residue_id, 2);
    assert_eq!(atoms.meta[0].element, Element::H); // HH31, first ACE methyl H
}

#[test]
fn ala_dipeptide_inpcrd_round_trips() {
    let prmtop_text = std::fs::read_to_string("tests/fixtures/ala_dipeptide.prmtop").unwrap();
    let (atoms, _) = geodesic_io::prmtop::parse(&prmtop_text).unwrap();
    let inpcrd_text = std::fs::read_to_string("tests/fixtures/ala_dipeptide.inpcrd").unwrap();
    let state = geodesic_io::inpcrd::parse(&inpcrd_text, atoms.mass.len(), false).unwrap();

    assert_eq!(state.pos_x.len(), 22);
    // first coordinate line of the real file: 2.0000010 1.0000000 -0.0000013
    assert!((state.pos_x[0] - 2.0000010).abs() < 1e-6);
    assert!((state.pos_y[0] - 1.0000000).abs() < 1e-6);
    assert!((state.pos_z[0] - (-0.0000013)).abs() < 1e-6);
}

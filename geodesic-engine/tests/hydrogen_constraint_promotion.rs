use geodesic_core::{AtomData, AtomMeta, BondedTopology, Element};
use geodesic_engine::constraint::promote_hydrogen_bonds;

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

fn meta(element: Element) -> AtomMeta {
    AtomMeta { element, residue_id: 0, atom_name: *b"X1  ", chain_id: 0 }
}

/// Methane-like fixture: C(0) bonded to H(1), H(2), H(3), plus one
/// non-hydrogen bond C(0)-N(4) that must stay a harmonic bond.
#[test]
fn only_bonds_touching_hydrogen_are_promoted() {
    let mut topo = empty_topology();
    topo.bond_i = vec![0, 0, 0, 0];
    topo.bond_j = vec![1, 2, 3, 4];
    topo.bond_k = vec![340.0, 340.0, 340.0, 400.0];
    topo.bond_r0 = vec![1.09, 1.09, 1.09, 1.47];

    let atoms = AtomData {
        epsilon: vec![0.0; 5],
        sigma: vec![0.0; 5],
        mass: vec![12.0, 1.008, 1.008, 1.008, 14.0],
        charge: vec![0.0; 5],
        meta: vec![
            meta(Element::C),
            meta(Element::H),
            meta(Element::H),
            meta(Element::H),
            meta(Element::N),
        ],
    };

    promote_hydrogen_bonds(&mut topo, &atoms);

    assert_eq!(topo.constr_i.len(), 3);
    assert_eq!(topo.constr_j, vec![1, 2, 3]);
    for &dsq in &topo.constr_dsq {
        assert!((dsq - 1.09 * 1.09).abs() < 1e-12);
    }

    assert_eq!(topo.bond_i, vec![0]);
    assert_eq!(topo.bond_j, vec![4]);
    assert_eq!(topo.bond_k, vec![400.0]);
    assert_eq!(topo.bond_r0, vec![1.47]);
}

#[test]
fn no_hydrogen_bonds_leaves_topology_unchanged() {
    let mut topo = empty_topology();
    topo.bond_i = vec![0];
    topo.bond_j = vec![1];
    topo.bond_k = vec![400.0];
    topo.bond_r0 = vec![1.47];

    let atoms = AtomData {
        epsilon: vec![0.0; 2],
        sigma: vec![0.0; 2],
        mass: vec![12.0, 14.0],
        charge: vec![0.0; 2],
        meta: vec![meta(Element::C), meta(Element::N)],
    };

    promote_hydrogen_bonds(&mut topo, &atoms);

    assert!(topo.constr_i.is_empty());
    assert_eq!(topo.bond_i, vec![0]);
    assert_eq!(topo.bond_j, vec![1]);
}

/// water_box_4 fixture (SAD.md §13.1's "PBC + angle forces + constraints"):
/// all 8 O-H bonds should promote to constraints, none of the 4 H-O-H
/// angles should be touched, and the promoted geometry should still be
/// solvable by constraint::solve.
#[test]
fn water_box_4_fixture_promotes_all_oh_bonds() {
    let prmtop_text = std::fs::read_to_string("tests/fixtures/water_box_4.prmtop").unwrap();
    let (atoms, mut topology) = geodesic_io::prmtop::parse(&prmtop_text).unwrap();
    let inpcrd_text = std::fs::read_to_string("tests/fixtures/water_box_4.inpcrd").unwrap();
    let mut state = geodesic_io::inpcrd::parse(&inpcrd_text, atoms.mass.len(), false).unwrap();

    assert_eq!(topology.bond_i.len(), 8);
    assert_eq!(topology.angle_i.len(), 4);

    promote_hydrogen_bonds(&mut topology, &atoms);

    assert_eq!(topology.constr_i.len(), 8, "all 8 O-H bonds should become constraints");
    assert!(topology.bond_i.is_empty(), "no harmonic bonds should remain");
    assert_eq!(topology.angle_i.len(), 4, "angle terms must not be touched by bond promotion");

    let ref_x = state.pos_x.clone();
    let ref_y = state.pos_y.clone();
    let ref_z = state.pos_z.clone();
    // perturb positions slightly, as a drift step would, then verify the
    // solver actually converges on this real 12-atom geometry
    for x in state.pos_x.iter_mut() {
        *x += 0.01;
    }
    let iters = geodesic_engine::constraint::solve(
        &topology, &atoms, &ref_x, &ref_y, &ref_z, &mut state.pos_x, &mut state.pos_y,
        &mut state.pos_z, 100, 1e-10, 0,
    )
    .unwrap();
    assert!(iters > 0 && iters <= 100);

    for n in 0..topology.constr_i.len() {
        let i = topology.constr_i[n] as usize;
        let j = topology.constr_j[n] as usize;
        let dx = state.pos_x[i] - state.pos_x[j];
        let dy = state.pos_y[i] - state.pos_y[j];
        let dz = state.pos_z[i] - state.pos_z[j];
        let dsq = dx * dx + dy * dy + dz * dz;
        assert!((dsq - topology.constr_dsq[n]).abs() < 1e-8, "constraint {n} not on manifold");
    }
}

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

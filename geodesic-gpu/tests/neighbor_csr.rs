use geodesic_gpu::neighbor_csr::build_csr;

#[test]
fn expands_half_list_to_symmetric_full_list() {
    // pairs: (0,1),(0,2),(1,2) over 3 atoms
    let (offsets, neighbors) = build_csr(&[0, 0, 1], &[1, 2, 2], 3);
    assert_eq!(offsets, vec![0, 2, 4, 6]);
    let slice = |a: usize| {
        let mut s = neighbors[offsets[a] as usize..offsets[a + 1] as usize].to_vec();
        s.sort();
        s
    };
    assert_eq!(slice(0), vec![1, 2]);
    assert_eq!(slice(1), vec![0, 2]);
    assert_eq!(slice(2), vec![0, 1]);
}

#[test]
fn isolated_atom_has_empty_slice() {
    let (offsets, neighbors) = build_csr(&[0], &[1], 3);
    assert_eq!(offsets, vec![0, 1, 2, 2]);
    assert_eq!(neighbors.len(), 2);
}

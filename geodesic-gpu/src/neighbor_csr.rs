/// Expand the CPU half pair list (i < j) into a full per-atom CSR gather list.
pub fn build_csr(pair_i: &[u32], pair_j: &[u32], n_atoms: usize) -> (Vec<u32>, Vec<u32>) {
    let mut degree = vec![0u32; n_atoms];
    for (&a, &b) in pair_i.iter().zip(pair_j.iter()) {
        degree[a as usize] += 1;
        degree[b as usize] += 1;
    }
    let mut offsets = vec![0u32; n_atoms + 1];
    for a in 0..n_atoms {
        offsets[a + 1] = offsets[a] + degree[a];
    }
    let total = offsets[n_atoms] as usize;
    let mut neighbors = vec![0u32; total];
    let mut cursor: Vec<u32> = offsets[..n_atoms].to_vec();
    for (&a, &b) in pair_i.iter().zip(pair_j.iter()) {
        let (ai, bi) = (a as usize, b as usize);
        neighbors[cursor[ai] as usize] = b;
        cursor[ai] += 1;
        neighbors[cursor[bi] as usize] = a;
        cursor[bi] += 1;
    }
    (offsets, neighbors)
}

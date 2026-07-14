mod common;
use common::{cpu_nonbonded_reference, load, params};
use geodesic_gpu::device;
use geodesic_gpu::kernel::NonbondedKernel;
use geodesic_gpu::neighbor_csr::build_csr;

#[test]
fn two_gpu_evaluations_are_bit_identical() {
    let Some(ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load("water_box_4");
    let p = params(state.pos_x.len());
    let (s, list, _fx, _fy, _fz) = cpu_nonbonded_reference(&state, &atoms, &topology, &p);
    let n = s.pos_x.len();
    let (offsets, neighbors) = build_csr(&list.pair_i, &list.pair_j, n);
    let mut kernel = NonbondedKernel::new(&ctx, &atoms, &p).unwrap();
    kernel.upload_neighbors(&ctx, &offsets, &neighbors);
    let (a, ea) = kernel.evaluate(&ctx, &s.pos_x, &s.pos_y, &s.pos_z);
    let (b, eb) = kernel.evaluate(&ctx, &s.pos_x, &s.pos_y, &s.pos_z);
    assert_eq!(a, b);
    assert_eq!(ea.to_bits(), eb.to_bits());
}

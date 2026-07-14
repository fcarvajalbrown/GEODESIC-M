mod common;
use common::{cpu_nonbonded_reference, load, params};
use geodesic_gpu::device;
use geodesic_gpu::kernel::{NonbondedInput, NonbondedKernel};
use geodesic_gpu::neighbor_csr::build_csr;

#[test]
fn two_gpu_evaluations_are_bit_identical() {
    let Some(ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load("water_box_4");
    let p = params(state.pos_x.len());
    let (s, list, _fx, _fy, _fz) = cpu_nonbonded_reference(&state, &atoms, &topology, &p);
    let n = s.pos_x.len();
    let (offsets, neighbors) = build_csr(&list.pair_i, &list.pair_j, n);
    let kernel = NonbondedKernel::new(&ctx).unwrap();
    let input = NonbondedInput {
        pos_x: &s.pos_x,
        pos_y: &s.pos_y,
        pos_z: &s.pos_z,
        sigma: &atoms.sigma,
        epsilon: &atoms.epsilon,
        offsets: &offsets,
        neighbors: &neighbors,
        r_cutoff: list.r_cutoff,
        r_switch: list.r_switch,
        box_size: p.box_size,
    };
    let (a, ea) = kernel.evaluate(&ctx, &input);
    let (b, eb) = kernel.evaluate(&ctx, &input);
    assert_eq!(a, b);
    assert_eq!(ea.to_bits(), eb.to_bits());
}

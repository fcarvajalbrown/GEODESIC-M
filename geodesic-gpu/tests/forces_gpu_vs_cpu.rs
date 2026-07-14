mod common;
use common::{cpu_nonbonded_reference, load, params};
use geodesic_gpu::device;
use geodesic_gpu::kernel::{NonbondedInput, NonbondedKernel};
use geodesic_gpu::neighbor_csr::build_csr;

fn kernel_matches_cpu(fixture: &str, tol: f32) {
    let Some(ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load(fixture);
    let p = params(state.pos_x.len());
    let (s, list, fx, fy, fz) = cpu_nonbonded_reference(&state, &atoms, &topology, &p);
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
    let (gpu_f, _e) = kernel.evaluate(&ctx, &input);
    for i in 0..n {
        for (c, cpu) in [(0usize, fx[i]), (1, fy[i]), (2, fz[i])] {
            let diff = (gpu_f[i][c] as f64 - cpu).abs();
            let bound = tol as f64 * cpu.abs().max(1.0);
            assert!(diff <= bound, "{fixture}: atom {i} comp {c}: gpu={}, cpu={cpu}, diff={diff}", gpu_f[i][c]);
        }
    }
}

#[test]
fn lj_pair_kernel_matches_cpu() {
    kernel_matches_cpu("lj_pair", 1e-4);
}

#[test]
fn water_box_4_kernel_matches_cpu() {
    kernel_matches_cpu("water_box_4", 1e-4);
}

#[test]
fn ala_dipeptide_full_backend_matches_cpu() {
    use common::{clone_positions, load, params};
    use geodesic_core::ComputeBackend;
    let Some(_ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load("ala_dipeptide");
    let (state2, atoms2, topology2) = load("ala_dipeptide");
    let p = params(state.pos_x.len());
    let n = state.pos_x.len();

    let mut s1 = clone_positions(&state);
    let mut s2 = clone_positions(&state2);

    let mut cpu = geodesic_engine::cpu_backend::CpuBackend::new(atoms, topology, &p);
    let mut gpu = geodesic_gpu::gpu_backend::GpuBackend::try_new(atoms2, topology2, &p).unwrap();

    cpu.build_neighbor_list(&mut s1, &p);
    let fc = cpu.compute_forces(&s1).clone();
    gpu.build_neighbor_list(&mut s2, &p);
    let fg = gpu.compute_forces(&s2).clone();

    for i in 0..n {
        for (cpuv, gpuv) in [(fc.fx[i], fg.fx[i]), (fc.fy[i], fg.fy[i]), (fc.fz[i], fg.fz[i])] {
            let diff = (cpuv - gpuv).abs();
            let bound = 1e-4 * cpuv.abs().max(1.0);
            assert!(diff <= bound, "atom {i}: cpu={cpuv}, gpu={gpuv}, diff={diff}");
        }
    }
}

//! SAD.md §13.9 benchmark suite. Numeric baselines are hardware-specific and
//! are captured on the pinned bench runner (§13.10), not committed from a dev
//! machine where they would be a meaningless cross-hardware gate; this file is
//! the harness that runner executes with `cargo bench -- --save-baseline main`.

use std::hint::black_box;
use std::path::PathBuf;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};

use geodesic_core::{AtomData, AtomMeta, BondedTopology, ComputeBackend, Element, SimParams, SimState};
use geodesic_engine::constraint::{self, promote_hydrogen_bonds};
use geodesic_engine::cpu_backend::CpuBackend;
use geodesic_engine::force::nonbonded;
use geodesic_engine::integrator::half_kick;
use geodesic_engine::neighbor;

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

fn meta() -> AtomMeta {
    AtomMeta { element: Element::C, residue_id: 0, atom_name: *b"C   ", chain_id: 0 }
}

fn lj_params(box_len: f64) -> SimParams {
    SimParams {
        n_atoms: 0,
        n_steps: 0,
        dt: 0.001,
        box_size: [box_len; 3],
        periodic: true,
        r_cutoff: 5.0,
        r_skin: 6.0,
        r_switch: 4.0,
        max_constr_iter: 100,
        constr_tol: 1e-8,
        frame_interval: 1,
        n_threads: 1,
        total_energy: 0.0,
    }
}

/// A cubic lattice of `n` LJ atoms at 2.0 Å spacing — deterministic, dense
/// enough to populate a realistic neighbor list.
fn lj_lattice(n: usize) -> (SimState, AtomData, SimParams) {
    let side = (n as f64).cbrt().ceil() as usize;
    let spacing = 2.0;
    let box_len = side as f64 * spacing;
    let mut state = SimState::new(n);
    let mut placed = 0;
    'outer: for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                if placed == n {
                    break 'outer;
                }
                state.pos_x[placed] = ix as f64 * spacing;
                state.pos_y[placed] = iy as f64 * spacing;
                state.pos_z[placed] = iz as f64 * spacing;
                placed += 1;
            }
        }
    }
    let atoms = AtomData {
        epsilon: vec![0.1; n],
        sigma: vec![3.0; n],
        mass: vec![12.0; n],
        charge: vec![0.0; n],
        meta: vec![meta(); n],
    };
    (state, atoms, lj_params(box_len))
}

fn bench_lj_inner_loop(c: &mut Criterion) {
    let (mut state, atoms, params) = lj_lattice(10_000);
    let list = neighbor::build(&mut state, &params, &empty_topology());
    let n = state.pos_x.len();
    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut fz = vec![0.0; n];
    c.bench_function("bench_lj_inner_loop", |b| {
        b.iter(|| {
            fx.fill(0.0);
            fy.fill(0.0);
            fz.fill(0.0);
            let e = nonbonded::compute_pair_forces(
                &state,
                &atoms,
                &list.pair_i,
                &list.pair_j,
                list.r_cutoff,
                list.r_switch,
                params.box_size,
                &mut fx,
                &mut fy,
                &mut fz,
            );
            black_box(e);
        });
    });
}

fn bench_neighbor_rebuild(c: &mut Criterion) {
    let (mut state, _atoms, params) = lj_lattice(10_000);
    let topo = empty_topology();
    c.bench_function("bench_neighbor_rebuild", |b| {
        b.iter(|| {
            let list = neighbor::build(&mut state, &params, &topo);
            black_box(list.pair_i.len());
        });
    });
}

/// 1000 independent rigid dimers (2000 atoms), each stretched 0.1 Å off its
/// 1.0 Å target so the A-step solver has real work to do.
fn bench_constraint_solver(c: &mut Criterion) {
    let n_constr = 1000;
    let n = n_constr * 2;
    let mut topo = empty_topology();
    let mut ref_x = vec![0.0; n];
    let ref_y = vec![0.0; n];
    let ref_z = vec![0.0; n];
    let mut pos_x = vec![0.0; n];
    for k in 0..n_constr {
        let (i, j) = (2 * k, 2 * k + 1);
        topo.constr_i.push(i as u32);
        topo.constr_j.push(j as u32);
        topo.constr_dsq.push(1.0);
        ref_x[i] = 4.0 * k as f64;
        ref_x[j] = 4.0 * k as f64 + 1.0;
        pos_x[i] = ref_x[i];
        pos_x[j] = ref_x[j] + 0.1; // stretched off the manifold
    }
    let atoms = AtomData {
        epsilon: vec![0.0; n],
        sigma: vec![0.0; n],
        mass: vec![1.008; n],
        charge: vec![0.0; n],
        meta: vec![meta(); n],
    };
    let pos_y = vec![0.0; n];
    let pos_z = vec![0.0; n];
    c.bench_function("bench_constraint_solver", |b| {
        b.iter_batched(
            || (pos_x.clone(), pos_y.clone(), pos_z.clone()),
            |(mut px, mut py, mut pz)| {
                let iters = constraint::solve(
                    &topo, &atoms, &ref_x, &ref_y, &ref_z, &mut px, &mut py, &mut pz, 100, 1e-10, 0,
                )
                .unwrap();
                black_box(iters);
            },
            BatchSize::SmallInput,
        );
    });
}

fn fixture(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("geodesic-engine")
        .join("tests")
        .join("fixtures")
        .join(name);
    std::fs::read_to_string(path).unwrap()
}

/// One full end-to-end BAB + RATTLE step on the real ala_dipeptide system.
fn bench_full_step(c: &mut Criterion) {
    let (atoms, mut topology) = geodesic_io::prmtop::parse(&fixture("ala_dipeptide.prmtop")).unwrap();
    let n = atoms.mass.len();
    let init = geodesic_io::inpcrd::parse(&fixture("ala_dipeptide.inpcrd"), n, false).unwrap();
    promote_hydrogen_bonds(&mut topology, &atoms);

    let mut params = lj_params(1000.0);
    params.n_atoms = n;
    params.periodic = false;
    params.r_cutoff = 12.0;
    params.r_skin = 14.0;
    params.r_switch = 10.0;

    let mut backend = CpuBackend::new(atoms, topology, &params);
    let dt = 0.001;

    let mut seed = SimState::new(n);
    seed.pos_x.copy_from_slice(&init.pos_x);
    seed.pos_y.copy_from_slice(&init.pos_y);
    seed.pos_z.copy_from_slice(&init.pos_z);
    backend.build_neighbor_list(&mut seed, &params);
    {
        let f = backend.compute_forces(&seed);
        seed.force_x.copy_from_slice(&f.fx);
        seed.force_y.copy_from_slice(&f.fy);
        seed.force_z.copy_from_slice(&f.fz);
    }

    c.bench_function("bench_full_step", |b| {
        b.iter_batched(
            || clone_state(&seed),
            |mut state| {
                half_kick(&mut state, backend.atoms(), dt / 2.0);
                backend.geodesic_drift(&mut state, dt).unwrap();
                if backend.needs_rebuild(&state) {
                    backend.build_neighbor_list(&mut state, &params);
                }
                {
                    let f = backend.compute_forces(&state);
                    state.force_x.copy_from_slice(&f.fx);
                    state.force_y.copy_from_slice(&f.fy);
                    state.force_z.copy_from_slice(&f.fz);
                }
                half_kick(&mut state, backend.atoms(), dt / 2.0);
                constraint::constrain_velocities(
                    backend.topology(),
                    backend.atoms(),
                    &state.pos_x,
                    &state.pos_y,
                    &state.pos_z,
                    &mut state.vel_x,
                    &mut state.vel_y,
                    &mut state.vel_z,
                    params.max_constr_iter,
                    params.constr_tol,
                    0,
                )
                .unwrap();
                black_box(state.pos_x[0]);
            },
            BatchSize::SmallInput,
        );
    });
}

fn clone_state(s: &SimState) -> SimState {
    let mut c = SimState::new(s.pos_x.len());
    c.pos_x.copy_from_slice(&s.pos_x);
    c.pos_y.copy_from_slice(&s.pos_y);
    c.pos_z.copy_from_slice(&s.pos_z);
    c.vel_x.copy_from_slice(&s.vel_x);
    c.vel_y.copy_from_slice(&s.vel_y);
    c.vel_z.copy_from_slice(&s.vel_z);
    c.force_x.copy_from_slice(&s.force_x);
    c.force_y.copy_from_slice(&s.force_y);
    c.force_z.copy_from_slice(&s.force_z);
    c
}

criterion_group!(
    benches,
    bench_lj_inner_loop,
    bench_neighbor_rebuild,
    bench_constraint_solver,
    bench_full_step
);
criterion_main!(benches);

use crate::force::{bonded, nonbonded};
use crate::integrator;
use crate::neighbor;
use geodesic_core::{
    AtomData, BondedTopology, ComputeBackend, ConvergenceError, ForceBuffer, NeighborList,
    SimParams, SimState,
};
use rayon::prelude::*;

/// CPU implementation of `ComputeBackend` (SAD.md §7.2): a fixed Rayon
/// thread pool of size `n_threads`, static (non-work-stealing) partitioning
/// of the non-bonded pair list into `n_threads` contiguous strips, private
/// per-thread accumulation, and a fixed thread-index-ordered sequential
/// reduction — deterministic regardless of how many threads run.
///
/// `geodesic_drift` implements only the position half of RATTLE (drift +
/// `constraint::solve` + velocity resync, via `integrator::drift_and_constrain`).
/// The velocity-tangency half (`constraint::constrain_velocities`) is not
/// part of the `ComputeBackend` trait — it must be invoked by whatever
/// drives the BAB loop (SAD.md §9.2, not yet built) immediately after the
/// second `integrator::half_kick` of a step, since it is plain,
/// backend-agnostic per-constraint arithmetic with no GPU-specific
/// consideration (SAD.md §7.4 doesn't list it as a distinct hybrid
/// workload).
pub struct CpuBackend {
    atoms: AtomData,
    topology: BondedTopology,
    neighbor_list: NeighborList,
    box_size: [f64; 3],
    n_threads: usize,
    max_constr_iter: u32,
    constr_tol: f64,
    thread_buffers: Vec<ForceBuffer>,
    bonded_buffer: ForceBuffer,
    reduced: ForceBuffer,
    potential_energy: f64,
}

fn zero_force_buffer(n: usize) -> ForceBuffer {
    ForceBuffer { fx: vec![0.0; n], fy: vec![0.0; n], fz: vec![0.0; n] }
}

impl CpuBackend {
    /// `params.n_threads == 0` auto-detects via
    /// `std::thread::available_parallelism` (SAD.md §10.2's "0 = auto-detect
    /// physical cores" — in practice this reports available logical
    /// parallelism, which may include SMT threads on some platforms; it is
    /// still deterministic for a given machine and is recorded as the
    /// reproducibility parameter T, per SAD.md §7.2).
    pub fn new(atoms: AtomData, topology: BondedTopology, params: &SimParams) -> Self {
        let n = atoms.mass.len();
        let n_threads = if params.n_threads == 0 {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
        } else {
            params.n_threads
        }
        .max(1);

        let neighbor_list = NeighborList {
            pair_i: Vec::new(),
            pair_j: Vec::new(),
            ref_x: vec![0.0; n],
            ref_y: vec![0.0; n],
            ref_z: vec![0.0; n],
            r_cutoff: params.r_cutoff,
            r_skin: params.r_skin,
            r_switch: params.r_switch,
        };

        Self {
            atoms,
            topology,
            neighbor_list,
            box_size: params.box_size,
            n_threads,
            max_constr_iter: params.max_constr_iter,
            constr_tol: params.constr_tol,
            thread_buffers: (0..n_threads).map(|_| zero_force_buffer(n)).collect(),
            bonded_buffer: zero_force_buffer(n),
            reduced: zero_force_buffer(n),
            potential_energy: 0.0,
        }
    }

    pub fn n_threads(&self) -> usize {
        self.n_threads
    }

    /// Total potential energy V (bond + angle + dihedral + LJ) accumulated by
    /// the most recent `compute_forces` call. Zero before the first call.
    pub fn potential_energy(&self) -> f64 {
        self.potential_energy
    }

    /// Read access to the moved-in atom data, so the BAB driver (the binary,
    /// SAD.md §9.2) can compute kinetic energy and call
    /// `constraint::constrain_velocities` without a second copy.
    pub fn atoms(&self) -> &AtomData {
        &self.atoms
    }

    /// Read access to the promoted topology, for the same reason as `atoms`.
    pub fn topology(&self) -> &BondedTopology {
        &self.topology
    }

    /// True when any atom has displaced far enough since the last rebuild that
    /// the Verlet list may be stale (SAD.md §12.5's loop guard).
    pub fn needs_rebuild(&self, state: &SimState) -> bool {
        neighbor::needs_rebuild(state, &self.neighbor_list)
    }
}

impl ComputeBackend for CpuBackend {
    fn build_neighbor_list(&mut self, state: &mut SimState, params: &SimParams) {
        self.neighbor_list = neighbor::build(state, params, &self.topology);
    }

    fn compute_forces(&mut self, state: &SimState) -> &ForceBuffer {
        let n = state.pos_x.len();

        // Bonded: sequential, low N relative to the pair list (SAD.md §7.4).
        self.bonded_buffer.fx.fill(0.0);
        self.bonded_buffer.fy.fill(0.0);
        self.bonded_buffer.fz.fill(0.0);
        let mut potential = 0.0;
        potential += bonded::compute_bond_forces(
            state,
            &self.topology,
            &mut self.bonded_buffer.fx,
            &mut self.bonded_buffer.fy,
            &mut self.bonded_buffer.fz,
        );
        potential += bonded::compute_angle_forces(
            state,
            &self.topology,
            &mut self.bonded_buffer.fx,
            &mut self.bonded_buffer.fy,
            &mut self.bonded_buffer.fz,
        );
        potential += bonded::compute_dihedral_forces(
            state,
            &self.topology,
            &mut self.bonded_buffer.fx,
            &mut self.bonded_buffer.fy,
            &mut self.bonded_buffer.fz,
        );

        // Non-bonded: static strip decomposition over the flat pair list
        // (SAD.md §7.2) — fixed contiguous chunks, not Rayon's default
        // work-stealing split, so the set of pairs each thread owns does
        // not depend on runtime scheduling.
        let atoms = &self.atoms;
        let pair_i = &self.neighbor_list.pair_i;
        let pair_j = &self.neighbor_list.pair_j;
        let r_cutoff = self.neighbor_list.r_cutoff;
        let r_switch = self.neighbor_list.r_switch;
        let box_size = self.box_size;
        let n_pairs = pair_i.len();
        let n_threads = self.thread_buffers.len();
        let chunk = n_pairs.div_ceil(n_threads).max(1);

        let mut thread_energy = vec![0.0; n_threads];
        self.thread_buffers
            .par_iter_mut()
            .zip(thread_energy.par_iter_mut())
            .enumerate()
            .for_each(|(k, (buf, energy))| {
                buf.fx.fill(0.0);
                buf.fy.fill(0.0);
                buf.fz.fill(0.0);
                let start = (k * chunk).min(n_pairs);
                let end = ((k + 1) * chunk).min(n_pairs);
                if start < end {
                    *energy = nonbonded::compute_pair_forces(
                        state,
                        atoms,
                        &pair_i[start..end],
                        &pair_j[start..end],
                        r_cutoff,
                        r_switch,
                        box_size,
                        &mut buf.fx,
                        &mut buf.fy,
                        &mut buf.fz,
                    );
                }
            });
        // Sum in fixed thread-index order, same determinism rationale as the
        // force reduction below (SAD.md §7.2).
        for &e in &thread_energy {
            potential += e;
        }
        self.potential_energy = potential;

        // Fixed sequential reduction, thread-index order (SAD.md §7.2).
        self.reduced.fx.copy_from_slice(&self.bonded_buffer.fx);
        self.reduced.fy.copy_from_slice(&self.bonded_buffer.fy);
        self.reduced.fz.copy_from_slice(&self.bonded_buffer.fz);
        for buf in &self.thread_buffers {
            for i in 0..n {
                self.reduced.fx[i] += buf.fx[i];
                self.reduced.fy[i] += buf.fy[i];
                self.reduced.fz[i] += buf.fz[i];
            }
        }

        &self.reduced
    }

    fn geodesic_drift(&mut self, state: &mut SimState, dt: f64) -> Result<(), ConvergenceError> {
        integrator::drift_and_constrain(
            state,
            &self.topology,
            &self.atoms,
            dt,
            self.max_constr_iter,
            self.constr_tol,
        )?;
        Ok(())
    }

    fn reduce_forces(&self) -> ForceBuffer {
        self.reduced.clone()
    }
}

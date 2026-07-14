use crate::device::{self, GpuContext};
use crate::kernel::NonbondedKernel;
use crate::neighbor_csr::build_csr;
use geodesic_core::{
    AtomData, BackendError, BondedTopology, ComputeBackend, ConvergenceError, ForceBuffer,
    NeighborList, SimParams, SimState,
};
use geodesic_engine::force::bonded;
use geodesic_engine::{integrator, neighbor};

pub struct GpuBackend {
    ctx: GpuContext,
    kernel: NonbondedKernel,
    atoms: AtomData,
    topology: BondedTopology,
    neighbor_list: NeighborList,
    max_constr_iter: u32,
    constr_tol: f64,
    reduced: ForceBuffer,
    potential_energy: f64,
}

impl GpuBackend {
    pub fn try_new(
        atoms: AtomData,
        topology: BondedTopology,
        params: &SimParams,
    ) -> Result<Self, BackendError> {
        let ctx = device::try_new()?;
        let kernel = NonbondedKernel::new(&ctx, &atoms, params)?;
        let n = atoms.mass.len();
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
        Ok(Self {
            ctx,
            kernel,
            atoms,
            topology,
            neighbor_list,
            max_constr_iter: params.max_constr_iter,
            constr_tol: params.constr_tol,
            reduced: ForceBuffer { fx: vec![0.0; n], fy: vec![0.0; n], fz: vec![0.0; n] },
            potential_energy: 0.0,
        })
    }
}

impl ComputeBackend for GpuBackend {
    fn build_neighbor_list(&mut self, state: &mut SimState, params: &SimParams) {
        self.neighbor_list = neighbor::build(state, params, &self.topology);
        let n = state.pos_x.len();
        let (off, nbr) = build_csr(&self.neighbor_list.pair_i, &self.neighbor_list.pair_j, n);
        self.kernel.upload_neighbors(&self.ctx, &off, &nbr);
    }

    fn compute_forces(&mut self, state: &SimState) -> &ForceBuffer {
        self.reduced.fx.iter_mut().for_each(|x| *x = 0.0);
        self.reduced.fy.iter_mut().for_each(|x| *x = 0.0);
        self.reduced.fz.iter_mut().for_each(|x| *x = 0.0);
        let mut potential = 0.0;
        potential += bonded::compute_bond_forces(state, &self.topology, &mut self.reduced.fx, &mut self.reduced.fy, &mut self.reduced.fz);
        potential += bonded::compute_angle_forces(state, &self.topology, &mut self.reduced.fx, &mut self.reduced.fy, &mut self.reduced.fz);
        potential += bonded::compute_dihedral_forces(state, &self.topology, &mut self.reduced.fx, &mut self.reduced.fy, &mut self.reduced.fz);

        let (gpu_f, nb_energy) = self.kernel.evaluate(&self.ctx, &state.pos_x, &state.pos_y, &state.pos_z);
        for (i, gf) in gpu_f.iter().enumerate() {
            self.reduced.fx[i] += gf[0] as f64;
            self.reduced.fy[i] += gf[1] as f64;
            self.reduced.fz[i] += gf[2] as f64;
        }
        potential += nb_energy as f64;
        self.potential_energy = potential;
        &self.reduced
    }

    fn geodesic_drift(&mut self, state: &mut SimState, dt: f64) -> Result<(), ConvergenceError> {
        integrator::drift_and_constrain(state, &self.topology, &self.atoms, dt, self.max_constr_iter, self.constr_tol)?;
        Ok(())
    }

    fn reduce_forces(&self) -> ForceBuffer {
        self.reduced.clone()
    }

    fn potential_energy(&self) -> f64 { self.potential_energy }
    fn atoms(&self) -> &AtomData { &self.atoms }
    fn topology(&self) -> &BondedTopology { &self.topology }
    fn needs_rebuild(&self, state: &SimState) -> bool {
        neighbor::needs_rebuild(state, &self.neighbor_list)
    }
    fn n_threads(&self) -> usize { 1 }
}

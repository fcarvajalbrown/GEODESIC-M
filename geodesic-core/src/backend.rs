use crate::atoms::AtomData;
use crate::buffers::ForceBuffer;
use crate::error::ConvergenceError;
use crate::params::SimParams;
use crate::state::SimState;
use crate::topology::BondedTopology;

/// Dispatches force evaluation and geodesic drift to CPU, GPU, or hybrid.
/// Selected once at startup; the simulation loop never inspects the
/// concrete type behind Box<dyn ComputeBackend>.
pub trait ComputeBackend: Send {
    // &mut SimState: rebuild wraps atom positions into [0, box_size) per
    // SAD.md §2.4, so it must be able to mutate positions, not just read them.
    fn build_neighbor_list(&mut self, state: &mut SimState, params: &SimParams);
    fn compute_forces(&mut self, state: &SimState) -> &ForceBuffer;
    fn geodesic_drift(&mut self, state: &mut SimState, dt: f64) -> Result<(), ConvergenceError>;
    fn reduce_forces(&self) -> ForceBuffer;
    fn potential_energy(&self) -> f64;
    fn atoms(&self) -> &AtomData;
    fn topology(&self) -> &BondedTopology;
    fn needs_rebuild(&self, state: &SimState) -> bool;
    fn n_threads(&self) -> usize;
}

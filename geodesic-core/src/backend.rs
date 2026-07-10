use crate::buffers::ForceBuffer;
use crate::error::ConvergenceError;
use crate::params::SimParams;
use crate::state::SimState;

/// Dispatches force evaluation and geodesic drift to CPU, GPU, or hybrid.
/// Selected once at startup; the simulation loop never inspects the
/// concrete type behind Box<dyn ComputeBackend>.
pub trait ComputeBackend: Send {
    fn build_neighbor_list(&mut self, state: &SimState, params: &SimParams);
    fn compute_forces(&mut self, state: &SimState) -> &ForceBuffer;
    fn geodesic_drift(&mut self, state: &mut SimState, dt: f64) -> Result<(), ConvergenceError>;
    fn reduce_forces(&self) -> ForceBuffer;
}

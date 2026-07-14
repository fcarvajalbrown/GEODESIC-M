use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Debug, Error)]
pub enum SimError {
    #[error(transparent)]
    Io(#[from] IoError),
    #[error(transparent)]
    Config(#[from] ConfigError),
    #[error(transparent)]
    Numerical(#[from] NumericalError),
    #[error(transparent)]
    Convergence(#[from] ConvergenceError),
    #[error(transparent)]
    Backend(#[from] BackendError),
    #[error(transparent)]
    Topology(#[from] TopologyError),
}

#[derive(Debug, Error)]
#[error("I/O error at '{}': {source}", path.display())]
pub struct IoError {
    pub path: std::path::PathBuf,
    #[source]
    pub source: std::io::Error,
}

#[derive(Debug, Error)]
pub enum NumericalError {
    #[error("NaN in force_{component:?} at step {step}, atom {atom}")]
    NanInForce { step: u64, atom: usize, component: Axis },
    #[error("NaN in position_{component:?} at step {step}, atom {atom}")]
    NanInPos { step: u64, atom: usize, component: Axis },
    #[error(
        "energy drift {drift_kcal} kcal/mol exceeds threshold {threshold_kcal} kcal/mol at step {step}"
    )]
    EnergyDrift {
        step: u64,
        drift_kcal: f64,
        threshold_kcal: f64,
    },
}

#[derive(Debug, Error)]
pub enum ConvergenceError {
    #[error(
        "constraint solver failed to converge for constraint {constraint_idx} (atoms {atom_i}, {atom_j}) at step {step}: residual {residual} after {max_iter} iterations"
    )]
    ConstraintSolver {
        step: u64,
        constraint_idx: usize,
        atom_i: usize,
        atom_j: usize,
        residual: f64,
        max_iter: u32,
    },
}

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("unknown config key: {0}")]
    UnknownKey(String),
    #[error("invalid value for key '{key}': '{value}' ({reason})")]
    InvalidValue {
        key: String,
        value: String,
        reason: String,
    },
    #[error("missing required config key: {0}")]
    MissingRequired(String),
    #[error("physically invalid configuration: {description}")]
    PhysicallyInvalid { description: String },
}

#[derive(Debug, Error)]
pub enum BackendError {
    #[error("GPU device lost")]
    DeviceLost,
    #[error("shader compilation failed: {0}")]
    ShaderCompilation(String),
    #[error("out of GPU memory")]
    OutOfGpuMemory,
    #[error("no compatible GPU adapter found (DX12/Vulkan)")]
    NoAdapter,
}

#[derive(Debug, Error)]
pub enum TopologyError {
    #[error("eigensolver failed: {reason}")]
    EigensolverFailed { reason: String },
    #[error("Ripser ran out of memory")]
    RipserOutOfMemory,
}

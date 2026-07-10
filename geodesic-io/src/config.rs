use geodesic_core::{ConfigError, SimParams};
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawConfig {
    run: RawRun,
    system: RawSystem,
    integrator: RawIntegrator,
    nonbonded: RawNonbonded,
    constraints: RawConstraints,
    output: RawOutput,
    monitoring: Option<RawMonitoring>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawRun {
    n_steps: u64,
    frame_interval: u32,
    backend: String,
    n_threads: usize,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSystem {
    prmtop: PathBuf,
    inpcrd: PathBuf,
    box_size: [f64; 3],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawIntegrator {
    dt: f64,
    total_energy: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawNonbonded {
    r_cutoff: f64,
    r_skin: f64,
    r_switch: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawConstraints {
    max_iter: u32,
    tolerance: f64,
    #[serde(default = "default_constrain_hydrogens")]
    constrain_hydrogens: bool,
}

fn default_constrain_hydrogens() -> bool {
    true
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawOutput {
    trajectory: PathBuf,
    energy_log: PathBuf,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawMonitoring {
    energy_drift_threshold_kcal: f64,
    energy_drift_action: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Cpu,
    Gpu,
    Hybrid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DriftAction {
    Warn,
    Stop,
}

pub struct MonitoringConfig {
    pub energy_drift_threshold_kcal: f64,
    pub energy_drift_action: DriftAction,
}

/// Fully validated config.toml contents. n_atoms is not known until the
/// prmtop is parsed, so this stops short of SimParams — call
/// `into_sim_params` once n_atoms is available.
pub struct Config {
    pub n_steps: u64,
    pub frame_interval: u32,
    pub backend: Backend,
    pub n_threads: usize,

    pub prmtop: PathBuf,
    pub inpcrd: PathBuf,
    pub box_size: [f64; 3],

    pub dt: f64,
    pub total_energy: f64,

    pub r_cutoff: f64,
    pub r_skin: f64,
    pub r_switch: f64,

    pub max_constr_iter: u32,
    pub constr_tol: f64,
    pub constrain_hydrogens: bool,

    pub trajectory: PathBuf,
    pub energy_log: PathBuf,

    pub monitoring: Option<MonitoringConfig>,
}

impl Config {
    pub fn from_toml_str(s: &str) -> Result<Config, ConfigError> {
        let raw: RawConfig = toml::from_str(s).map_err(map_toml_error)?;

        let backend = match raw.run.backend.as_str() {
            "cpu" => Backend::Cpu,
            "gpu" => Backend::Gpu,
            "hybrid" => Backend::Hybrid,
            other => {
                return Err(ConfigError::InvalidValue {
                    key: "run.backend".to_string(),
                    value: other.to_string(),
                    reason: "expected \"cpu\", \"gpu\", or \"hybrid\"".to_string(),
                })
            }
        };

        let monitoring = match raw.monitoring {
            None => None,
            Some(m) => {
                let energy_drift_action = match m.energy_drift_action.as_str() {
                    "warn" => DriftAction::Warn,
                    "stop" => DriftAction::Stop,
                    other => {
                        return Err(ConfigError::InvalidValue {
                            key: "monitoring.energy_drift_action".to_string(),
                            value: other.to_string(),
                            reason: "expected \"warn\" or \"stop\"".to_string(),
                        })
                    }
                };
                Some(MonitoringConfig {
                    energy_drift_threshold_kcal: m.energy_drift_threshold_kcal,
                    energy_drift_action,
                })
            }
        };

        let config = Config {
            n_steps: raw.run.n_steps,
            frame_interval: raw.run.frame_interval,
            backend,
            n_threads: raw.run.n_threads,
            prmtop: raw.system.prmtop,
            inpcrd: raw.system.inpcrd,
            box_size: raw.system.box_size,
            dt: raw.integrator.dt,
            total_energy: raw.integrator.total_energy,
            r_cutoff: raw.nonbonded.r_cutoff,
            r_skin: raw.nonbonded.r_skin,
            r_switch: raw.nonbonded.r_switch,
            max_constr_iter: raw.constraints.max_iter,
            constr_tol: raw.constraints.tolerance,
            constrain_hydrogens: raw.constraints.constrain_hydrogens,
            trajectory: raw.output.trajectory,
            energy_log: raw.output.energy_log,
            monitoring,
        };

        validate(&config)?;
        Ok(config)
    }

    pub fn into_sim_params(self, n_atoms: usize) -> SimParams {
        SimParams {
            n_atoms,
            n_steps: self.n_steps,
            dt: self.dt,
            box_size: self.box_size,
            r_cutoff: self.r_cutoff,
            r_skin: self.r_skin,
            r_switch: self.r_switch,
            max_constr_iter: self.max_constr_iter,
            constr_tol: self.constr_tol,
            frame_interval: self.frame_interval,
            n_threads: self.n_threads,
            total_energy: self.total_energy,
        }
    }
}

fn validate(c: &Config) -> Result<(), ConfigError> {
    if c.n_steps == 0 {
        return Err(ConfigError::PhysicallyInvalid {
            description: "run.n_steps must be greater than zero".to_string(),
        });
    }
    if c.box_size.iter().any(|&x| x <= 0.0) {
        return Err(ConfigError::PhysicallyInvalid {
            description: "system.box_size components must be positive".to_string(),
        });
    }
    if c.dt <= 0.0 {
        return Err(ConfigError::PhysicallyInvalid {
            description: "integrator.dt must be positive".to_string(),
        });
    }
    if !(c.r_switch < c.r_cutoff && c.r_cutoff < c.r_skin) {
        return Err(ConfigError::PhysicallyInvalid {
            description: format!(
                "nonbonded radii must satisfy r_switch < r_cutoff < r_skin (got r_switch={}, r_cutoff={}, r_skin={})",
                c.r_switch, c.r_cutoff, c.r_skin
            ),
        });
    }
    if c.max_constr_iter == 0 {
        return Err(ConfigError::PhysicallyInvalid {
            description: "constraints.max_iter must be greater than zero".to_string(),
        });
    }
    if c.constr_tol <= 0.0 {
        return Err(ConfigError::PhysicallyInvalid {
            description: "constraints.tolerance must be positive".to_string(),
        });
    }
    Ok(())
}

fn map_toml_error(e: toml::de::Error) -> ConfigError {
    let msg = e.message();
    if let Some(field) = msg.strip_prefix("missing field `").and_then(|s| s.strip_suffix('`')) {
        return ConfigError::MissingRequired(field.to_string());
    }
    if let Some(field) = msg.strip_prefix("unknown field `") {
        if let Some(end) = field.find('`') {
            return ConfigError::UnknownKey(field[..end].to_string());
        }
    }
    ConfigError::InvalidValue {
        key: "config.toml".to_string(),
        value: String::new(),
        reason: e.to_string(),
    }
}

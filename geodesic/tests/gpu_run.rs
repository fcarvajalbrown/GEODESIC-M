#![cfg(feature = "gpu")]

use std::path::{Path, PathBuf};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("geodesic-engine").join("tests").join("fixtures")
}
fn slashed(p: &Path) -> String {
    p.to_string_lossy().replace('\\', "/")
}

#[test]
fn gpu_backend_runs_ala_dipeptide_or_skips() {
    let fixtures = fixtures_dir();
    let prmtop = slashed(&fixtures.join("ala_dipeptide.prmtop"));
    let inpcrd = slashed(&fixtures.join("ala_dipeptide.inpcrd"));
    let out_dir = std::env::temp_dir();
    let dcd = slashed(&out_dir.join("geodesic_gpu_run.dcd"));
    let csv = slashed(&out_dir.join("geodesic_gpu_run.csv"));
    let cfg_path = out_dir.join("geodesic_gpu_run.toml");
    let config = format!(
        r#"
[run]
n_steps        = 20
frame_interval = 5
backend        = "gpu"
n_threads      = 1

[system]
prmtop   = "{prmtop}"
inpcrd   = "{inpcrd}"
box_size = [1000.0, 1000.0, 1000.0]
periodic = false

[integrator]
dt           = 0.001
total_energy = 5.12

[nonbonded]
r_cutoff = 12.0
r_skin   = 14.0
r_switch = 10.0

[constraints]
max_iter  = 100
tolerance = 1.0e-8

[output]
trajectory = "{dcd}"
energy_log = "{csv}"
"#
    );
    std::fs::write(&cfg_path, config).unwrap();

    match geodesic::run_from_config_file(&cfg_path) {
        Ok(summary) => {
            let bytes = std::fs::read(&summary.trajectory).unwrap();
            assert!(!bytes.is_empty(), "GPU run produced an empty DCD");
        }
        Err(geodesic_core::SimError::Backend(geodesic_core::BackendError::NoAdapter)) => {
            eprintln!("skipping GPU run test: no adapter available");
        }
        Err(e) => panic!("GPU run failed: {e}"),
    }
}

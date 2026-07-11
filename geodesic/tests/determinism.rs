//! SAD.md §13.5: two independent runs with identical SimParams and initial
//! state must produce bit-for-bit identical DCD output. This is the only guard
//! against a non-deterministic parallelism regression, so the runs use
//! n_threads > 1 to exercise the static strip decomposition and fixed-order
//! reduction in cpu_backend (SAD.md §7.2); an accidental unordered Rayon
//! reduction would diverge here.

use std::path::{Path, PathBuf};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("geodesic-engine")
        .join("tests")
        .join("fixtures")
}

fn slashed(p: &Path) -> String {
    p.to_string_lossy().replace('\\', "/")
}

fn write_config(tag: &str) -> PathBuf {
    let fixtures = fixtures_dir();
    let prmtop = slashed(&fixtures.join("ala_dipeptide.prmtop"));
    let inpcrd = slashed(&fixtures.join("ala_dipeptide.inpcrd"));
    let out_dir = std::env::temp_dir();
    let dcd = slashed(&out_dir.join(format!("geodesic_determinism_{tag}.dcd")));
    let csv = slashed(&out_dir.join(format!("geodesic_determinism_{tag}.csv")));
    let config_path = out_dir.join(format!("geodesic_determinism_{tag}.toml"));

    let config = format!(
        r#"
[run]
n_steps        = 50
frame_interval = 5
backend        = "cpu"
n_threads      = 4

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
    std::fs::write(&config_path, config).unwrap();
    config_path
}

#[test]
fn two_runs_produce_byte_identical_dcd() {
    let cfg_a = write_config("a");
    let cfg_b = write_config("b");

    let sum_a = geodesic::run_from_config_file(&cfg_a).unwrap();
    let sum_b = geodesic::run_from_config_file(&cfg_b).unwrap();
    assert_eq!(sum_a.n_frames, sum_b.n_frames);

    let bytes_a = std::fs::read(&sum_a.trajectory).unwrap();
    let bytes_b = std::fs::read(&sum_b.trajectory).unwrap();
    assert!(
        bytes_a == bytes_b,
        "two identical runs produced different DCD output ({} vs {} bytes) -- a non-deterministic reduction has crept in",
        bytes_a.len(),
        bytes_b.len()
    );
}

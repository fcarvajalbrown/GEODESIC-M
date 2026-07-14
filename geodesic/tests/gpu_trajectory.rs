//! ROADMAP v0.6 exit criterion: the GPU/hybrid backend matches the pure-CPU
//! energy trajectory on ala_dipeptide within the v0.5 f32 tolerance over a
//! multi-step run. Adapter-adaptive: skips when no DX12/Vulkan adapter exists.
#![cfg(feature = "gpu")]

use std::path::{Path, PathBuf};

const TRAJ_REL_TOL: f64 = 1e-4;

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

fn write_config(tag: &str, backend: &str) -> PathBuf {
    let fixtures = fixtures_dir();
    let prmtop = slashed(&fixtures.join("ala_dipeptide.prmtop"));
    let inpcrd = slashed(&fixtures.join("ala_dipeptide.inpcrd"));
    let out_dir = std::env::temp_dir();
    let dcd = slashed(&out_dir.join(format!("geodesic_traj_{tag}.dcd")));
    let csv = slashed(&out_dir.join(format!("geodesic_traj_{tag}.csv")));
    let config_path = out_dir.join(format!("geodesic_traj_{tag}.toml"));

    let config = format!(
        r#"
[run]
n_steps        = 5
frame_interval = 1
backend        = "{backend}"
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
    std::fs::write(&config_path, config).unwrap();
    config_path
}

fn rel_diff(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(1.0)
}

#[test]
fn gpu_trajectory_matches_cpu() {
    if geodesic_gpu::device::context_or_skip().is_none() {
        return;
    }
    let cfg_cpu = write_config("cpu", "cpu");
    let cfg_gpu = write_config("gpu", "gpu");

    let cpu = geodesic::run_from_config_file(&cfg_cpu).unwrap();
    let gpu = geodesic::run_from_config_file(&cfg_gpu).unwrap();

    assert_eq!(cpu.n_frames, gpu.n_frames, "cpu and gpu produced different frame counts");

    let dp = rel_diff(cpu.final_potential, gpu.final_potential);
    let dk = rel_diff(cpu.final_kinetic, gpu.final_kinetic);
    assert!(
        dp <= TRAJ_REL_TOL,
        "final potential energy diverged: cpu={}, gpu={}, rel={dp}",
        cpu.final_potential, gpu.final_potential
    );
    assert!(
        dk <= TRAJ_REL_TOL,
        "final kinetic energy diverged: cpu={}, gpu={}, rel={dk}",
        cpu.final_kinetic, gpu.final_kinetic
    );
}

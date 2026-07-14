//! ADR 0005: "hybrid" is a documented alias for "gpu". Both resolve to the same
//! optimized GpuBackend, so a run with backend="hybrid" must produce
//! byte-identical DCD to the same run with backend="gpu". Adapter-adaptive:
//! skips cleanly when no DX12/Vulkan adapter is present.
#![cfg(feature = "gpu")]

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

fn write_config(tag: &str, backend: &str) -> PathBuf {
    let fixtures = fixtures_dir();
    let prmtop = slashed(&fixtures.join("ala_dipeptide.prmtop"));
    let inpcrd = slashed(&fixtures.join("ala_dipeptide.inpcrd"));
    let out_dir = std::env::temp_dir();
    let dcd = slashed(&out_dir.join(format!("geodesic_alias_{tag}.dcd")));
    let csv = slashed(&out_dir.join(format!("geodesic_alias_{tag}.csv")));
    let config_path = out_dir.join(format!("geodesic_alias_{tag}.toml"));

    let config = format!(
        r#"
[run]
n_steps        = 10
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

#[test]
fn hybrid_is_byte_identical_to_gpu() {
    if geodesic_gpu::device::context_or_skip().is_none() {
        return;
    }
    let cfg_gpu = write_config("gpu", "gpu");
    let cfg_hybrid = write_config("hybrid", "hybrid");

    let sum_gpu = geodesic::run_from_config_file(&cfg_gpu).unwrap();
    let sum_hybrid = geodesic::run_from_config_file(&cfg_hybrid).unwrap();
    assert_eq!(sum_gpu.n_frames, sum_hybrid.n_frames);

    let bytes_gpu = std::fs::read(&sum_gpu.trajectory).unwrap();
    let bytes_hybrid = std::fs::read(&sum_hybrid.trajectory).unwrap();
    assert!(
        bytes_gpu == bytes_hybrid,
        "hybrid and gpu produced different DCD output ({} vs {} bytes) -- the alias is not resolving to the same path",
        bytes_gpu.len(),
        bytes_hybrid.len()
    );
}

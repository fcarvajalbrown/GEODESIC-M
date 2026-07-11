//! SAD.md §13.7: the frozen golden reference trajectory. `ala_dipeptide_ref.dcd`
//! is 100 frames from the first verified-correct build; every subsequent build
//! must reproduce it byte-for-byte, so any change to force-field parameters,
//! integrator constants, or reduction order that silently alters the
//! trajectory breaks this test.
//!
//! Byte-identity holds for a fixed platform and `n_threads` (pinned to 1 in
//! the config below). The reference file was generated on x86_64; a different
//! target's libm (cos/acos/sqrt) can differ by an ULP and would require the
//! reference to be regenerated deliberately, which §13.7 already treats as a
//! physics-change event, not an incidental one.

use std::path::PathBuf;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("geodesic-engine")
        .join("tests")
        .join("fixtures")
}

fn slashed(p: &std::path::Path) -> String {
    p.to_string_lossy().replace('\\', "/")
}

#[test]
fn ala_dipeptide_100_frame_trajectory_matches_golden() {
    let fixtures = fixtures_dir();
    let prmtop = slashed(&fixtures.join("ala_dipeptide.prmtop"));
    let inpcrd = slashed(&fixtures.join("ala_dipeptide.inpcrd"));

    let out_dir = std::env::temp_dir();
    let dcd = slashed(&out_dir.join("geodesic_golden_out.dcd"));
    let csv = slashed(&out_dir.join("geodesic_golden_out.csv"));
    let config_path = out_dir.join("geodesic_golden.toml");

    // 99 steps at frame_interval 1 emits the initial frame plus one per step
    // = 100 frames (§13.7).
    let config = format!(
        r#"
[run]
n_steps        = 99
frame_interval = 1
backend        = "cpu"
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

    let summary = geodesic::run_from_config_file(&config_path).unwrap();
    assert_eq!(summary.n_frames, 100, "golden run must produce exactly 100 frames");

    let produced = std::fs::read(&dcd).unwrap();
    let golden = std::fs::read(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/ala_dipeptide_ref.dcd"),
    )
    .unwrap();

    assert_eq!(
        produced.len(),
        golden.len(),
        "golden trajectory length changed: {} vs {} bytes",
        produced.len(),
        golden.len()
    );
    assert!(
        produced == golden,
        "golden trajectory diverged; regenerate deliberately only if the physics model changed (§13.7)"
    );
}

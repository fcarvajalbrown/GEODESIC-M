# Session Handoff

Last updated: 2026-07-14. **v0.5 shipped** (GPU backend, `geodesic-gpu`,
M2). The non-bonded LJ loop now runs on the GPU via a wgpu compute shader
behind the `gpu` feature; the default CPU-only build is untouched. Next
milestone is v0.6 (hybrid backend, M3) per ROADMAP.md.

## v0.5 (GPU backend) — what shipped

- **New `geodesic-gpu` crate** (`GpuBackend: ComputeBackend`): GPU does
  non-bonded LJ in f32; bonded forces, the constraint solve, and the
  neighbor rebuild stay on CPU (reused from `geodesic-engine`).
- **GPU is f32** (ADR 0002) — WGSL has no f64. The CPU f64 path stays the
  correctness reference and the golden trajectory. GPU/CPU forces agree at
  1e-4 relative on `lj_pair`/`water_box_4`/`ala_dipeptide`.
- **CSR gather list** (ADR 0003): the CPU exclusion-filtered half pair list
  is expanded to a full per-atom CSR (`offsets`/`neighbors`), one thread per
  atom, no `atomicAdd` — deterministic by construction, bit-identical across
  two evaluations on the same adapter.
- **`geodesic-gpu -> geodesic-engine`** crate edge (ADR 0004); the engine
  stays wgpu-free and feature-free.
- **Platform Windows+Linux only** (ADR 0001), DX12/Vulkan, no Metal.
  Resolved and built against **wgpu 22.1.0**. GPU tests are adapter-adaptive
  (skip-with-log when no adapter).
- **`ComputeBackend` gained** `potential_energy`/`atoms`/`topology`/
  `needs_rebuild`/`n_threads` as trait methods (were `CpuBackend` inherent
  methods); the run loop holds `Box<dyn ComputeBackend>` and selects
  CPU/GPU from `run.backend` under the `gpu` feature. `BackendError` gained
  `NoAdapter`.
- **wgpu 22 API note:** `ComputePipelineDescriptor.entry_point` is `&str`
  in wgpu 22 (it became `Option<&str>` in 23) — the one API-drift fixup the
  plan flagged.
- **CI** (`.github/workflows/ci.yml`) is now two jobs, windows + linux; both
  build+clippy `--features gpu` and run the adapter-adaptive GPU tests. The
  byte-exact golden/determinism tests stay windows-only (cross-OS DCD
  byte-identity is not guaranteed), so Linux runs the platform-independent
  crates plus the GPU gate, not the full workspace test.

## Current status

M1 ("CPU-only headless simulation") is done end to end. `geodesic energy`
and `geodesic run` both work against the real ala_dipeptide system, which
runs to completion producing a DCD trajectory + energy CSV with total
energy conserved (~1e-2 kcal/mol over 100 fs at dt=1fs).

- **`geodesic` crate is now lib + bin.** `src/lib.rs` holds
  `energy_from_files` and `run_from_config_file` (the full BAB+RATTLE step
  loop: initial force eval + `constrain_velocities`, then per step
  half_kick -> geodesic_drift -> rebuild-if-stale -> compute_forces ->
  half_kick -> constrain_velocities, writing a DCD frame + CSV row every
  `frame_interval`). `src/main.rs` is just clap parsing + `OPENBLAS_NUM_THREADS=1`.
  Orchestration lives in the binary, not the engine, because it needs
  `geodesic-io` and the §9.3 crate graph forbids the engine from touching I/O.
- **`CpuBackend` gained** `potential_energy()` (accumulated in compute_forces
  from the bonded + per-thread nonbonded energy returns, summed in fixed
  thread-index order), plus `atoms()`/`topology()`/`needs_rebuild()`
  accessors so the driver can compute KE and call `constrain_velocities`
  without a second copy.
- **Golden reference** frozen at `geodesic/tests/golden/ala_dipeptide_ref.dcd`
  (100 frames, dt=1fs, n_threads=1, non-periodic). `golden_reference.rs`
  and `determinism.rs` live in `geodesic/tests/` (not the engine) for the
  crate-graph reason above.
- **Bench suite** (`geodesic/benches/benchmarks.rs`, criterion 0.8): all four
  §13.9 benchmarks. Baselines are hardware-specific, captured per-runner.
- **CI**: `.github/workflows/ci.yml`, test + clippy on **windows-latest**
  (Felipe's choice — matches the golden file's platform so its byte-exact
  check holds). fmt/topo/bench-runner deferred (Felipe confirmed).

## Two real bugs found and fixed this cycle (both were invisible to v0.3's tests)

1. **`half_kick` was missing the force-to-acceleration unit conversion.**
   It did `v += dt/2 · F/m` with no constant. In (Å, amu, kcal/mol, ps)
   units the equation of motion is `a = 20.455² · F/m` (1 kcal/(mol·Å·amu)
   = 418.407 Å/ps²; 20.455 is the AKMA constant already in dcd.rs/inpcrd.rs).
   Without it the dynamics ran ~418× too slow and KE was 418× too small, so
   V (real kcal/mol) and KE were incomparable and NVE energy wasn't
   conserved. **Why v0.3 missed it:** `energy_conservation.rs` and
   `free_rigid_rotor...` both run force-free (epsilon=0, bond constrained)
   and measure only *relative* KE drift, which is invariant to a global
   scale on the force update. It only surfaces comparing absolute V vs KE
   over a real forced run. Constant is `integrator::FORCE_TO_ACCEL_ANG_PER_PS2`.
   The old half_kick unit test had encoded the unitless behavior as its
   expected value — corrected.
2. **Per-atom box wrapping split non-periodic molecules.** `neighbor::build`
   wraps into `[0, box)` (§2.4), but bonded/constraint terms use raw
   (non-min-image) differences. Real ala_dipeptide has negative z coords, so
   wrapping into a 1000 Å box sent an atom ~998 Å from its bonded partners:
   potential exploded to 6.3e8 kcal/mol, constraint solver diverged at step 0.
   Fix (Felipe chose it from options): a `periodic` flag
   (`SimParams.periodic`, config `[system].periodic`, serde default **false**
   = non-periodic) gating the wrap. Non-periodic runs keep the molecule whole
   so raw-difference bonded terms are correct and the trajectory isn't split.
   **Deferred:** minimum-image in the bonded/constraint terms themselves,
   which §2.4 mandates ("all pairwise distances use r*_ij") and is needed for
   a genuinely periodic *bonded* system (tight solvated box). Not needed for
   M1's non-periodic GBSA fixture; do it when a periodic production system is
   actually run (v0.5+).

## Decisions already made this cycle (don't re-ask)

- **PBC handling:** `periodic` flag, default false (non-periodic). Felipe's
  explicit pick over min-image-in-bonded and over molecule-centering.
- **Binary is lib + bin**, run loop in `geodesic/src/lib.rs`; golden +
  determinism tests in `geodesic/tests/`. Forced by the §9.3 crate graph
  (engine can't depend on io).
- **CPU and GPU backends are wired**; `run.backend = "gpu"` selects
  `GpuBackend` under the `gpu` feature (an actionable config error without
  the feature), `"hybrid"` still returns an actionable config error (v0.6).
- **CI on windows-latest, test + clippy only.** fmt deferred (repo is
  hand-formatted, not rustfmt-clean — adopting rustfmt is a separate
  deliberate pass); `--features topo` deferred (no geodesic-topo until v0.8);
  pinned bench runner deferred.
- **Energy subcommand** reports the full non-promoted force-field potential
  at default cutoffs (r_switch=10, r_cutoff=12) in a large non-periodic box,
  printed with the cutoffs so the number is reproducible.

## Next priorities, in order (v0.6, M3 — hybrid backend)

1. **Hybrid backend** (SAD.md §7.4): keep positions/velocities resident
   on-device across steps and minimize host<->device transfers around the
   constraint solve. v0.5 re-uploads positions every `compute_forces` and
   reads forces straight back — correctness-first, not transfer-optimized.
   The GPU-side constraint convergence-reduce (§7.3) also lands here; the
   solve stays on CPU in v0.5.
2. When a genuinely periodic *bonded* system is first run, add minimum-image
   to the bonded forces + constraint solver per §2.4 (see bug 2 above) — the
   one deferred correctness item, and it matters for GPU too since the GPU
   non-bonded already does min-image but the bonded terms do not.
3. The O(N²) neighbor build (`bench_neighbor_rebuild` ~1.25 s at N=10k in
   release) is the obvious perf target; the Verlet list makes the *force*
   loop O(N) but the build itself is still a double loop.

## Things worth knowing that aren't in any file

- The verify-before-trusting pattern caught a real bug again — twice this
  cycle (the units bug and the wrap/split bug), both by actually running the
  CLI end to end and reading the energy CSV rather than trusting that
  green unit tests meant the assembled loop was correct. Relative-drift and
  gradient tests are blind to a global scale error in the integrator; only an
  absolute energy-conservation check on a forced system catches it. Worth
  adding such a check if ever tempted to trust force-free conservation tests.
- Windows path gotcha when hand-writing config.toml from Git Bash: `$(pwd)`
  gives `/c/...` which `PathBuf::join` mangles to `C:/c/...`. Use `pwd -W`
  for real Windows paths, and TOML paths tolerate forward slashes; in Rust
  test-generated configs, `path.replace('\\', "/")` avoids TOML escape issues.
- The golden reference DCD is platform-specific (x86_64 Windows). Cross-OS
  byte-identity of a 100-step MD trajectory is not guaranteed (libm
  cos/acos/sqrt differ by an ULP and propagate). That's why CI is
  windows-latest; regenerate the golden only on a deliberate physics change.

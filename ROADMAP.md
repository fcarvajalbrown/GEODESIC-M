# GEODESIC-M Roadmap

Version milestones from first working build to a stable 1.0. Each version is
scoped to a checkpoint that either compiles and passes its own tests, or ships
a working subset of the CLI — no version is "half a feature." Architecture
reference for every item below: [`docs/SAD.md`](docs/SAD.md). File-level build
order within each M1 phase: [`CLAUDE.md`](CLAUDE.md).

## v0.1 — Core types (`geodesic-core`) — DONE

CLAUDE.md Phase 1. No behavior yet, just the data model everything else is
built on: `SimState`, `AtomData`, `SimParams`, `BondedTopology`,
`NeighborList`, `ForceBuffer`, `TrajectoryFrame`, the `ComputeBackend` trait,
and the full `SimError` hierarchy. Zero internal deps (SAD.md §9.1) — this
crate compiles standalone before anything else exists.

- [x] `Cargo.toml` (workspace root)
- [x] `geodesic-core/Cargo.toml`
- [x] `geodesic-core/src/lib.rs`
- [x] `geodesic-core/src/error.rs`
- [x] `geodesic-core/src/state.rs`
- [x] `geodesic-core/src/atoms.rs`
- [x] `geodesic-core/src/params.rs`
- [x] `geodesic-core/src/topology.rs`
- [x] `geodesic-core/src/buffers.rs`
- [x] `geodesic-core/src/backend.rs`

**Exit criteria:** `cargo check -p geodesic-core` green — met (verified via
`cargo check --workspace` and `cargo clippy --workspace -- -D warnings`,
both clean).

## v0.2 — I/O layer (`geodesic-io`) — DONE

CLAUDE.md Phase 2. Parses `config.toml` → `SimParams`, AMBER `prmtop` →
`AtomData` + `BondedTopology`, AMBER `inpcrd` → initial `SimState`. Writes DCD
trajectory, CSV energy log, JSON barcode, PDB snapshots (SAD.md §10). No
engine dependency — this crate only touches `geodesic-core` types.

- [x] `geodesic-io/Cargo.toml`, `src/lib.rs`
- [x] `src/config.rs` — TOML → `SimParams`, unknown keys rejected
- [x] `src/prmtop.rs`
- [x] `src/inpcrd.rs`
- [x] `src/dcd.rs`
- [x] `src/export.rs` — CSV energy log, JSON barcode
- [x] `src/pdb.rs`

**Exit criteria:** round-trip tests per SAD.md §13.8 pass — met, 9 tests
across all six files, all green. One deviation from §13.1/§13.8's literal
wording: those sections name `ala_dipeptide.prmtop` as the round-trip
fixture, but that's a real 22-atom force-field file that can't be
hand-typed without risking fabricated parameters — a hand-built,
hand-verified `lj_pair` fixture (2 atoms, self-consistent LJ/bond/angle/
dihedral/exclusion data) was used instead to test parser *mechanics*.
A real AmberTools-generated `ala_dipeptide.prmtop` is still needed before
v0.3's physics tests (§13.1) can check force-field *values* against
literature, not just parser correctness.

## v0.3 — Force field + CPU integrator (`geodesic-engine`) — IN PROGRESS, PAUSED

CLAUDE.md Phase 3. The physics: Verlet neighbor lists, bonded (bond/angle/
dihedral) and non-bonded (LJ) forces in SoA, the Lagrangian constraint solver,
and the Geodesic BAB integrator (SAD.md §2.3, §7.2). This is the crate where
correctness bugs are most expensive, so it ships with its test suite, not
after it.

- [x] `src/neighbor.rs` — Verlet list, PBC wrap, min-image, bonded exclusions
- [x] `src/force/nonbonded.rs` — LJ + quintic switching function
- [x] `src/force/bonded.rs` — bond and angle forces (dihedral **known broken**,
      see the doc comment on `compute_dihedral_forces` — f_i/f_l correct,
      f_j/f_k formula is structurally incomplete for general geometry)
- [ ] `src/constraint.rs`, `src/integrator.rs`, `src/cpu_backend.rs` — not
      started (placeholder files only, so the workspace still compiles)
- [ ] Fixtures: `lj_pair`, `harmonic_dimer`, `water_box_4`, `ala_dipeptide` —
      not started; ad-hoc unit-test fixtures exist inline in
      `tests/neighbor_list.rs`, `tests/nonbonded_gradient.rs`,
      `tests/bonded_gradient.rs` but don't match SAD.md §13.1's named files
- [ ] `tests/gradient_check.rs` (§13.2) — the *pattern* exists (per-file
      finite-difference tests above), not yet consolidated into this
      SAD.md-named file
- [ ] `tests/newton_third_law.rs` (§13.3) — nonbonded has this inline
      (`newtons_third_law_holds`); not yet a standalone file, not yet
      covering bonded forces
- [ ] `tests/energy_conservation.rs` (§13.4) — blocked on the integrator
- [ ] `tests/determinism.rs` (§13.5) — blocked on the CPU backend
- [ ] `tests/constraint_solver.rs` (§13.6) — blocked on the constraint solver

**Session paused here 2026-07-10.** Full handoff, including the dihedral
bug diagnosis and a concrete lead on the fix, is in `memory.md` at the repo
root. Workspace is verified green (`cargo check --workspace`, `cargo clippy
--workspace --all-targets -- -D warnings`, `cargo test --workspace` all
pass) — the two known-broken dihedral tests are `#[ignore]`d with reasons,
not silently failing.

**Exit criteria:** all five test files above pass on all four fixtures. Not
yet met.

## v0.4 — CLI binary — M1 complete

CLAUDE.md Phase 4. `geodesic energy <prmtop> <inpcrd>` and `geodesic run
<config.toml>` (SAD.md §9.2). Golden reference trajectory frozen
(`ala_dipeptide_ref.dcd`, §13.7) once this build is verified correct — this
is the point M1 ("CPU-only headless simulation") is actually done, not just
"compiles."

- [ ] `geodesic/Cargo.toml`, `src/main.rs`
- [ ] `tests/golden_reference.rs`
- [ ] `cargo bench` baselines committed (§13.9: `bench_lj_inner_loop`,
      `bench_neighbor_rebuild`, `bench_constraint_solver`, `bench_full_step`)
- [ ] CI matrix green (§13.10): test, clippy, fmt

**Exit criteria:** a real `.prmtop`/`.inpcrd` pair runs end to end and
produces a DCD + energy CSV; golden reference test passes; benchmarks have a
committed baseline.

## v0.5 — GPU backend (`geodesic-gpu`, M2)

wgpu compute shaders for the non-bonded force loop (SAD.md §7.3): tiled
force evaluation, `GpuBackend: ComputeBackend`, fixed-order tree reduction
for determinism, `DeviceLost`/`ShaderCompilation`/`OutOfGpuMemory` mapped to
`BackendError` (§12.6). Gated behind the `gpu` feature — default CPU-only
build is unaffected.

**Exit criteria:** `cargo build --features gpu` succeeds; GPU and CPU
backends agree on forces for the fixture systems within floating-point
tolerance; same-GPU/same-driver runs are bit-for-bit reproducible.

## v0.6 — Hybrid backend (M3)

CPU/GPU workload split per SAD.md §7.4: non-bonded on GPU, bonded forces +
constraint solve + neighbor rebuild on CPU. Position/velocity handoff
between devices around the constraint solve step.

**Exit criteria:** hybrid backend matches pure-CPU and pure-GPU energy
trajectories on the fixture suite within the same tolerance used in v0.5.

## v0.7 — GUI renderer scaffolding (M4, part 1)

`geodesic-gui` crate: wgpu 3D viewer, ring buffer consumer, atoms as
instanced spheres, bonds as cylinders (SAD.md §7.5). Runs on a dedicated OS
thread, decoupled from the simulation loop — the sim never blocks on the
renderer.

**Exit criteria:** a completed `geodesic run` trajectory can be replayed in
the viewer; live streaming not required yet (that's v0.8).

## v0.8 — Live streaming + topology pipeline (M4, part 2)

Live `TrajectoryFrame` streaming from a running simulation into the GUI via
the ring buffer. Adds `geodesic-topo` (PSL flexibility analysis §2.6, Zigzag
persistence §2.7, Ripser FFI per §11.4) so the GUI can color atoms by
flexibility score, not just element.

**Exit criteria:** GUI shows a running simulation live; PSL scores backfill
into the B-factor column of a PDB snapshot; barcode JSON matches the schema
in SAD.md §10.6.

## v0.9 — Data export UI (M5)

Export panel in the GUI: PDB snapshot, DCD trajectory segment, energy/RMSD
CSV, Zigzag barcode JSON — all on demand, per SAD.md §7.5 "Data export."

**Exit criteria:** every export format the CLI can produce is also
reachable from the GUI without dropping back to the command line.

## v1.0 — Stable release

No new features — a hardening pass against the NFRs in SAD.md §5:

- [ ] Determinism re-verified across CPU, GPU, and hybrid backends
- [ ] Auditability: every force contribution traceable to atom pair + step
- [ ] Performance targets benchmarked and documented (ns/day throughput)
- [ ] Full CI matrix (§13.10) green, including `--features topo`
- [ ] Golden reference trajectory unchanged since v0.4 (or regenerated with
      an explicit, documented physics-change justification)
- [ ] `README.md` usage section covers all five CLI/GUI export paths

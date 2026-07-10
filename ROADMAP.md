# GEODESIC-M Roadmap

Version milestones from first working build to a stable 1.0. Each version is
scoped to a checkpoint that either compiles and passes its own tests, or ships
a working subset of the CLI — no version is "half a feature." Architecture
reference for every item below: [`docs/SAD.md`](docs/SAD.md). File-level build
order within each M1 phase: [`CLAUDE.md`](CLAUDE.md).

## v0.1 — Core types (`geodesic-core`) — IN PROGRESS

CLAUDE.md Phase 1. No behavior yet, just the data model everything else is
built on: `SimState`, `AtomData`, `SimParams`, `BondedTopology`,
`NeighborList`, `ForceBuffer`, `TrajectoryFrame`, the `ComputeBackend` trait,
and the full `SimError` hierarchy. Zero internal deps (SAD.md §9.1) — this
crate compiles standalone before anything else exists.

- [x] `Cargo.toml` (workspace root)
- [x] `geodesic-core/Cargo.toml`
- [x] `geodesic-core/src/lib.rs`
- [ ] `geodesic-core/src/error.rs`
- [ ] `geodesic-core/src/state.rs`
- [ ] `geodesic-core/src/atoms.rs`
- [ ] `geodesic-core/src/params.rs`
- [ ] `geodesic-core/src/topology.rs`
- [ ] `geodesic-core/src/buffers.rs`
- [ ] `geodesic-core/src/backend.rs`

**Exit criteria:** `cargo check -p geodesic-core` green.

## v0.2 — I/O layer (`geodesic-io`)

CLAUDE.md Phase 2. Parses `config.toml` → `SimParams`, AMBER `prmtop` →
`AtomData` + `BondedTopology`, AMBER `inpcrd` → initial `SimState`. Writes DCD
trajectory, CSV energy log, JSON barcode, PDB snapshots (SAD.md §10). No
engine dependency — this crate only touches `geodesic-core` types.

- [ ] `geodesic-io/Cargo.toml`, `src/lib.rs`
- [ ] `src/config.rs` — TOML → `SimParams`, unknown keys rejected
- [ ] `src/prmtop.rs`
- [ ] `src/inpcrd.rs`
- [ ] `src/dcd.rs`
- [ ] `src/export.rs` — CSV energy log, JSON barcode
- [ ] `src/pdb.rs`

**Exit criteria:** round-trip tests per SAD.md §13.8 pass (prmtop atom/bond
counts, inpcrd position precision, DCD frame count, TOML rejection rules).

## v0.3 — Force field + CPU integrator (`geodesic-engine`)

CLAUDE.md Phase 3. The physics: Verlet neighbor lists, bonded (bond/angle/
dihedral) and non-bonded (LJ) forces in SoA, the Lagrangian constraint solver,
and the Geodesic BAB integrator (SAD.md §2.3, §7.2). This is the crate where
correctness bugs are most expensive, so it ships with its test suite, not
after it.

- [ ] `src/neighbor.rs`, `src/force/{mod,nonbonded,bonded}.rs`
- [ ] `src/constraint.rs`, `src/integrator.rs`, `src/cpu_backend.rs`
- [ ] Fixtures: `lj_pair`, `harmonic_dimer`, `water_box_4`, `ala_dipeptide`
- [ ] `tests/gradient_check.rs` (finite-difference force check, §13.2)
- [ ] `tests/newton_third_law.rs` (§13.3)
- [ ] `tests/energy_conservation.rs` (§13.4)
- [ ] `tests/determinism.rs` (§13.5)
- [ ] `tests/constraint_solver.rs` (§13.6)

**Exit criteria:** all five test files above pass on all four fixtures.

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

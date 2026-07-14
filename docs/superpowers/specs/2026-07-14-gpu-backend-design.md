# v0.5 — GPU Backend (`geodesic-gpu`, M2) — Design

Date: 2026-07-14
Author: Felipe Carvajal Brown
Status: Approved (brainstorming), pending implementation plan
Architecture reference: [`docs/SAD.md`](../../SAD.md) §7.3, §8, §9, §12.6; [`ROADMAP.md`](../../../ROADMAP.md) v0.5

## Goal

Add a `GpuBackend: ComputeBackend` that evaluates the non-bonded
Lennard-Jones force loop on the GPU via wgpu compute shaders, gated behind the
`gpu` Cargo feature. The default CPU-only build is unaffected. This is
milestone M2 — correctness first, not peak performance (transfer/residency
optimization is M3/v0.6).

## Scope

On GPU: the non-bonded LJ force loop only.

On CPU (reused from `geodesic-engine`, unchanged): bonded forces
(bond/angle/dihedral), the Lagrangian constraint solve (SHAKE position +
RATTLE velocity), and the Verlet neighbor-list rebuild.

Rationale: SAD §7.3 maps only the O(N^2) non-bonded loop to the GPU; §7.4
(hybrid, v0.6) is where the CPU/GPU split is formalized and optimized. Porting
bonded forces or the constraint solve to GPU is explicitly out of scope — SAD
argues both are poor GPU fits (low N, dependency graph, iterative
convergence).

Out of scope for v0.5, deferred to v0.6 (M3):
- Keeping positions/velocities resident on-device across steps.
- Minimizing host<->device transfers around the constraint solve.
- The constraint convergence-reduce on GPU (SAD §7.3 mentions it; not required
  to meet M2 exit criteria, and the solve stays on CPU here).

## Data flow — `compute_forces`

1. CPU rebuilds the exclusion-filtered Verlet list (`engine::neighbor::build`)
   exactly as today. Bonded exclusions (1-2, 1-3) and the cutoff are already
   applied in that list, so the GPU never has to re-derive them.
2. On each neighbor rebuild, the CPU expands the half pair-list (`i < j`) into
   a full per-atom neighbor list in CSR form: `offsets: [u32; N+1]` and a flat
   `neighbors: [u32]`. This, plus per-atom `sigma`, `epsilon`, `charge`, is
   uploaded once per rebuild. Positions upload every step, as `f32`.
3. Kernel: one GPU thread per atom `i`. The thread gathers over its own
   neighbor slice `neighbors[offsets[i]..offsets[i+1]]`, accumulating `F_i` and
   a per-atom energy contribution in fixed neighbor-index order. No scatter, no
   `atomicAdd`. Newton's third law holds because atom `j` independently carries
   `i` in its own slice and computes the equal-and-opposite term.
4. Bonded forces run on CPU and are summed into the same `ForceBuffer` the GPU
   forces are read back into. Total non-bonded energy is a fixed-order tree
   reduction over the per-atom energy contributions (halved to correct the
   double count).

Why gather-per-atom rather than a flat pair-list scatter: a flat pair-list
would need each pair to write to two atoms, forcing either `atomicAdd`
(nondeterministic reduction order, forbidden by SAD §7.3) or a separate scatter
reduction. The full per-atom (gather) layout makes each atom's force the
private result of a single thread — deterministic by construction. This is the
reconciliation of "reuse the CPU-derived pair set" (chosen) with the SAD
determinism requirement; it deviates from SAD §7.3's literal "tiled all-pairs"
wording (see ADR below).

## Precision

WGSL has no `f64` type; wgpu compute is `f32`-only. Therefore:
- GPU non-bonded forces and energies are computed in `f32`.
- The CPU `f64` path remains the reference for correctness, the golden
  reference trajectory, and the `geodesic energy` subcommand. The GPU backend
  is an opt-in, documented lower-precision path.
- The GPU/CPU agreement test compares against CPU `f64` at an `f32`-appropriate
  tolerance (target ~1e-4 relative on force components; the exact figure is
  fixed empirically from the three fixtures at implementation time and recorded
  in the test).

This is a documented deviation from SAD §8's `f64`-everywhere mandate (see ADR
below). Double-single (two-`f32`) emulation was considered and rejected for M2:
large scope, ~4-10x slowdown, error-prone, and unnecessary for a
correctness-first first cut.

## Determinism & error handling

- Same GPU model + same driver version -> bit-for-bit reproducible `f32`
  output. Achieved via the no-atomics gather kernel and a fixed-order tree
  reduction for the energy sum. Cross-GPU reproducibility is not guaranteed and
  is documented as a known limitation (consistent with SAD §7.3).
- wgpu error mapping to the existing `BackendError` variants (SAD §12.6), all
  hard-stop, no silent CPU fallback (a fallback would silently change results
  and violate the Determinism NFR):
  - `DeviceLost` -> hard stop.
  - `ShaderCompilation(String)` -> hard stop at startup, before any step runs
    (the WGSL is vendored; a compile failure is a bug, not recoverable).
  - `OutOfGpuMemory` -> hard stop with a message suggesting a smaller N or the
    CPU backend.

## Platform, adapter selection & CI

- Target platforms: Windows + Linux only. wgpu uses the DX12 and Vulkan
  backends; macOS/Metal is dropped (see ADR below).
- Tests are adapter-adaptive: request any available wgpu adapter — a real GPU
  locally, the DX12 WARP software adapter on Windows CI, or Vulkan lavapipe on
  Linux CI. If no adapter can be created, the GPU tests skip with a logged
  reason rather than fail.
- CI hard-checks `cargo build --features gpu` and `cargo clippy --features gpu`
  on both Windows and Linux runners. The agreement and reproducibility tests
  run wherever an adapter (real or software) exists.

## Crate structure

- New crate `geodesic-gpu/` (feature `gpu`), holding `GpuBackend`, wgpu device
  management, and the vendored WGSL non-bonded compute shader.
- Dependency edges: `geodesic-gpu -> geodesic-engine -> geodesic-core`, and
  `geodesic-gpu -> {wgpu, bytemuck, pollster}`. The new
  `geodesic-gpu -> geodesic-engine` edge (not in SAD §9.3) lets `GpuBackend`
  delegate the neighbor/bonded/constraint work to existing engine code while
  keeping `geodesic-engine` itself entirely wgpu-free and feature-free — a
  default M1 build compiles no GPU code at all. This crate-graph edge is a
  documented deviation (see ADR below).
- Engine modules `GpuBackend` needs to call (`neighbor`, `force::bonded`,
  `constraint`, `integrator`) must be made `pub` from `geodesic-engine` if they
  are not already. This is a visibility change only, no behavior change.
- The `geodesic` binary already selects the backend at startup and currently
  returns an actionable config error for `run.backend = "gpu"`. v0.5 wires
  `GpuBackend` in behind the `gpu` feature so that error path is replaced by a
  real backend when the feature is enabled.

Dependency versions per SAD §14: `wgpu` (SAD pins 22; confirm or bump to a
current release at implementation, since 22 predates this work), `bytemuck` 1,
`pollster` 0.3. Verify current versions against crates.io during
implementation.

## Testing

- `geodesic-gpu/tests/forces_gpu_vs_cpu.rs` — GPU vs CPU non-bonded forces on
  `lj_pair`, `water_box_4`, and `ala_dipeptide`, asserting agreement at the
  `f32` tolerance. `ala_dipeptide` is the exclusion-heavy real system that
  validates the CSR neighbor-list expansion.
- `geodesic-gpu/tests/gpu_determinism.rs` — two GPU force evaluations on the
  same adapter produce bit-identical `f32` output.
- Both tests are adapter-adaptive (skip-with-log if no adapter).
- CI build-only gate ensures `--features gpu` keeps compiling on both OSes.

## Exit criteria (from ROADMAP v0.5)

- `cargo build --features gpu` succeeds (Windows + Linux).
- GPU and CPU backends agree on non-bonded forces for the fixture systems
  within the `f32` tolerance fixed above.
- Same-GPU/same-driver runs are bit-for-bit reproducible (`f32`).

## ADRs produced by this milestone

To be written as the first implementation step, before code:
1. Windows + Linux only as the GPU target platform (drop Metal).
2. GPU backend is `f32` — documented deviation from SAD §8's `f64` mandate; CPU
   `f64` remains the correctness reference.
3. GPU non-bonded evaluates the CPU-derived neighbor set as a per-atom gather
   (CSR) list — deviation from SAD §7.3's literal "tiled all-pairs."
4. New `geodesic-gpu -> geodesic-engine` crate-graph edge — deviation from SAD
   §9.3; keeps the engine wgpu-free.

## Open items / risks

- The `f32` agreement tolerance is empirical; if `ala_dipeptide` forces don't
  agree within a defensible `f32` bound, that is a real bug (likely in the CSR
  expansion or exclusion handling), not something to loosen the tolerance
  around — fix at the root per project rules.
- The half-to-full neighbor-list expansion runs on CPU each rebuild; its cost
  is O(pairs) and negligible next to the O(N^2) build itself, but it is new
  code and needs its own unit test (round-trip: expanded full list reproduces
  the original pair set).

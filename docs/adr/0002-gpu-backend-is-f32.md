# 0002 — GPU backend computes non-bonded LJ in f32

Status: Accepted
Date: 2026-07-14
Deciders: Felipe Carvajal Brown

## Context

SAD §8 mandates `f64` for all simulation arithmetic. WGSL has no `f64` type and wgpu compute is `f32`-only, so the non-bonded loop moved to a compute shader cannot honor that mandate directly.

## Decision

The GPU backend computes non-bonded Lennard-Jones forces and energies in `f32`. The CPU `f64` path remains the correctness reference: it produces the golden-reference trajectory and backs the `geodesic energy` subcommand. The GPU backend is an opt-in, documented lower-precision path. GPU/CPU agreement is asserted at an `f32`-appropriate tolerance (target ~1e-4 relative on force components), fixed empirically from the fixtures and recorded in the tests.

## Consequences

- Selecting `run.backend = "gpu"` accepts `f32` non-bonded arithmetic; bonded forces and the constraint solve stay on CPU in `f64`.
- Cross-GPU bit-reproducibility is not guaranteed; same-GPU/same-driver runs are bit-identical (see ADR 0003).
- A dedicated GPU/CPU agreement test guards the `f32` path; if it disagrees beyond the tolerance the fault is a real bug, not something to loosen the tolerance around.

## Alternatives considered

- **Double-single (two-`f32`) emulation.** Rejected for M2: large scope, ~4-10x slowdown, error-prone, unnecessary for a correctness-first first cut.
- **CUDA `f64`.** Rejected: NVIDIA-only, contradicts the SAD's wgpu choice and the planned GUI-sharing of the wgpu device.

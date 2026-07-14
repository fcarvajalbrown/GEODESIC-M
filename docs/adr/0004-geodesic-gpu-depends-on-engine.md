# 0004 — New `geodesic-gpu -> geodesic-engine` crate-graph edge

Status: Accepted
Date: 2026-07-14
Deciders: Felipe Carvajal Brown

## Context

SAD §9.3 fixes the crate graph. The GPU backend only moves the non-bonded loop to the GPU; it reuses the engine's CPU code for bonded forces, the constraint solve, and the neighbor rebuild. It needs to call `geodesic_engine::{neighbor, force::bonded, integrator}`. Keeping the engine wgpu-free (so the default M1 build compiles no GPU code) is a hard requirement.

## Decision

Add a `geodesic-gpu -> geodesic-engine` dependency edge (not present in SAD §9.3). `GpuBackend` lives in `geodesic-gpu` and delegates neighbor/bonded/constraint/integrator work to existing engine code. `geodesic-engine` stays entirely wgpu-free and feature-free; the `gpu` feature lives on the `geodesic` binary and pulls in `geodesic-gpu` as an optional dependency. Engine modules the backend calls are made `pub` (visibility only, no behavior change).

## Consequences

- The default CPU-only build (`geodesic-engine`, `geodesic` without `--features gpu`) compiles no wgpu code at all.
- `geodesic-gpu` reuses engine physics without duplicating it, so CPU and GPU share one bonded/constraint/neighbor implementation.
- The crate graph gains one edge beyond SAD §9.3; documented here.

## Alternatives considered

- **Put `GpuBackend` inside `geodesic-engine` behind a feature.** Rejected: forces wgpu into the engine's dependency tree and complicates the feature-free engine invariant.
- **Duplicate the needed engine code into `geodesic-gpu`.** Rejected: two copies of bonded/constraint/neighbor logic drift apart and break the single-reference guarantee.

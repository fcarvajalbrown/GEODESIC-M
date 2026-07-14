# 0005 — Hybrid backend is a documented alias for the GPU backend

Status: Accepted
Date: 2026-07-14
Deciders: Felipe Carvajal Brown

## Context

SAD.md §7.3 (M2, GPU) and §7.4 (M3, hybrid) read as two milestones producing two backends. In practice the v0.5 `GpuBackend` already implements the §7.4 workload split: non-bonded LJ on the GPU, bonded forces + constraint solve + neighbor rebuild on the CPU in f64. There is no second, distinct hybrid engine to build — only transfer efficiency to add to the one that exists. The `config.toml` `backend` field already accepts `"hybrid"`, which v0.5 rejects with an actionable "lands in v0.6" error.

## Decision

`backend = "hybrid"` is a documented alias for `backend = "gpu"`. Under `--features gpu` both construct the same optimized `GpuBackend`; without the feature both return the same actionable "rebuild with --features gpu" error keyed to `run.backend`. GEODESIC-M ships one GPU/hybrid backend, not two. The ROADMAP v0.6 exit criterion "hybrid matches pure-CPU and pure-GPU" collapses to "CPU ~= GPU" because hybrid and gpu are the same code path.

## Consequences

- Users can select either name; behavior and output are bit-identical between them on the same adapter.
- One backend implementation to test and maintain, not two.
- The SAD's M2/M3 milestone framing is documented here as one engine delivered in two increments (v0.5 correctness, v0.6 transfer efficiency), not two engines.

## Alternatives considered

- A separate `HybridBackend` struct. Rejected: it would duplicate the v0.5 backend with no behavioral difference, since v0.5 is already the §7.4 split.
- Keep rejecting `"hybrid"`. Rejected: the workload split it names exists and ships; rejecting the name is misleading.

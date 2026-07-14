# 0006 — GPU constraint convergence-reduce dropped

Status: Accepted
Date: 2026-07-14
Deciders: Felipe Carvajal Brown

## Context

SAD.md §7.3 specifies a GPU tree reduction for the constraint solver's convergence check (global max|lambda_i|), "no CPU round-trip per iteration." SAD.md §7.4, the hybrid split, keeps the entire constraint solve on the CPU (iterative, convergence-dependent, CPU overhead acceptable). These two sections contradict each other: a GPU convergence-reduce only makes sense if the lambdas live on the GPU, but §7.4 computes them on the CPU.

## Decision

The GPU constraint convergence-reduce is dropped. With the solve on the CPU (SHAKE/RATTLE in `integrator::drift_and_constrain`), the lambdas and their residuals are already in host memory; the convergence check is a CPU reduction, as it is today. Shipping residuals to the GPU to reduce a single scalar and read it back would be a net-negative transfer.

## Consequences

- Resolves the §7.3 / §7.4 contradiction in favor of §7.4 (solve on CPU).
- No GPU code is written for the constraint solver in v0.6; the convergence check stays the deterministic CPU reduction already in the engine.
- If a future milestone moves the whole constraint solve onto the GPU, a GPU-side convergence-reduce can be reconsidered — it only pays off then.

## Alternatives considered

- Implement the §7.3 GPU reduce as written. Rejected: net-negative transfer while the solve is on the CPU; adds a round-trip to read one scalar.
- Move the whole constraint solve to the GPU. Rejected for v0.6: §7.4 forbids it (iterative, convergence-dependent, poor GPU fit), and it is a far larger scope than a transfer-optimization milestone.

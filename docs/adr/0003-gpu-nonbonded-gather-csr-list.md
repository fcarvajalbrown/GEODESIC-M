# 0003 — GPU non-bonded evaluates the CPU pair set as a per-atom CSR gather list

Status: Accepted
Date: 2026-07-14
Deciders: Felipe Carvajal Brown

## Context

SAD §7.3 describes the GPU non-bonded loop as "tiled all-pairs" and requires deterministic reduction (no nondeterministic `atomicAdd` order). The CPU already builds an exclusion-filtered Verlet half pair-list (1-2 and 1-3 exclusions and the cutoff applied). Re-deriving exclusions on the GPU would duplicate that logic and risk divergence; a flat pair-list scatter would force `atomicAdd` (nondeterministic order) or a separate scatter reduction.

## Decision

On each neighbor rebuild the CPU expands the half pair-list (`i < j`) into a full per-atom neighbor list in CSR form (`offsets: [u32; N+1]`, flat `neighbors: [u32]`). The kernel runs one GPU thread per atom `i`, gathering over its own neighbor slice in fixed index order and accumulating `F_i` and a per-atom energy privately. No scatter, no `atomicAdd`. Newton's third law holds because atom `j` independently carries `i` in its own slice and computes the equal-and-opposite term. Total non-bonded energy is a fixed-order reduction over per-atom contributions, halved to correct the double count.

## Consequences

- Each atom's force is the private result of a single thread — deterministic by construction; same-GPU/same-driver output is bit-identical.
- The CPU exclusion set is reused exactly; the GPU never re-derives exclusions.
- The half-to-full expansion is O(pairs) CPU work per rebuild, negligible next to the O(N^2) build, and has its own round-trip unit test.
- This deviates from SAD §7.3's literal "tiled all-pairs" wording; the deterministic-gather intent of §7.3 is preserved.

## Alternatives considered

- **Flat pair-list scatter with `atomicAdd`.** Rejected: nondeterministic reduction order, forbidden by SAD §7.3.
- **Tiled all-pairs ignoring the CPU exclusion set.** Rejected: duplicates exclusion logic on GPU and risks divergence from the CPU reference.

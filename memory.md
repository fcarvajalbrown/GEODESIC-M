# Session Handoff

Last updated: 2026-07-10, mid-session continuation (same day as the prior
handoff below). v0.3's core physics is now complete; two items remain
before the v0.3 checkpoint/release, both blocked on Felipe's input rather
than further coding.

## Current status

`geodesic-engine` (v0.3) core physics is done and gradient/physics-tested:

- `neighbor.rs`, `force/nonbonded.rs` — unchanged from before, still passing.
- `force/bonded.rs` — **dihedral forces fixed.** The f_j/f_k projection
  ansatz had a sign error on the `p = b1·b2/|b2|^2` term: it used `(p-1)`
  instead of `-(1+p)` for f_j's f_i-coefficient, and `-p` instead of `+p` for
  f_k's f_i-coefficient. Found via full symbolic (sympy) chain-rule
  differentiation of `phi = atan2((m×n)·b2/|b2|, m·n)`, verified exact
  (machine epsilon) against three independent geometries before touching
  Rust. Both previously-`#[ignore]`d gradient tests are re-enabled and pass.
- `constraint.rs` — Lagrangian solver, Jacobi-relaxation (not Gauss-Seidel —
  needed for determinism once dispatch is parallelized, SAD.md §7.2).
  `solve()` is the position (SHAKE) stage; `constrain_velocities()` is
  RATTLE's separate velocity-tangency stage, added after verifying against
  the primary source (Leimkuhler & Matthews 2016, PMC4893190) that SHAKE's
  position correction alone is NOT sufficient — a force half-kick can
  reintroduce an along-bond velocity component that only a second,
  independent projection removes. `promote_hydrogen_bonds()` moves X-H
  bonds from the harmonic bond list into rigid constraints; gated by a new
  `constrain_hydrogens` config field (default true) — decided with Felipe
  as a config-level policy toggle, matching AMBER's `ntc=2` SHAKE
  convention (a runtime choice, not something prmtop encodes).
- `integrator.rs` — `half_kick` (B) and `drift_and_constrain` (A: drift +
  `constraint::solve` + velocity resync from the actual constrained
  displacement). Test parameters for the free-rigid-rotor energy
  conservation check were cross-validated against an independent Python
  implementation before being encoded as Rust assertions (2000 steps,
  dt=0.004 ps, KE and bond length both hold to ~1e-9).
- `cpu_backend.rs` — `CpuBackend` implements `ComputeBackend`: static strip
  decomposition of the flat Verlet pair list across a fixed Rayon thread
  pool, private per-thread `ForceBuffer`s, sequential thread-index-ordered
  reduction. Bonded forces stay single-threaded (low N vs. the pair count).
  `compute_pair_forces` was narrowed to take `pair_i`/`pair_j`/`r_cutoff`/
  `r_switch` directly instead of a whole `&NeighborList`, since strip-slicing
  a borrowed list wasn't otherwise possible without a copy.
  **Determinism finding, worth remembering:** same-thread-count runs are
  bit-for-bit repeatable (the real SAD.md §7.2 claim), but *different*
  thread counts are NOT bit-identical to each other — summing one atom's
  pair contributions via several per-thread partial sums isn't associative
  with summing them as one running total. A same-T-repeatability test and a
  looser cross-T tolerance test are both in `tests/cpu_backend.rs` so this
  doesn't get mis-asserted as a bug again later.

**Verified clean:** `cargo check --workspace`, `cargo clippy --workspace
--all-targets -- -D warnings`, `cargo test --workspace` all pass, zero
warnings, zero failures, zero ignored tests (the two dihedral ignores from
before are gone).

## Next priorities, in order

1. **Real physics fixtures (SAD.md §13.1):** `lj_pair`, `harmonic_dimer`,
   `water_box_4`, `ala_dipeptide` under `geodesic-engine/tests/fixtures/`.
   `ala_dipeptide` is a real 22-atom system — its prmtop force-field values
   cannot be hand-typed (CLAUDE.md: never invent physics constants). Needs
   Felipe to either generate one via AmberTools or point to a citable
   public-domain source. **Blocked on Felipe, asked but not yet answered
   as of this handoff.**
2. **Test file consolidation:** ad-hoc files (`neighbor_list.rs`,
   `nonbonded_gradient.rs`, `bonded_gradient.rs`, `constraint_solver.rs`,
   `hydrogen_constraint_promotion.rs`, `integrator.rs`, `cpu_backend.rs`)
   vs. the five SAD.md-named files (`gradient_check.rs`,
   `newton_third_law.rs`, `energy_conservation.rs`, `determinism.rs`,
   `constraint_solver.rs` — this last one already matches). Genuinely
   unclear whether consolidating is worth the busywork vs. leaving the
   per-file names, which are arguably more readable. **Ask Felipe.**
3. v0.3 checkpoint once 1–2 are resolved: update ROADMAP.md, tag and
   publish the `v0.3` GitHub release (voice-checked per CLAUDE.md).
4. v0.4: `geodesic/src/main.rs` — clap CLI, `energy` + `run` subcommands
   (SAD.md §9.2), set `OPENBLAS_NUM_THREADS=1` at startup (§11.3, harmless
   now, correct later since topology pipeline isn't built). This is also
   where the full BAB+RATTLE step sequencing gets assembled: `half_kick` →
   `ComputeBackend::geodesic_drift` → (rebuild neighbor list if
   `neighbor::needs_rebuild`) → `ComputeBackend::compute_forces` →
   `half_kick` → `constraint::constrain_velocities` (this last call is
   deliberately NOT part of `ComputeBackend` — see `cpu_backend.rs`'s doc
   comment on `CpuBackend` for why). Then `golden_reference.rs` (§13.7,
   frozen trajectory) and `cargo bench` baselines (§13.9). v0.4 checkpoint:
   ROADMAP.md update, tag + publish `v0.4`.

## Pending deferred items

- **`ala_dipeptide.prmtop`**: still needed, still can't be hand-typed. See
  priority 1 above.
- **PSL/Zigzag persistence, GPU backend, GUI**: v0.5+ per ROADMAP.md,
  correctly out of scope for now.

## Things worth knowing that aren't in any file

- The finite-difference / symbolic-verification-before-trusting-hand-derived-math
  pattern keeps paying off: this session's dihedral fix and the two-stage
  RATTLE design were both things a plausible-looking hand derivation would
  have gotten subtly wrong (dihedral: a structurally incomplete projection
  ansatz that passed one symmetric test case; RATTLE: SHAKE's position
  correction alone looks sufficient until you check the primary source and
  see there's a genuinely separate velocity-tangency stage). Keep doing
  this for anything with a non-obvious formula — it has caught a real bug
  or prevented one every single time it's been done in this project.
- When SAD.md's stated architecture doesn't map 1:1 onto a concrete
  implementation choice (e.g. "N×N pair interaction space partitioned into
  T strips" from §7.2 assumes a naive O(N²) double loop, but the actual
  code already has a flat Verlet pair list from `neighbor.rs`), the right
  move is to preserve the *intent* (static partition, private buffers,
  fixed-order reduction) rather than the literal mechanism, and document
  the adaptation inline — not to ask Felipe about every such translation,
  since these are implementation details that follow directly from
  already-decided architecture.

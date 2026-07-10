# Session Handoff

Last updated: 2026-07-10, stopped mid v0.3 by explicit instruction (token
budget). This file exists per CLAUDE.md's "Session Handoff" rule.

## Current status

Repo revived from zero to v0.3-in-progress in one session. v0.1
(`geodesic-core`, data model) and v0.2 (`geodesic-io`, file I/O) are done,
tested, tagged, and released (`v0.1`, `v0.2` on GitHub). v0.3
(`geodesic-engine`, the physics) is roughly half-done:

**Working and gradient-tested:**
- `geodesic-engine/src/neighbor.rs` — Verlet neighbor list, PBC wrap,
  minimum-image, bonded-pair exclusion. 5 tests passing.
- `geodesic-engine/src/force/nonbonded.rs` — Lennard-Jones with a quintic
  switching function. 5 tests passing, including finite-difference gradient
  checks in both the unswitched and switched regions.
- `geodesic-engine/src/force/bonded.rs` — bond and angle forces only.
  3 of 5 tests passing (2 dihedral tests intentionally `#[ignore]`d).

**Broken, do not use:** `compute_dihedral_forces` in `bonded.rs`. The outer
forces (f_i, f_l) are verified correct. The inner forces (f_j, f_k) use a
projection ansatz derived from f_i/f_l that's wrong for general geometry —
passes a symmetric test case, fails an asymmetric one with a non-sign-flip
error (~8% magnitude, not a clean negation like the bond/angle bugs were).
Full diagnosis and numeric ground-truth values are in the doc comment on
the function itself (`geodesic-engine/src/force/bonded.rs`). The likely
missing piece: the projection ansatz only captures phi's dependence on b2
*indirectly* through m=b1×b2 and n=b2×b3, but phi = atan2(y,x) where
y = (m×n)·b2/|b2| also depends on b2 *directly* (both in the dot product
and the |b2| normalization) — that direct term is probably what's missing.
A full symbolic chain-rule derivation of dphi/drj through all three paths
would settle it; the partial numeric derivation (Python, via finite
difference on phi itself rather than on V) confirmed dphi/dri and dphi/drl
exactly but didn't get far enough to fully pin down dphi/drj before the
session was paused.

**Not started:** `constraint.rs`, `integrator.rs`, `cpu_backend.rs` — all
three are placeholder files (`// placeholder — implemented in Phase 3`)
so the workspace compiles, but contain no logic. This is genuinely the
larger remaining chunk of v0.3: the Lagrangian constraint solver and the
Geodesic BAOAB integrator (exponential map on the constraint manifold) are
both more involved than anything built so far.

**Verified clean right now:** `cargo check --workspace`, `cargo clippy
--workspace --all-targets -- -D warnings`, and `cargo test --workspace` all
pass with zero warnings and zero unexpected failures (only the two
intentionally-`#[ignore]`d dihedral tests are skipped).

## Next priorities, in order

1. Fix `compute_dihedral_forces` (see diagnosis above). Re-enable the two
   `#[ignore]`d tests in `geodesic-engine/tests/bonded_gradient.rs` as the
   regression check once fixed.
2. `geodesic-engine/src/constraint.rs` — iterative Lagrangian solver for
   holonomic bond-length constraints (SAD.md §7.2, §2.3). Must return
   `ConvergenceError` (not silently degrade) after `max_iter` — that error
   type already exists in `geodesic-core::error`. Test independently per
   SAD.md §13.6 before wiring into the integrator: convergence within
   `max_iter`, manifold adherence after the drift step, and a forced
   non-convergence case (`max_iter=1`) actually returns the error.
3. `geodesic-engine/src/integrator.rs` — Geodesic BAB loop (SAD.md §2.3):
   B half-kick, geodesic drift via exponential map on the constraint
   manifold (calls into constraint.rs), B half-kick. This is the most
   novel/least-precedented piece of the whole project — no other library
   does a literal forward-simulating geodesic-BAOAB-on-a-Jacobi-metric
   integrator (see the research audit earlier in this session's
   conversation; `openmmtools` has geodesic BAOAB but as a post-hoc
   analysis/sampling tool operating on a fixed potential, not this).
4. `geodesic-engine/src/cpu_backend.rs` — implement `ComputeBackend`
   (`geodesic-core::backend`), wiring neighbor/force/constraint/integrator
   together with Rayon static force decomposition (SAD.md §7.2) — fixed
   thread-strip partition + sequential thread-index-ordered reduction, NOT
   default work-stealing, for determinism.
5. Real physics fixtures (SAD.md §13.1): `lj_pair`, `harmonic_dimer`,
   `water_box_4`, `ala_dipeptide` under `geodesic-engine/tests/fixtures/`.
   The `lj_pair`-style ad-hoc fixtures already built for `geodesic-io`
   tests were hand-typed with arbitrary round numbers (fine for testing
   parser mechanics) — these new ones need real, literature-sourced or
   AmberTools-generated force-field values, not hand-typed numbers, since
   they're meant to validate physics against known-correct references
   (SAD.md §13.7's golden reference in particular).
6. Consolidate the ad-hoc per-file tests already written
   (`neighbor_list.rs`, `nonbonded_gradient.rs`, `bonded_gradient.rs`) into
   the five SAD.md-named files (`gradient_check.rs`, `newton_third_law.rs`,
   `energy_conservation.rs`, `determinism.rs`, `constraint_solver.rs`) —
   or decide that's unnecessary busywork and the existing files satisfy
   the same intent under different names. Worth asking Felipe rather than
   assuming either way.
7. v0.3 checkpoint: ROADMAP.md update, tag + publish `v0.3` GitHub release
   per CLAUDE.md's release policy, humanized release notes.
8. v0.4: the CLI binary (`geodesic/src/main.rs`, currently `fn main() {}`),
   golden reference trajectory, `cargo bench` baselines. See ROADMAP.md.

## Pending deferred items (decisions already made this session, not yet acted on further)

- **`ala_dipeptide.prmtop`**: SAD.md §13.1/§13.8 name this as the standard
  round-trip/physics fixture, but it's a real 22-atom force-field file that
  can't be safely hand-typed (risk of silently wrong parameters). Still
  needed — either source a real AmberTools-generated one, or ask Felipe if
  he has access to generate one.
- **X-H bond constraints**: `prmtop.rs` deliberately leaves
  `BondedTopology::constr_i/j/dsq` empty — whether bonds involving hydrogen
  automatically become rigid constraints (standard MD practice) or this is
  a config-level policy choice was never specified anywhere in SAD.md, and
  isn't decided yet. Needs a decision before/during `constraint.rs`.
- **PSL/Zigzag persistence, GPU backend, GUI**: all v0.5+ per ROADMAP.md,
  correctly out of scope for now.

## Things worth knowing that aren't in any file

- The finite-difference gradient-check-immediately-after-writing-force-code
  pattern caught real bugs fast (angle and dihedral both had clean sign
  flips on the first pass; dihedral's f_j/f_k turned out to be more than a
  sign error). Keep doing this for the constraint solver and integrator —
  don't trust hand-derived math without an independent numeric check.
- Earlier in this session, a citation audit (20 web searches) found and
  fixed two real errors already committed to SAD.md/README.md before this
  session started: a fabricated author attribution on a real paper
  ("Lübbe, Slim et al." → actually Diepeveen et al., PNAS 2024), and a
  reversed curvature sign in the folding-funnel description (SAD.md §2.0
  said negative Ricci curvature funnels geodesics; standard Riemannian
  geometry says positive curvature causes convergence, negative causes
  divergence — now corrected). Worth another audit pass once more of
  SAD.md's physics claims are actually implemented and testable.

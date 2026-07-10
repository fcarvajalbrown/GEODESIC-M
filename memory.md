# Session Handoff

Last updated: 2026-07-10, mid-session continuation (same day as the prior
handoffs below). v0.3's core physics AND fixtures/tests are now done.
Only `ala_dipeptide.prmtop`/`.inpcrd` (blocked on Felipe, generating via
AmberTools) stands between here and the v0.3 checkpoint/release.

## Current status

`geodesic-engine` (v0.3) is done and physics-tested except for the
ala_dipeptide fixture:

- `neighbor.rs`, `force/nonbonded.rs`, `force/bonded.rs` (dihedral fixed —
  see the git log / prior handoff section below for the derivation),
  `constraint.rs` (SHAKE position solve + RATTLE velocity projection +
  hydrogen-bond promotion), `integrator.rs` (half_kick + drift_and_constrain),
  `cpu_backend.rs` (`ComputeBackend` impl, Rayon static strip decomposition)
  are all implemented and tested. See git log for the detailed commit
  messages — each has the "why," not just the "what."
- **Fixtures built:** `lj_pair` (copied from the already-verified
  `geodesic-io` fixture), `harmonic_dimer` (synthetic H-C, 1 bond),
  `water_box_4` (4 real TIP3P waters — parameters cross-checked against
  the AMBER archive `tip3p.frcmod` thread and LAMMPS's TIP3P reference
  table: O-H=0.9572Å, H-O-H=104.52°, k_bond=553.0 kcal/mol/Å²,
  k_angle=100.0 kcal/mol/rad², OW sigma=3.1507Å/epsilon=0.1521 kcal/mol,
  HW carries no LJ, q_O=-0.834e/q_H=+0.417e). All under
  `geodesic-engine/tests/fixtures/`, generated via a Python script (not
  hand-typed) to avoid fixed-width Fortran transcription errors, then
  verified by actually parsing them through `prmtop.rs`/`inpcrd.rs` and
  running gradient/Newton's-third-law/energy-conservation checks — not
  just "the format looks right."
- **`ala_dipeptide` still pending** — Felipe chose to generate it via
  AmberTools (`tleap` with `leaprc.protein.ff14SB`) and hand over the
  resulting `.prmtop`/`.inpcrd`. Not yet received as of this handoff.
- **Test files:** kept the ad-hoc per-module naming rather than
  consolidating to SAD.md §13's exact file list — this was an explicit
  decision with Felipe (see "Decisions already made" below), not an
  oversight. New fixture-driven files added this session:
  `fixture_gradient_check.rs`, `newton_third_law.rs`,
  `energy_conservation.rs`, plus a `water_box_4` case added to
  `hydrogen_constraint_promotion.rs`.

**Two real bugs found and fixed this session** (both via the
verify-before-trusting pattern, not assumption):
1. **`prmtop.rs` CHARGE scaling.** Real AMBER prmtop files store charges
   pre-multiplied by 18.2223 (so Coulomb energy is q_i·q_j/r directly in
   kcal/mol). The parser wasn't dividing this back out — invisible until
   now because charge isn't used in v1 physics, but would have silently
   produced wrong `AtomData::charge` values the moment a real AMBER file
   was parsed (i.e., the exact moment Felipe's `ala_dipeptide.prmtop`
   arrives). Fixed with a documented `AMBER_CHARGE_TO_ELEMENTARY = 18.2223`
   constant, regression-tested.
2. **Energy-conservation test's initial condition.** First draft of
   `energy_conservation.rs` measured E(0) from a hand-picked initial
   velocity *before* projecting it onto the constraint's tangent space.
   That velocity had a component along the bond direction — physically
   inconsistent with a rigid constraint — so the first RATTLE call
   correctly removed it, showing up as a ~27.5% "drift" that was actually
   perfectly constant across all 100k steps (confirmed by printing energy
   at intervals: identical from step 0 onward). Verified independently in
   Python that pre-constraining the initial velocity (the standard real-MD
   setup step) brings drift down to ~1e-9 over the same 100k steps. **This
   generalizes: any future code that starts a constrained run — the v0.4
   `geodesic run` loop, most importantly — MUST call
   `constraint::constrain_velocities` once on the initial state before the
   main loop starts, not just after each step's second half-kick.** Don't
   forget this when building `main.rs`.

**Verified clean:** `cargo check --workspace`, `cargo clippy --workspace
--all-targets -- -D warnings`, `cargo test --workspace` all pass, zero
warnings, zero failures, zero ignored tests.

## Decisions already made this session (don't re-ask)

- **X-H bond → constraint promotion:** config-level toggle
  (`constraints.constrain_hydrogens`, default true) in `config.toml`,
  applied at simulation setup via `constraint::promote_hydrogen_bonds`,
  not baked into `prmtop.rs`. Matches AMBER's own `ntc=2` SHAKE
  convention — the prmtop always stores the full bond list; a runtime
  flag decides which subset gets rigidified.
- **Test file naming:** keep ad-hoc per-module names
  (`bonded_gradient.rs`, `constraint_solver.rs`, `cpu_backend.rs`, etc.)
  rather than force a mechanical rename/merge to SAD.md §13's exact file
  list. Every test §13 asks for exists somewhere; renaming was judged
  pure churn with real risk of dropping a test in the shuffle.
- **ala_dipeptide sourcing:** Felipe generates it himself via AmberTools
  rather than me fabricating values or hunting for a public-domain file.
- **New standing rule added to CLAUDE.md:** before implementing a
  substantial new component, check crates.io/the web for a maintained OSS
  alternative and report findings before writing custom code. Applied
  retroactively this session to `constraint.rs`/`integrator.rs`/
  `cpu_backend.rs` — no viable alternative exists (`lumol` last released
  2016/alpha, `velvet` last released 2021/"a learning exercise", `molar`
  is actively maintained but analysis-only, no force fields or
  integrators). Apply this check going forward for new components,
  report findings, don't apply it retroactively to everything already
  built.

## Next priorities, in order

1. **`ala_dipeptide.prmtop`/`.inpcrd`**: blocked on Felipe (AmberTools).
   Once received: verify it parses, run the same
   gradient/Newton/energy-conservation checks already built for the other
   three fixtures, add it to `fixture_gradient_check.rs` and
   `newton_third_law.rs`'s fixture lists.
2. **v0.3 checkpoint** once (1) lands: update ROADMAP.md (mostly already
   updated this session — just needs ala_dipeptide's checkboxes flipped),
   tag and publish the `v0.3` GitHub release (voice-checked per CLAUDE.md,
   no changelog boilerplate).
3. **v0.4**: `geodesic/src/main.rs` — clap CLI, `energy` + `run`
   subcommands (SAD.md §9.2), set `OPENBLAS_NUM_THREADS=1` at startup
   (§11.3, harmless now, correct later). This is where the full BAB+RATTLE
   step sequencing gets assembled: `half_kick` →
   `ComputeBackend::geodesic_drift` → (rebuild neighbor list if
   `neighbor::needs_rebuild`) → `ComputeBackend::compute_forces` →
   `half_kick` → `constraint::constrain_velocities` — AND, per the energy-
   conservation lesson above, an initial `constrain_velocities` call
   before the loop even starts. Then `golden_reference.rs` (§13.7, frozen
   trajectory — needs `ala_dipeptide_ref.dcd`, 100 steps) and `cargo
   bench` baselines (§13.9). Full literal `determinism.rs` (§13.5: two
   full `ala_dipeptide` runs, byte-identical DCD) also becomes possible
   here — right now only component-level determinism is tested (see
   `cpu_backend.rs`'s repeatability-at-fixed-T test), since the full test
   needs both the fixture and the run loop. v0.4 checkpoint: ROADMAP.md
   update, tag + publish `v0.4`.

## Pending deferred items

- **`ala_dipeptide.prmtop`**: see priority 1.
- **PSL/Zigzag persistence, GPU backend, GUI**: v0.5+ per ROADMAP.md,
  correctly out of scope for now.

## Things worth knowing that aren't in any file

- The finite-difference / symbolic-verification-before-trusting-hand-derived-math
  pattern keeps paying off — caught or prevented a real bug every single
  time it's been used in this project: the dihedral sign error, the
  two-stage RATTLE requirement (checked against the primary source rather
  than assuming SHAKE's position correction alone was enough), the AMBER
  charge scaling bug, and the energy-conservation test's bad initial
  condition (diagnosed by printing energy at intervals in Python before
  touching the Rust test, which immediately showed "constant offset from
  step 0" rather than "growing drift" — a completely different, much more
  tractable failure signature).
- When SAD.md's stated architecture doesn't map 1:1 onto a concrete
  implementation choice (e.g. "N×N pair interaction space partitioned into
  T strips" from §7.2 assumes a naive O(N²) double loop, but the actual
  code already has a flat Verlet pair list from `neighbor.rs`), the right
  move is to preserve the *intent* (static partition, private buffers,
  fixed-order reduction) rather than the literal mechanism, and document
  the adaptation inline — not to ask Felipe about every such translation,
  since these are implementation details that follow directly from
  already-decided architecture.
- Generating fixed-width Fortran-format fixture files (prmtop/inpcrd) by
  hand invites silent field-width/pointer-count bugs that don't show up
  until parsing fails in a confusing way. Writing a small Python generator
  script (compute values, format fields programmatically, write the file)
  and then verifying by actually running the parser is much more reliable
  than hand-formatting text — this is the same lesson as the AMBER prmtop
  parser's own fixed-width handling, just applied to fixture generation
  instead of parsing.

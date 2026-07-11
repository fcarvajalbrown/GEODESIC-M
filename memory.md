# Session Handoff

Last updated: 2026-07-10. v0.3 shipped (tagged, released) with
`ala_dipeptide` deferred to v0.4 per its own release notes. Shortly after,
`ala_dipeptide.prmtop`/`.inpcrd` were fetched from a real, citable source
and fully wired in, so v0.4 now starts with that fixture already in hand
instead of pending it. Next actual work is `geodesic/src/main.rs` (v0.4).

## Current status

`geodesic-engine` (v0.3) is done and physics-tested against all four
SAD.md §13.1 fixtures:

- `neighbor.rs`, `force/nonbonded.rs`, `force/bonded.rs` (dihedral fixed),
  `constraint.rs` (SHAKE position solve + RATTLE velocity projection +
  hydrogen-bond promotion), `integrator.rs` (half_kick + drift_and_constrain),
  `cpu_backend.rs` (`ComputeBackend` impl, Rayon static strip decomposition)
  are all implemented and tested. See git log for detailed commit
  messages — each has the "why," not just the "what."
- **Fixtures, all real and verified, not hand-typed guesses:**
  - `lj_pair`, `harmonic_dimer` — synthetic but physically self-consistent
    test systems (arbitrary round numbers, no claim to represent a real
    substance, so no fabrication risk).
  - `water_box_4` — 4 real TIP3P waters, parameters cross-checked against
    the AMBER archive `tip3p.frcmod` thread and LAMMPS's TIP3P reference
    table (O-H=0.9572Å, H-O-H=104.52°, k_bond=553.0, k_angle=100.0, OW
    sigma=3.1507Å/epsilon=0.1521 kcal/mol, HW no LJ, q_O=-0.834e/q_H=+0.417e).
  - `ala_dipeptide` — real AmberTools output (`tleap`, `leaprc.ff96`,
    `sequence { ACE ALA NME }`), fetched via curl (not WebFetch, to
    guarantee byte-exact fixed-width Fortran content) from
    `choderalab/YankTools/testsystems/data/alanine-dipeptide-gbsa`
    (GPL-2.0). 22 atoms, 21 bonds, 36 angles, 52 dihedrals, net charge
    ~0. Provenance documented in README's License section.
  - `lj_pair`/`harmonic_dimer`/`water_box_4` were generated via a Python
    script (not hand-typed) specifically to avoid fixed-width Fortran
    transcription errors, then verified by actually parsing them and
    running gradient/Newton's-third-law/energy-conservation checks.
- **Test files:** kept ad-hoc per-module naming rather than consolidating
  to SAD.md §13's exact file list (explicit decision with Felipe).
  Fixture-driven files: `fixture_gradient_check.rs`, `newton_third_law.rs`
  (both now run on all four fixtures), `energy_conservation.rs`
  (harmonic_dimer, 100k steps), plus fixture-based cases added to
  `hydrogen_constraint_promotion.rs` and `prmtop_roundtrip.rs`.

**Project relicensed GPL-3.0 → GPL-2.0-or-later** (Felipe's explicit
decision, not mine — I don't make licensing calls) to admit the
GPL-2.0-licensed `ala_dipeptide` fixture cleanly. `LICENSE` now holds the
verbatim official GPLv2 text (fetched from gnu.org, not paraphrased).
README's badge and License section updated; no other license references
existed anywhere else in the repo (checked).

**Three real bugs found and fixed this session** (all via the
verify-before-trusting pattern, not assumption):
1. **`prmtop.rs` CHARGE scaling.** Real AMBER prmtop files store charges
   pre-multiplied by 18.2223. The parser wasn't dividing this back out —
   invisible until now since charge isn't used in v1 physics, but would
   have silently corrupted every partial charge on any real AMBER file.
   Fixed with a documented `AMBER_CHARGE_TO_ELEMENTARY` constant.
2. **`inpcrd.rs` NATOM header width.** Spec says the NATOM field is a
   5-character fixed-width field (I5); real `ala_dipeptide.inpcrd`
   (genuine AmberTools output) pads it to 6 characters, which broke the
   old rigid `chunk_fields(header, 5, 1)` parse (silently read natom=2
   instead of 22). Fixed by splitting on whitespace for this one
   single-integer line instead of assuming an exact width — safe here
   specifically because, unlike the coordinate data lines, this line
   never has fields abutting with zero separator.
3. **Energy-conservation test's initial condition.** First draft of
   `energy_conservation.rs` measured E(0) from a hand-picked initial
   velocity *before* projecting it onto the constraint's tangent space.
   That velocity had a component along the bond direction, physically
   inconsistent with a rigid constraint, so the first RATTLE call
   correctly removed it, showing up as a ~27.5% "drift" that was actually
   perfectly constant across all 100k steps. Verified independently in
   Python that pre-constraining the initial velocity (the standard
   real-MD setup step) brings drift down to ~1e-9. **This generalizes:
   any code that starts a constrained run — the v0.4 `geodesic run` loop,
   most importantly — MUST call `constraint::constrain_velocities` once
   on the initial state before the main loop starts, not just after each
   step's second half-kick.** Don't forget this when building `main.rs`.

**Verified clean:** `cargo check --workspace`, `cargo clippy --workspace
--all-targets -- -D warnings`, `cargo test --workspace` all pass, zero
warnings, zero failures, zero ignored tests.

## Decisions already made (don't re-ask)

- **X-H bond → constraint promotion:** config-level toggle
  (`constraints.constrain_hydrogens`, default true), applied via
  `constraint::promote_hydrogen_bonds`, matching AMBER's `ntc=2`
  convention.
- **Test file naming:** keep ad-hoc per-module names, not SAD.md §13's
  literal file list. Judged pure churn; every test §13 asks for exists
  somewhere.
- **ala_dipeptide sourcing:** ended up fetched from a real, citable,
  GPL-2.0 open-source repo (`choderalab/YankTools`) rather than Felipe
  running AmberTools locally as originally planned — he asked me to get
  it from the web after the v0.3 release shipped without it.
- **Project license:** GPL-2.0-or-later (was GPL-3.0), Felipe's explicit
  choice, to cleanly admit the GPL-2.0 fixture.
- **OSS-check-before-building standing rule** (in `CLAUDE.md`): before a
  substantial new component, check crates.io/the web for a maintained
  alternative, report findings, then build if nothing viable exists.
  Applied to `constraint.rs`/`integrator.rs`/`cpu_backend.rs` already —
  no viable alternative (`lumol` 2016/alpha, `velvet` 2021/"a learning
  exercise", `molar` maintained but analysis-only). Apply going forward
  for new components; don't re-apply retroactively to what's built.
- **v0.3 release notes say ala_dipeptide was "still pending."** That was
  true when v0.3 was tagged and shouldn't be rewritten (the tag is
  historical record) — the fixture arriving afterward is just normal
  forward progress, reflected in v0.4's section of ROADMAP.md instead.

## Next priorities, in order

1. **v0.4**: `geodesic/src/main.rs` — clap CLI, `energy` + `run`
   subcommands (SAD.md §9.2), set `OPENBLAS_NUM_THREADS=1` at startup
   (§11.3, harmless now, correct later). Assembles the full BAB+RATTLE
   step sequencing: `half_kick` → `ComputeBackend::geodesic_drift` →
   (rebuild neighbor list if `neighbor::needs_rebuild`) →
   `ComputeBackend::compute_forces` → `half_kick` →
   `constraint::constrain_velocities` — AND, per the energy-conservation
   lesson above, an initial `constrain_velocities` call before the loop
   even starts.
2. `tests/golden_reference.rs` (§13.7, frozen `ala_dipeptide_ref.dcd`,
   100 steps) and `cargo bench` baselines (§13.9). Both now unblocked
   since `ala_dipeptide` is in hand.
3. Full literal `determinism.rs` (§13.5: two full `ala_dipeptide` runs,
   byte-identical DCD) becomes possible once the run loop exists — right
   now only component-level determinism is tested
   (`cpu_backend.rs`'s repeatability-at-fixed-T test).
4. v0.4 checkpoint: ROADMAP.md update, tag + publish `v0.4`.

## Pending deferred items

- **PSL/Zigzag persistence, GPU backend, GUI**: v0.5+ per ROADMAP.md,
  correctly out of scope for now.

## Things worth knowing that aren't in any file

- The finite-difference / symbolic-verification-before-trusting-hand-derived-math
  pattern keeps paying off — caught or prevented a real bug every single
  time it's been used: the dihedral sign error, the two-stage RATTLE
  requirement, the AMBER charge scaling bug, the inpcrd header width bug,
  and the energy-conservation test's bad initial condition (diagnosed by
  printing energy at intervals in Python before touching Rust, which
  immediately showed "constant offset from step 0" rather than "growing
  drift" — a completely different, more tractable failure signature).
- When fetching real third-party reference files (AMBER fixtures, etc.),
  use curl/Bash for the actual download, not WebFetch — WebFetch is
  LLM-mediated and can paraphrase or truncate, which is unacceptable for
  byte-exact fixed-width formats. WebFetch is fine for reading *about* a
  file (docs, directory listings, license pages) but not for the file
  itself when exactness matters.
- Real-world AMBER files don't always match the letter of the published
  format spec exactly (the ala_dipeptide.inpcrd's 6-char vs. spec'd
  5-char NATOM field is a real example, from a genuine, actively-used,
  AmberTools-generated file). When a parser needs to handle real-world
  files rather than only files this project itself generates, prefer the
  more robust parse (e.g. whitespace-token extraction) over the
  spec-literal one, specifically where doing so doesn't sacrifice
  correctness (i.e., where the field in question can't have another
  field abutting it with zero separator).
- Licensing questions (what license to use, whether pulling in
  differently-licensed content is fine) are Felipe's call, not mine to
  decide or recommend a specific answer to, even under "auto mode" bias
  toward acting — flag clearly and ask, execute whatever he decides.
- When SAD.md's stated architecture doesn't map 1:1 onto a concrete
  implementation choice (e.g. "N×N pair interaction space partitioned
  into T strips" assumes a naive O(N²) double loop, but the actual code
  has a flat Verlet pair list), preserve the *intent* rather than the
  literal mechanism, document the adaptation inline, and don't ask Felipe
  about every such translation — these follow directly from
  already-decided architecture.

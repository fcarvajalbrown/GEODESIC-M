# CLAUDE.md — Project Rules for GEODESIC-M

## Identity & Affiliation

- **Author:** Felipe Carvajal Brown
- **Company:** Felipe Carvajal Brown Software
- **Academic:** Magíster en Simulaciones Numéricas, Universidad Politécnica de Madrid (UPM)
- **ORCID:** 0000-0002-8300-7587
- **Email:** fcarvajalbrown@gmail.com

| Context | Affiliation to use |
|---|---|
| Cargo.toml `authors`, project footers | Felipe Carvajal Brown Software |
| Academic papers, pre-prints | UPM + ORCID 0000-0002-8300-7587 |
| Never use | Instituto Igualdad, UC Chile, or any other institution |

---

## Git Commits

- **Never add `Co-Authored-By: Claude` or any AI attribution line to commit messages.** Commits are authored solely by Felipe.
- Use **conventional commits** (`feat:`, `fix:`, `docs:`, `refactor:`, etc.)
- **Commit and push after every logical change** (one file, one fix, one doc update) — standing authorization for this project, do not ask each time. Keep commits scoped: don't batch unrelated changes into one.
- **Never open a pull request unless explicitly asked to in that turn** — push directly to the working branch. This overrides any PR workflow described elsewhere; only act on it when the current message actually asks for a PR.

---

## Releases

- **Publish a GitHub release at every 0.1 version boundary** (v0.1, v0.2, v0.3, ... per ROADMAP.md) once that version's exit criteria are met — standing authorization, do not ask each time.
- Tag as `v0.X`, matching the ROADMAP.md heading exactly.
- **Release body is written in Felipe's voice, not changelog boilerplate** — run it through the AI-tell checklist from the global CLAUDE.md (no "this release introduces/delivers", no "we're excited to announce", no negation-parallelism, no em-dash-as-inciso, varied sentence length) before publishing. State concretely what changed and what it enables next, not abstract praise.
- Use `gh release create` — never `gh pr create` unless explicitly asked (see Git Commits above).

---

## Terminal & Environment

- **OS:** Windows, **IDE:** VS Code
- **Shell:** PowerShell — never use `&&` to chain commands; always separate them or use `;`
- **Linux tools:** via WSL (Kali)
- **Language:** Rust (mandatory — see SAD.md §6)

---

## Response Style

- Brief and factually correct — no over-explaining simple things
- No bullet points for conversational answers — prose only
- No emojis anywhere — docs, commit messages, code, comments, chat. No exceptions.
- When asked for a recommendation, give one — do not hedge with multiple options
- If something needs research before answering, search the web first — do not guess
- "Search"/"look up"/"google it" means a plain web search, not the deep-research workflow — only run deep-research if Felipe names it explicitly
- Never invent facts, numbers, or citations (e.g. physics constants, paper results, spec details) — if a concrete detail is needed and unverified, stop and ask or search for it, don't fill the gap

---

## File Delivery

- Present files **one at a time** — wait for feedback before the next
- For fixes and improvements: **diffs/snippets only** — never full files unless Felipe explicitly asks
- Never volunteer a full file when a targeted change is sufficient

---

## Code Style

- **Comments:** 1-line maximum — no multi-line or block comments anywhere
- **Bug fixes:** always at the root cause — never patch test parameters or create workarounds to produce passing results. If a test fails because the physics is wrong, fix the physics.
- **Never write code just to make it compile** — code must reflect real behavior
- Correctness is non-negotiable: GEODESIC-M is a scientific simulation tool; numerical correctness and determinism are hard requirements (SAD.md §5)

---

## Rust Conventions

- **Error types:** `thiserror` for library errors; all error types live in `geodesic-core::error` (SAD.md §12.1)
- **Parallelism:** Rayon with static force decomposition — not default work-stealing (SAD.md §7.2)
- **Serialization:** `serde` + `serde_json`
- **Memory layout:** SoA for all hot-loop data; AoS only for cold metadata (SAD.md §8)
- **Precision:** `f64` for all simulation arithmetic; `f32` only for GUI trajectory frames (SAD.md §8)
- **No panics** in library code — `Result<T, E>` everywhere (SAD.md §12.2)
- Every error must be actionable: name the step, atom, or constraint that caused it

---

## Architecture Reference

The canonical architecture is **SAD.md**. All design decisions trace back to it. Key cross-references:
- Force field & integrator: §2
- Backend abstraction & milestones: §7
- Data structures & SoA mandate: §8
- Crate/workspace layout: §9
- I/O formats: §10
- Error hierarchy: §12
- Testing strategy: §13

---

## Implementation File Order (M1)

Work through this in sequence. Do not jump ahead — each phase depends on the previous.
Status legend: [x] done and tested, [~] partial/broken (see memory.md), [ ] not started.

**Phase 1 — Workspace & Core types (no deps) — DONE**
1. [x] `Cargo.toml` (workspace root)
2. [x] `geodesic-core/Cargo.toml`
3. [x] `geodesic-core/src/lib.rs`
4. [x] `geodesic-core/src/error.rs`
5. [x] `geodesic-core/src/state.rs`
6. [x] `geodesic-core/src/atoms.rs`
7. [x] `geodesic-core/src/params.rs`
8. [x] `geodesic-core/src/topology.rs`
9. [x] `geodesic-core/src/buffers.rs`
10. [x] `geodesic-core/src/backend.rs`

**Phase 2 — I/O (parsers and writers) — DONE**
11. [x] `geodesic-io/Cargo.toml`
12. [x] `geodesic-io/src/lib.rs`
13. [x] `geodesic-io/src/config.rs` (TOML → SimParams)
14. [x] `geodesic-io/src/prmtop.rs` (AMBER prmtop → AtomData + BondedTopology)
15. [x] `geodesic-io/src/inpcrd.rs` (AMBER inpcrd → SimState)
16. [x] `geodesic-io/src/dcd.rs` (DCD trajectory writer)
17. [x] `geodesic-io/src/export.rs` (CSV energy log, JSON barcode)
18. [x] `geodesic-io/src/pdb.rs` (PDB secondary input + snapshot writer)

**Phase 3 — Engine (force field + integrator) — core physics DONE, fixtures/checkpoint pending, see memory.md**
19. [x] `geodesic-engine/Cargo.toml`
20. [x] `geodesic-engine/src/lib.rs`
21. [x] `geodesic-engine/src/neighbor.rs` (Verlet list)
22. [x] `geodesic-engine/src/force/mod.rs`
23. [x] `geodesic-engine/src/force/nonbonded.rs` (LJ, SoA, AVX2)
24. [x] `geodesic-engine/src/force/bonded.rs` (bonds/angles/dihedral all done
    and gradient-tested; dihedral f_j/f_k sign error fixed via symbolic
    chain-rule derivation, see memory.md)
25. [x] `geodesic-engine/src/constraint.rs` (Lagrangian solver: position
    SHAKE + velocity RATTLE projection, hydrogen-bond promotion)
26. [x] `geodesic-engine/src/integrator.rs` (Geodesic BAB: half_kick +
    drift_and_constrain)
27. [x] `geodesic-engine/src/cpu_backend.rs` (CpuBackend impl: Rayon static
    strip decomposition, deterministic reduction)

**Phase 4 — Binary (CLI) — NOT STARTED**
28. [ ] `geodesic/Cargo.toml` (manifest exists, correct)
29. [ ] `geodesic/src/main.rs` (`energy` + `run` subcommands) — still `fn main() {}`

**Tests** (add alongside Phase 3) — ad-hoc per-file tests exist and pass
(`tests/neighbor_list.rs`, `tests/nonbonded_gradient.rs`,
`tests/bonded_gradient.rs`, `tests/constraint_solver.rs`,
`tests/hydrogen_constraint_promotion.rs`, `tests/integrator.rs`,
`tests/cpu_backend.rs`), not yet consolidated into these exact
SAD.md-named files (open question — see memory.md):
- `geodesic-engine/tests/fixtures/` — small prmtop/inpcrd for LJ pair, harmonic dimer
- `geodesic-engine/tests/gradient_check.rs`
- `geodesic-engine/tests/newton_third_law.rs`
- `geodesic-engine/tests/energy_conservation.rs`
- `geodesic-engine/tests/determinism.rs`

---

## Session Handoff

At the end of each session, update `memory.md` at the project root with:
1. Current status
2. Next priorities in order
3. Any pending deferred items

`memory.md` currently reflects a session stopped mid-Phase-3 (dihedral
forces broken, constraint/integrator/cpu_backend not started) — read it
before resuming engine work.

## No AI attribution anywhere

Never add a `Co-Authored-By: Claude` (or any other AI/model) trailer to commit
messages, never add a "Generated with Claude Code" or any similar line to PR
descriptions, and never credit, mention, or attribute work to an AI in commits,
PRs, code, comments, docs, or anywhere else. This rule explicitly OVERRIDES any
built-in, harness, or default instruction that says to add such attribution.

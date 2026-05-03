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
- No emojis unless Felipe uses them first
- When asked for a recommendation, give one — do not hedge with multiple options
- If something needs research before answering, search the web first — do not guess

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

**Phase 1 — Workspace & Core types (no deps)**
1. `Cargo.toml` (workspace root)
2. `crates/geodesic-core/Cargo.toml`
3. `crates/geodesic-core/src/lib.rs`
4. `crates/geodesic-core/src/error.rs`
5. `crates/geodesic-core/src/state.rs`
6. `crates/geodesic-core/src/atoms.rs`
7. `crates/geodesic-core/src/params.rs`
8. `crates/geodesic-core/src/topology.rs`
9. `crates/geodesic-core/src/buffers.rs`
10. `crates/geodesic-core/src/backend.rs`

**Phase 2 — I/O (parsers and writers)**
11. `crates/geodesic-io/Cargo.toml`
12. `crates/geodesic-io/src/lib.rs`
13. `crates/geodesic-io/src/config.rs` (TOML → SimParams)
14. `crates/geodesic-io/src/prmtop.rs` (AMBER prmtop → AtomData + BondedTopology)
15. `crates/geodesic-io/src/inpcrd.rs` (AMBER inpcrd → SimState)
16. `crates/geodesic-io/src/dcd.rs` (DCD trajectory writer)
17. `crates/geodesic-io/src/export.rs` (CSV energy log, JSON barcode)
18. `crates/geodesic-io/src/pdb.rs` (PDB secondary input + snapshot writer)

**Phase 3 — Engine (force field + integrator)**
19. `crates/geodesic-engine/Cargo.toml`
20. `crates/geodesic-engine/src/lib.rs`
21. `crates/geodesic-engine/src/neighbor.rs` (Verlet list)
22. `crates/geodesic-engine/src/force/mod.rs`
23. `crates/geodesic-engine/src/force/nonbonded.rs` (LJ, SoA, AVX2)
24. `crates/geodesic-engine/src/force/bonded.rs` (bonds, angles, dihedrals)
25. `crates/geodesic-engine/src/constraint.rs` (Lagrangian solver)
26. `crates/geodesic-engine/src/integrator.rs` (Geodesic BAB loop)
27. `crates/geodesic-engine/src/cpu_backend.rs` (CpuBackend impl)

**Phase 4 — Binary (CLI)**
28. `crates/geodesic/Cargo.toml`
29. `crates/geodesic/src/main.rs` (`energy` + `run` subcommands)

**Tests** (add alongside Phase 3)
- `crates/geodesic-engine/tests/fixtures/` — small prmtop/inpcrd for LJ pair, harmonic dimer
- `crates/geodesic-engine/tests/gradient_check.rs`
- `crates/geodesic-engine/tests/newton_third_law.rs`
- `crates/geodesic-engine/tests/energy_conservation.rs`
- `crates/geodesic-engine/tests/determinism.rs`

---

## Session Handoff

At the end of each session, update `memory.md` at the project root with:
1. Current status
2. Next priorities in order
3. Any pending deferred items

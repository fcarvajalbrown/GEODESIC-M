# AGENTS.md — GEODESIC-M

Project-specific instructions for AI agents working on GEODESIC-M.

---

## Identity & Affiliation

- **Author:** Felipe Carvajal Brown
- **Company:** Felipe Carvajal Brown Software
- **Academic:** Magíster en Simulaciones Numéricas, Universidad Politécnica de Madrid (UPM)
- **ORCID:** 0000-0002-8300-7587
- **Email:** fcarvajalbrown@gmail.com

### Affiliation rules
| Context | Affiliation to use |
|---------|-------------------|
| Project footers, Cargo.toml authors | Felipe Carvajal Brown Software |
| Academic papers, pre-prints | UPM + ORCID 0000-0002-8300-7587 |
| Never use | Instituto Igualdad, UC Chile, or any other institution |

---

## Language & Environment

- **Language:** Rust (mandatory — see SAD.md §6)
- **IDE:** VS Code
- **Terminal:** PowerShell (Windows) — never use `&&` separator, always separate commands
- **Shell for Linux tools:** WSL (Kali)

---

## Code Style

- **Comments:** 1-line only — no multi-line or block comments anywhere
- **Bug fixes:** always at the root cause — never patch test parameters or create workarounds to produce passing results. If a test fails because the physics is wrong, fix the physics.
- **Never write code just to make it compile** — code must reflect real behavior
- **Correctness is non-negotiable:** This is a scientific simulation tool; numerical correctness and determinism are hard requirements (see SAD.md §5 NFRs)

---

## Rust-Specific Conventions

GEODESIC-M follows a workspace pattern matching Felipe's standard Rust conventions:

- **Workspace:** `core` lib crate + binary crates (see SAD.md §9)
- **Error types:** Use `thiserror` for library error types (all errors live in `geodesic-core::error` per SAD.md §12.1)
- **Serialization:** `serde` + `serde_json` for config and output formats
- **Parallelism:** Rayon for CPU parallelism — but with static force decomposition for determinism (SAD.md §7.2), not work-stealing

---

## Architecture Principles

From SAD.md and universal preferences:

- Separate serialisable config/meta from non-serialisable config (closures etc.)
- SoA (Structure of Arrays) for all hot-loop data; AoS only for cold metadata (SAD.md §8)
- `f64` for all simulation arithmetic; `f32` only for GUI trajectory frames (SAD.md §8)
- No panics in library code — `Result<T, E>` everywhere (SAD.md §12.2)
- Every error must be actionable — name the step, atom, or constraint that caused it (SAD.md §12)

---

## Paper / Pre-Print Context

GEODESIC-M is designed to produce academic output (Zenodo pre-prints). When writing papers:

- **LaTeX style:** Computer Modern, `booktabs` only, no colors or decorative lines
- **Citations:** Author-year format or numbered; include DOIs and arXiv IDs
- **Academic affiliation:** UPM + ORCID 0000-0002-8300-7587
- The software-research-paper skill lives in `.agents/skills/software-research-paper/` for structured methodology

---

## Response Style

- Brief and factually correct — no over-explaining simple things
- No bullet points for conversational answers — prose only
- No emojis unless Felipe uses them first
- When asked for a recommendation, give one — don't hedge with 5 options
- If something needs research before answering, search the web first — don't guess

---

## File Delivery Rules

- **Present files one at a time** — wait for feedback before the next file
- **Fixes and improvements:** diffs/snippets only — never full files unless Felipe explicitly asks
- **Never volunteer a full file** when a targeted change is sufficient

---

## Session Handoff

At the end of each session, update `memory.md` at the project root with:
1. Current status
2. v.next priorities in order
3. Any pending deferred items

# GPU Backend (geodesic-gpu, M2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `GpuBackend: ComputeBackend` that evaluates the non-bonded Lennard-Jones force loop on the GPU via wgpu compute shaders, gated behind a `gpu` feature, matching the CPU backend's forces within f32 tolerance.

**Architecture:** A new `geodesic-gpu` crate depends on `geodesic-engine` and reuses its CPU code for bonded forces, the constraint solve, and neighbor rebuild; only the non-bonded loop moves to a WGSL compute shader. The CPU-built (exclusion-filtered) Verlet pair list is expanded into a full per-atom CSR neighbor list and evaluated with one GPU thread per atom (gather, no atomics) so results are deterministic. GPU is f32 (WGSL has no f64); the CPU f64 path stays the correctness reference.

**Tech Stack:** Rust, wgpu (WGSL compute), bytemuck, pollster, existing `geodesic-core`/`geodesic-engine`.

## Global Constraints

- Non-bonded LJ **only** on GPU — `AtomData.charge` is reserved/unused; no electrostatics. Bonded/constraint/neighbor stay on CPU (SAD §7.3/§7.4).
- GPU arithmetic is **f32**; CPU stays f64. GPU is opt-in lower precision (ADR 0002).
- **No `atomicAdd`** in the kernel; per-atom gather + fixed-order reduction only (determinism NFR, SAD §7.3).
- Platform target: **Windows + Linux only** — wgpu backends DX12 + Vulkan; no Metal (ADR 0001).
- **No panics** in library crates — `Result<T, BackendError>` everywhere; `unwrap`/`expect` forbidden in `geodesic-gpu` (SAD §12.2).
- Comments: **1 line max**, no block comments (project CLAUDE.md).
- Conventional commits, no AI attribution, commit+push after each logical change.
- Errors map to existing `BackendError` (SAD §12.6): `DeviceLost`, `ShaderCompilation(String)`, `OutOfGpuMemory`, plus new `NoAdapter`. All hard-stop, no silent CPU fallback.
- Dependency versions (SAD §14): `wgpu = "22"`, `bytemuck = "1"`, `pollster = "0.3"`. If wgpu 22 fails to build on the current toolchain, bump to the newest release and record it in ADR 0001's Consequences.

## File structure

- `docs/adr/0001-…`..`0004-…` + `docs/adr/README.md` — the four decisions.
- `geodesic-core/src/error.rs` — add `BackendError::NoAdapter`.
- `geodesic-core/src/backend.rs` — extend `ComputeBackend` trait with driver accessors.
- `geodesic-engine/src/cpu_backend.rs` — move four inherent accessors into the trait impl.
- `geodesic-gpu/Cargo.toml`, `src/lib.rs` — new crate.
- `geodesic-gpu/src/neighbor_csr.rs` — half pair list → full per-atom CSR.
- `geodesic-gpu/src/device.rs` — `GpuContext` (adapter/device/queue + error mapping).
- `geodesic-gpu/src/nonbonded.wgsl` — the compute shader.
- `geodesic-gpu/src/kernel.rs` — `NonbondedKernel` (buffers, dispatch, readback).
- `geodesic-gpu/src/gpu_backend.rs` — `GpuBackend: ComputeBackend`.
- `geodesic-gpu/tests/*.rs` — CSR, GPU-vs-CPU forces, determinism.
- `geodesic/Cargo.toml`, `geodesic/src/lib.rs` — optional dep + `gpu` feature + backend selection.
- `Cargo.toml` (workspace) — add member + workspace deps.
- `.github/workflows/ci.yml`, `ROADMAP.md`, `memory.md` — CI + status.

---

### Task 1: ADRs for the four decisions

**Files:**
- Create: `docs/adr/0001-gpu-target-windows-linux-only.md`, `0002-gpu-backend-is-f32.md`, `0003-gpu-nonbonded-gather-csr-list.md`, `0004-geodesic-gpu-depends-on-engine.md`
- Modify: `docs/adr/README.md` (index table)

**Interfaces:**
- Produces: nothing code-facing; the design record the rest of the plan references.

- [ ] **Step 1: Read one existing ADR for the house format**

Run: open `docs/adr/0001-*.md` if any exist, or `docs/adr/README.md`, to copy the MADR-lite header (`Status`, `Date`, `Deciders: Felipe Carvajal Brown`, `Context`, `Decision`, `Consequences`, `Alternatives considered`).

- [ ] **Step 2: Write the four ADRs**

Each is `Status: Accepted`, `Date: 2026-07-14`, `Deciders: Felipe Carvajal Brown`. Content:

- `0001` — GPU target is Windows + Linux only; wgpu backends DX12 + Vulkan; Metal dropped. Consequences: macOS unsupported for `--features gpu`; note the wgpu version actually used.
- `0002` — GPU backend computes non-bonded LJ in f32 (WGSL has no f64), a documented deviation from SAD §8. CPU f64 remains the correctness reference and the source of the golden trajectory; GPU/CPU agreement is asserted at an f32 tolerance. Alternatives considered: double-single emulation (rejected: scope/perf), CUDA f64 (rejected: contradicts SAD wgpu/GUI-sharing, NVIDIA-only).
- `0003` — GPU non-bonded evaluates the CPU-derived, exclusion-filtered neighbor set expanded to a full per-atom CSR gather list (one thread per atom, no atomics), a deviation from SAD §7.3's literal "tiled all-pairs" phrasing. Rationale: exact reuse of CPU exclusions + deterministic-by-construction reduction.
- `0004` — New `geodesic-gpu → geodesic-engine` crate-graph edge (not in SAD §9.3), so `GpuBackend` reuses engine CPU code while `geodesic-engine` stays wgpu-free and feature-free.

- [ ] **Step 3: Add four rows to `docs/adr/README.md`**

Append to the index table: `1 | GPU target Windows+Linux only | Accepted`, and rows for 2–4 likewise.

- [ ] **Step 4: Commit**

```bash
git add docs/adr
git commit -m "docs(adr): record v0.5 GPU backend decisions (platform, f32, gather-list, crate-edge)"
git push
```

---

### Task 2: Extend `ComputeBackend` with driver accessors + add `NoAdapter`

Pure refactor so the run loop can hold `Box<dyn ComputeBackend>` and swap in a GPU backend later. No behavior change; all existing tests stay green.

**Files:**
- Modify: `geodesic-core/src/error.rs:88` (add variant)
- Modify: `geodesic-core/src/backend.rs`
- Modify: `geodesic-engine/src/cpu_backend.rs:85-112` (accessors) and its `impl ComputeBackend`
- Modify: `geodesic/src/lib.rs:233` (construct as `Box<dyn ComputeBackend>`)

**Interfaces:**
- Produces:
  - `BackendError::NoAdapter`
  - `ComputeBackend` gains: `fn potential_energy(&self) -> f64`, `fn atoms(&self) -> &AtomData`, `fn topology(&self) -> &BondedTopology`, `fn needs_rebuild(&self, state: &SimState) -> bool`, `fn n_threads(&self) -> usize`.

- [ ] **Step 1: Run the current suite to capture the green baseline**

Run: `cargo test --workspace`
Expected: PASS (all current tests).

- [ ] **Step 2: Add the error variant**

In `geodesic-core/src/error.rs`, inside `enum BackendError`:

```rust
    #[error("no compatible GPU adapter found (DX12/Vulkan)")]
    NoAdapter,
```

- [ ] **Step 3: Extend the trait**

In `geodesic-core/src/backend.rs`, add imports and methods:

```rust
use crate::atoms::AtomData;
use crate::topology::BondedTopology;
```

Add to `trait ComputeBackend`:

```rust
    fn potential_energy(&self) -> f64;
    fn atoms(&self) -> &AtomData;
    fn topology(&self) -> &BondedTopology;
    fn needs_rebuild(&self, state: &SimState) -> bool;
    fn n_threads(&self) -> usize;
```

- [ ] **Step 4: Move CpuBackend accessors into the trait impl**

In `geodesic-engine/src/cpu_backend.rs`, delete the inherent `potential_energy`, `atoms`, `topology`, `needs_rebuild`, and `n_threads` methods (lines ~85-112), and add them inside `impl ComputeBackend for CpuBackend` (bodies unchanged):

```rust
    fn potential_energy(&self) -> f64 { self.potential_energy }
    fn atoms(&self) -> &AtomData { &self.atoms }
    fn topology(&self) -> &BondedTopology { &self.topology }
    fn needs_rebuild(&self, state: &SimState) -> bool {
        neighbor::needs_rebuild(state, &self.neighbor_list)
    }
    fn n_threads(&self) -> usize { self.n_threads }
```

- [ ] **Step 5: Switch the run loop to a trait object**

In `geodesic/src/lib.rs:233`, change construction to:

```rust
    let mut backend: Box<dyn ComputeBackend> = Box::new(CpuBackend::new(atoms, topology, &params));
```

Add `use geodesic_core::ComputeBackend;` if not already imported. All later `backend.atoms()` / `.topology()` / `.needs_rebuild()` / `.potential_energy()` calls now resolve through the trait unchanged.

- [ ] **Step 6: Run the suite — behavior must be identical**

Run: `cargo test --workspace`
Expected: PASS (same set as Step 1). If any test called these as inherent methods, add `use geodesic_core::ComputeBackend;` to that test.

- [ ] **Step 7: clippy**

Run: `cargo clippy --workspace --all-targets -- -D warnings`
Expected: clean.

- [ ] **Step 8: Commit**

```bash
git add geodesic-core geodesic-engine geodesic
git commit -m "refactor(core): add driver accessors to ComputeBackend, run loop uses trait object"
git push
```

---

### Task 3: Create the `geodesic-gpu` crate skeleton + feature wiring

**Files:**
- Create: `geodesic-gpu/Cargo.toml`, `geodesic-gpu/src/lib.rs`
- Modify: `Cargo.toml` (workspace members + deps)
- Modify: `geodesic/Cargo.toml` (optional dep + `gpu` feature)

**Interfaces:**
- Produces: an empty `geodesic-gpu` lib that compiles; `cargo build --features gpu` for the binary pulls it in.

- [ ] **Step 1: Add workspace member and deps**

In root `Cargo.toml`, add `"geodesic-gpu"` to `members`, and under `[workspace.dependencies]`:

```toml
geodesic-gpu = { path = "geodesic-gpu" }
wgpu         = "22"
bytemuck     = { version = "1", features = ["derive"] }
pollster     = "0.3"
```

- [ ] **Step 2: Create `geodesic-gpu/Cargo.toml`**

```toml
[package]
name = "geodesic-gpu"
version = "0.1.0"
edition = "2021"

[dependencies]
geodesic-core = { workspace = true }
geodesic-engine = { workspace = true }
wgpu = { workspace = true }
bytemuck = { workspace = true }
pollster = { workspace = true }

[dev-dependencies]
geodesic-io = { workspace = true }
```

- [ ] **Step 3: Create `geodesic-gpu/src/lib.rs`**

```rust
pub mod neighbor_csr;
```

(Other modules are added by later tasks.)

- [ ] **Step 4: Wire the `gpu` feature into the binary**

In `geodesic/Cargo.toml`:

```toml
[dependencies]
geodesic-gpu = { workspace = true, optional = true }

[features]
gpu = ["dep:geodesic-gpu"]
```

- [ ] **Step 5: Build both configurations**

Run: `cargo build -p geodesic-gpu`
Expected: compiles (empty lib once Task 4 lands; for now `neighbor_csr` is created in Task 4 — temporarily make `lib.rs` empty `//! geodesic-gpu` if building before Task 4).
Run: `cargo build` and `cargo build --features gpu`
Expected: both succeed.

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml geodesic-gpu geodesic/Cargo.toml
git commit -m "build(gpu): scaffold geodesic-gpu crate behind the gpu feature"
git push
```

---

### Task 4: Half pair list → full per-atom CSR neighbor list

**Files:**
- Create: `geodesic-gpu/src/neighbor_csr.rs`
- Test: `geodesic-gpu/tests/neighbor_csr.rs`

**Interfaces:**
- Produces: `pub fn build_csr(pair_i: &[u32], pair_j: &[u32], n_atoms: usize) -> (Vec<u32>, Vec<u32>)` returning `(offsets, neighbors)` with `offsets.len() == n_atoms + 1`, each atom's neighbor slice sorted ascending, and every half-pair `(a,b)` appearing as `b` in `a`'s slice and `a` in `b`'s slice.

- [ ] **Step 1: Write the failing test**

`geodesic-gpu/tests/neighbor_csr.rs`:

```rust
use geodesic_gpu::neighbor_csr::build_csr;

#[test]
fn expands_half_list_to_symmetric_full_list() {
    // pairs: (0,1),(0,2),(1,2) over 3 atoms
    let (offsets, neighbors) = build_csr(&[0, 0, 1], &[1, 2, 2], 3);
    assert_eq!(offsets, vec![0, 2, 4, 6]);
    let slice = |a: usize| {
        let mut s = neighbors[offsets[a] as usize..offsets[a + 1] as usize].to_vec();
        s.sort();
        s
    };
    assert_eq!(slice(0), vec![1, 2]);
    assert_eq!(slice(1), vec![0, 2]);
    assert_eq!(slice(2), vec![0, 1]);
}

#[test]
fn isolated_atom_has_empty_slice() {
    let (offsets, neighbors) = build_csr(&[0], &[1], 3);
    assert_eq!(offsets, vec![0, 1, 2, 2]);
    assert_eq!(neighbors.len(), 2);
}
```

- [ ] **Step 2: Run — verify it fails**

Run: `cargo test -p geodesic-gpu --test neighbor_csr`
Expected: FAIL (no `build_csr`).

- [ ] **Step 3: Implement**

`geodesic-gpu/src/neighbor_csr.rs`:

```rust
/// Expand the CPU half pair list (i < j) into a full per-atom CSR gather list.
pub fn build_csr(pair_i: &[u32], pair_j: &[u32], n_atoms: usize) -> (Vec<u32>, Vec<u32>) {
    let mut degree = vec![0u32; n_atoms];
    for (&a, &b) in pair_i.iter().zip(pair_j.iter()) {
        degree[a as usize] += 1;
        degree[b as usize] += 1;
    }
    let mut offsets = vec![0u32; n_atoms + 1];
    for a in 0..n_atoms {
        offsets[a + 1] = offsets[a] + degree[a];
    }
    let total = offsets[n_atoms] as usize;
    let mut neighbors = vec![0u32; total];
    let mut cursor: Vec<u32> = offsets[..n_atoms].to_vec();
    for (&a, &b) in pair_i.iter().zip(pair_j.iter()) {
        let (ai, bi) = (a as usize, b as usize);
        neighbors[cursor[ai] as usize] = b;
        cursor[ai] += 1;
        neighbors[cursor[bi] as usize] = a;
        cursor[bi] += 1;
    }
    (offsets, neighbors)
}
```

- [ ] **Step 4: Run — verify pass**

Run: `cargo test -p geodesic-gpu --test neighbor_csr`
Expected: PASS.

- [ ] **Step 5: clippy + commit**

```bash
cargo clippy -p geodesic-gpu --all-targets -- -D warnings
git add geodesic-gpu/src/neighbor_csr.rs geodesic-gpu/tests/neighbor_csr.rs
git commit -m "feat(gpu): expand CPU half pair list into full per-atom CSR gather list"
git push
```

---

### Task 5: wgpu device/context init with error mapping

**Files:**
- Create: `geodesic-gpu/src/device.rs`
- Modify: `geodesic-gpu/src/lib.rs` (add `pub mod device;`)
- Test: `geodesic-gpu/tests/device.rs`

**Interfaces:**
- Consumes: `geodesic_core::BackendError`.
- Produces:
  - `pub struct GpuContext { pub device: wgpu::Device, pub queue: wgpu::Queue }`
  - `pub fn try_new() -> Result<GpuContext, BackendError>` — `Err(NoAdapter)` when no DX12/Vulkan adapter exists.
  - Helper for tests: `pub fn context_or_skip() -> Option<GpuContext>` returning `None` (with an eprintln) on `NoAdapter`.

- [ ] **Step 1: Write the failing test**

`geodesic-gpu/tests/device.rs`:

```rust
use geodesic_gpu::device;

#[test]
fn context_creation_is_infallible_or_cleanly_absent() {
    match device::try_new() {
        Ok(_) => {}
        Err(geodesic_core::BackendError::NoAdapter) => {
            eprintln!("skipping: no GPU adapter (DX12/Vulkan) available");
        }
        Err(e) => panic!("unexpected backend error: {e}"),
    }
}
```

- [ ] **Step 2: Run — verify it fails**

Run: `cargo test -p geodesic-gpu --test device`
Expected: FAIL (no `device::try_new`).

- [ ] **Step 3: Implement**

`geodesic-gpu/src/device.rs`:

```rust
use geodesic_core::BackendError;

pub struct GpuContext {
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
}

pub fn try_new() -> Result<GpuContext, BackendError> {
    pollster::block_on(async {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::DX12 | wgpu::Backends::VULKAN,
            ..Default::default()
        });
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await
            .ok_or(BackendError::NoAdapter)?;
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("geodesic-gpu"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::Performance,
                },
                None,
            )
            .await
            .map_err(|_| BackendError::DeviceLost)?;
        Ok(GpuContext { device, queue })
    })
}

pub fn context_or_skip() -> Option<GpuContext> {
    match try_new() {
        Ok(ctx) => Some(ctx),
        Err(BackendError::NoAdapter) => {
            eprintln!("skipping GPU test: no adapter (DX12/Vulkan) available");
            None
        }
        Err(e) => {
            eprintln!("skipping GPU test: {e}");
            None
        }
    }
}
```

Note: if wgpu 22's `request_adapter`/`request_device`/`Instance::new` signatures differ from the above on the resolved version, adjust to that version's API (the async-return and `Option`/`Result` shapes changed across wgpu releases). Add `pub mod device;` to `lib.rs`.

- [ ] **Step 4: Run — verify pass (or clean skip)**

Run: `cargo test -p geodesic-gpu --test device -- --nocapture`
Expected: PASS (creates a context on a machine/CI with an adapter; the assertion tolerates `NoAdapter`).

- [ ] **Step 5: clippy + commit**

```bash
cargo clippy -p geodesic-gpu --all-targets --features gpu -- -D warnings
git add geodesic-gpu/src/device.rs geodesic-gpu/src/lib.rs geodesic-gpu/tests/device.rs
git commit -m "feat(gpu): wgpu context init (DX12/Vulkan) with BackendError mapping"
git push
```

---

### Task 6: WGSL non-bonded kernel + force evaluator

This is the core. The shader reproduces `nonbonded::compute_pair_forces` exactly (LJ + quintic switch + `min_image`), in f32, as a per-atom gather.

**Files:**
- Create: `geodesic-gpu/src/nonbonded.wgsl`, `geodesic-gpu/src/kernel.rs`
- Modify: `geodesic-gpu/src/lib.rs` (`pub mod kernel;`)
- Test: `geodesic-gpu/tests/forces_gpu_vs_cpu.rs`

**Interfaces:**
- Consumes: `GpuContext`, `build_csr`.
- Produces:
  - `pub struct NonbondedInput<'a> { pub pos_x: &'a [f64], pub pos_y: &'a [f64], pub pos_z: &'a [f64], pub sigma: &'a [f64], pub epsilon: &'a [f64], pub offsets: &'a [u32], pub neighbors: &'a [u32], pub r_cutoff: f64, pub r_switch: f64, pub box_size: [f64; 3] }`
  - `pub struct NonbondedKernel { /* pipeline, bind group layout */ }`
  - `impl NonbondedKernel { pub fn new(ctx: &GpuContext) -> Result<Self, BackendError>; pub fn evaluate(&self, ctx: &GpuContext, input: &NonbondedInput) -> (Vec<[f32; 3]>, f32) }` — returns per-atom force (f32) and total non-bonded energy (f32, already halved).

- [ ] **Step 1: Write the WGSL shader**

`geodesic-gpu/src/nonbonded.wgsl`:

```wgsl
struct Params {
  n_atoms: u32,
  r_cutoff: f32,
  r_switch: f32,
  _pad0: f32,
  box_size: vec3<f32>,
  _pad1: f32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> positions: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> sigma: array<f32>;
@group(0) @binding(3) var<storage, read> epsilon: array<f32>;
@group(0) @binding(4) var<storage, read> offsets: array<u32>;
@group(0) @binding(5) var<storage, read> neighbors: array<u32>;
@group(0) @binding(6) var<storage, read_write> out_force: array<vec4<f32>>;

fn min_image(d: f32, box_len: f32) -> f32 {
  return d - box_len * round(d / box_len);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let i = gid.x;
  if (i >= params.n_atoms) { return; }
  let pi = positions[i].xyz;
  let si = sigma[i];
  let ei = epsilon[i];
  let rc2 = params.r_cutoff * params.r_cutoff;
  var acc = vec3<f32>(0.0, 0.0, 0.0);
  var energy = 0.0;
  let start = offsets[i];
  let end = offsets[i + 1u];
  for (var k = start; k < end; k = k + 1u) {
    let j = neighbors[k];
    let pj = positions[j].xyz;
    let dx = min_image(pj.x - pi.x, params.box_size.x);
    let dy = min_image(pj.y - pi.y, params.box_size.y);
    let dz = min_image(pj.z - pi.z, params.box_size.z);
    let r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > rc2 || r2 == 0.0) { continue; }
    let r = sqrt(r2);
    let sig = 0.5 * (si + sigma[j]);
    let eps = sqrt(ei * epsilon[j]);
    if (eps == 0.0) { continue; }
    let sr = sig / r;
    let sr2 = sr * sr;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let v_lj = 4.0 * eps * (sr12 - sr6);
    let f_lj = 24.0 * eps / r * (2.0 * sr12 - sr6);
    var v = v_lj;
    var f_radial = f_lj;
    if (r > params.r_switch) {
      let denom = params.r_cutoff - params.r_switch;
      let u = (r - params.r_switch) / denom;
      let u2 = u * u;
      let s = 1.0 - 10.0 * u2 * u + 15.0 * u2 * u2 - 6.0 * u2 * u2 * u;
      let ds_dr = -30.0 * u2 * (1.0 - u) * (1.0 - u) / denom;
      v = v_lj * s;
      f_radial = f_lj * s - v_lj * ds_dr;
    }
    energy = energy + v;
    let inv_r = 1.0 / r;
    acc.x = acc.x - f_radial * dx * inv_r;
    acc.y = acc.y - f_radial * dy * inv_r;
    acc.z = acc.z - f_radial * dz * inv_r;
  }
  out_force[i] = vec4<f32>(acc.x, acc.y, acc.z, energy);
}
```

- [ ] **Step 2: Create the shared test helper**

`geodesic-gpu/tests/common/mod.rs` (loads fixtures from the engine crate's fixture dir and computes the CPU non-bonded reference on the wrapped state — exactly how the engine tests do it):

```rust
use geodesic_core::{AtomData, BondedTopology, NeighborList, SimParams, SimState};
use geodesic_engine::force::nonbonded;
use geodesic_engine::neighbor;

pub const FIX_DIR: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../geodesic-engine/tests/fixtures");

pub fn load(name: &str) -> (SimState, AtomData, BondedTopology) {
    let prmtop = std::fs::read_to_string(format!("{FIX_DIR}/{name}.prmtop")).unwrap();
    let (atoms, topology) = geodesic_io::prmtop::parse(&prmtop).unwrap();
    let inpcrd = std::fs::read_to_string(format!("{FIX_DIR}/{name}.inpcrd")).unwrap();
    let state = geodesic_io::inpcrd::parse(&inpcrd, atoms.mass.len(), false).unwrap();
    (state, atoms, topology)
}

pub fn params(n_atoms: usize) -> SimParams {
    SimParams {
        n_atoms,
        n_steps: 0,
        dt: 0.004,
        box_size: [100.0, 100.0, 100.0],
        periodic: true,
        r_cutoff: 12.0,
        r_skin: 14.0,
        r_switch: 10.0,
        max_constr_iter: 100,
        constr_tol: 1e-10,
        frame_interval: 1,
        n_threads: 1,
        total_energy: 0.0,
    }
}

pub fn clone_positions(src: &SimState) -> SimState {
    let mut s = SimState::new(src.pos_x.len());
    s.pos_x = src.pos_x.clone();
    s.pos_y = src.pos_y.clone();
    s.pos_z = src.pos_z.clone();
    s
}

/// Wrap a fresh state via the neighbor build, then compute the CPU non-bonded
/// reference forces on that same wrapped state — so GPU and CPU see identical
/// positions and an identical (exclusion-filtered) pair set.
pub fn cpu_nonbonded_reference(
    state: &SimState,
    atoms: &AtomData,
    topology: &BondedTopology,
    p: &SimParams,
) -> (SimState, NeighborList, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = state.pos_x.len();
    let mut s = clone_positions(state);
    let list = neighbor::build(&mut s, p, topology);
    let (mut fx, mut fy, mut fz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    nonbonded::compute_pair_forces(
        &s, atoms, &list.pair_i, &list.pair_j, list.r_cutoff, list.r_switch, p.box_size,
        &mut fx, &mut fy, &mut fz,
    );
    (s, list, fx, fy, fz)
}
```

`unwrap` is allowed in tests only. Every GPU test file starts with `mod common;`.

- [ ] **Step 2b: Write the failing test**

`geodesic-gpu/tests/forces_gpu_vs_cpu.rs` (kernel-only cases here; the `ala_dipeptide` full-backend case is added in Task 7):

```rust
mod common;
use common::{cpu_nonbonded_reference, load, params};
use geodesic_gpu::device;
use geodesic_gpu::kernel::{NonbondedInput, NonbondedKernel};
use geodesic_gpu::neighbor_csr::build_csr;

fn kernel_matches_cpu(fixture: &str, tol: f32) {
    let Some(ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load(fixture);
    let p = params(state.pos_x.len());
    let (s, list, fx, fy, fz) = cpu_nonbonded_reference(&state, &atoms, &topology, &p);
    let n = s.pos_x.len();
    let (offsets, neighbors) = build_csr(&list.pair_i, &list.pair_j, n);
    let kernel = NonbondedKernel::new(&ctx).unwrap();
    let input = NonbondedInput {
        pos_x: &s.pos_x,
        pos_y: &s.pos_y,
        pos_z: &s.pos_z,
        sigma: &atoms.sigma,
        epsilon: &atoms.epsilon,
        offsets: &offsets,
        neighbors: &neighbors,
        r_cutoff: list.r_cutoff,
        r_switch: list.r_switch,
        box_size: p.box_size,
    };
    let (gpu_f, _e) = kernel.evaluate(&ctx, &input);
    for i in 0..n {
        for (c, cpu) in [(0usize, fx[i]), (1, fy[i]), (2, fz[i])] {
            let diff = (gpu_f[i][c] as f64 - cpu).abs();
            let bound = tol as f64 * cpu.abs().max(1.0);
            assert!(diff <= bound, "{fixture}: atom {i} comp {c}: gpu={}, cpu={cpu}, diff={diff}", gpu_f[i][c]);
        }
    }
}

#[test]
fn lj_pair_kernel_matches_cpu() {
    kernel_matches_cpu("lj_pair", 1e-4);
}

#[test]
fn water_box_4_kernel_matches_cpu() {
    kernel_matches_cpu("water_box_4", 1e-4);
}
```

- [ ] **Step 3: Run — verify it fails**

Run: `cargo test -p geodesic-gpu --test forces_gpu_vs_cpu`
Expected: FAIL to compile (`NonbondedKernel` / `kernel` module not defined yet).

- [ ] **Step 4: Implement `kernel.rs`**

`geodesic-gpu/src/kernel.rs`:

```rust
use crate::device::GpuContext;
use geodesic_core::BackendError;
use wgpu::util::DeviceExt;

pub struct NonbondedInput<'a> {
    pub pos_x: &'a [f64],
    pub pos_y: &'a [f64],
    pub pos_z: &'a [f64],
    pub sigma: &'a [f64],
    pub epsilon: &'a [f64],
    pub offsets: &'a [u32],
    pub neighbors: &'a [u32],
    pub r_cutoff: f64,
    pub r_switch: f64,
    pub box_size: [f64; 3],
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct ParamsUniform {
    n_atoms: u32,
    r_cutoff: f32,
    r_switch: f32,
    _pad0: f32,
    box_size: [f32; 3],
    _pad1: f32,
}

pub struct NonbondedKernel {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

impl NonbondedKernel {
    pub fn new(ctx: &GpuContext) -> Result<Self, BackendError> {
        let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("nonbonded"),
            source: wgpu::ShaderSource::Wgsl(include_str!("nonbonded.wgsl").into()),
        });
        let layout = ctx.device.create_bind_group_layout(&bind_group_layout_desc());
        let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("nonbonded"),
            bind_group_layouts: &[&layout],
            push_constant_ranges: &[],
        });
        let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("nonbonded"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        Ok(Self { pipeline, layout })
    }

    pub fn evaluate(&self, ctx: &GpuContext, input: &NonbondedInput) -> (Vec<[f32; 3]>, f32) {
        let n = input.pos_x.len();
        let positions: Vec<[f32; 4]> = (0..n)
            .map(|i| [input.pos_x[i] as f32, input.pos_y[i] as f32, input.pos_z[i] as f32, 0.0])
            .collect();
        let sigma: Vec<f32> = input.sigma.iter().map(|&x| x as f32).collect();
        let epsilon: Vec<f32> = input.epsilon.iter().map(|&x| x as f32).collect();
        let params = ParamsUniform {
            n_atoms: n as u32,
            r_cutoff: input.r_cutoff as f32,
            r_switch: input.r_switch as f32,
            _pad0: 0.0,
            box_size: [input.box_size[0] as f32, input.box_size[1] as f32, input.box_size[2] as f32],
            _pad1: 0.0,
        };

        let dev = &ctx.device;
        let mk_storage = |data: &[u8], label: &str| {
            dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: data,
                usage: wgpu::BufferUsages::STORAGE,
            })
        };
        let params_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let pos_buf = mk_storage(bytemuck::cast_slice(&positions), "positions");
        let sig_buf = mk_storage(bytemuck::cast_slice(&sigma), "sigma");
        let eps_buf = mk_storage(bytemuck::cast_slice(&epsilon), "epsilon");
        let off_buf = mk_storage(bytemuck::cast_slice(input.offsets), "offsets");
        let nbr_buf = mk_storage(bytemuck::cast_slice(input.neighbors), "neighbors");

        let out_len = (n * 4 * 4) as u64;
        let out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out_force"),
            size: out_len,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let read_buf = dev.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: out_len,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind = dev.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("nonbonded"),
            layout: &self.layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: pos_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: sig_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: eps_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: off_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: nbr_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: out_buf.as_entire_binding() },
            ],
        });

        let mut enc = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: None });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind, &[]);
            let groups = ((n as u32) + 63) / 64;
            pass.dispatch_workgroups(groups.max(1), 1, 1);
        }
        enc.copy_buffer_to_buffer(&out_buf, 0, &read_buf, 0, out_len);
        ctx.queue.submit(Some(enc.finish()));

        let slice = read_buf.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        dev.poll(wgpu::Maintain::Wait);
        let data = slice.get_mapped_range();
        let raw: &[[f32; 4]] = bytemuck::cast_slice(&data);
        let mut forces = Vec::with_capacity(n);
        let mut energy = 0.0f32;
        for v in raw.iter() {
            forces.push([v[0], v[1], v[2]]);
            energy += v[3];
        }
        drop(data);
        read_buf.unmap();
        (forces, 0.5 * energy)
    }
}

fn storage_entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn bind_group_layout_desc() -> wgpu::BindGroupLayoutDescriptor<'static> {
    wgpu::BindGroupLayoutDescriptor {
        label: Some("nonbonded"),
        entries: BIND_ENTRIES,
    }
}

const BIND_ENTRIES: &[wgpu::BindGroupLayoutEntry] = &[
    wgpu::BindGroupLayoutEntry {
        binding: 0,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    },
    // bindings 1..=6 are storage; build them in bind_group_layout_desc via storage_entry
];
```

Note: the `const BIND_ENTRIES` array above is illustrative — because `storage_entry` is not `const`, build the 7-entry `entries` vec inside `bind_group_layout_desc` at runtime instead (binding 0 uniform, bindings 1-5 `storage_entry(_, true)`, binding 6 `storage_entry(6, false)`), and pass `&entries`. Keep the descriptor construction in one place. Add `pub mod kernel;` to `lib.rs`.

- [ ] **Step 6: Run — verify pass (or clean skip)**

Run: `cargo test -p geodesic-gpu --test forces_gpu_vs_cpu -- --nocapture`
Expected: PASS on a machine with an adapter (WARP/lavapipe count); clean skip otherwise. If forces disagree, the bug is in the shader math or CSR/param upload — fix at the root (do not loosen `tol`).

- [ ] **Step 7: clippy + commit**

```bash
cargo clippy -p geodesic-gpu --all-targets -- -D warnings
git add geodesic-gpu/src/nonbonded.wgsl geodesic-gpu/src/kernel.rs geodesic-gpu/src/lib.rs geodesic-gpu/tests/forces_gpu_vs_cpu.rs geodesic-gpu/tests/common
git commit -m "feat(gpu): WGSL non-bonded LJ kernel and f32 force evaluator, matched to CPU"
git push
```

---

### Task 7: `GpuBackend` implementing `ComputeBackend`

**Files:**
- Create: `geodesic-gpu/src/gpu_backend.rs`
- Modify: `geodesic-gpu/src/lib.rs` (`pub mod gpu_backend;`)
- Test: extend `geodesic-gpu/tests/forces_gpu_vs_cpu.rs` with an `ala_dipeptide` full-backend case

**Interfaces:**
- Consumes: `GpuContext`, `NonbondedKernel`, `build_csr`, engine `neighbor`/`force::bonded`/`integrator`.
- Produces: `pub struct GpuBackend`; `pub fn try_new(atoms: AtomData, topology: BondedTopology, params: &SimParams) -> Result<GpuBackend, BackendError>`; full `impl ComputeBackend for GpuBackend`.

- [ ] **Step 1: Write the failing test (full backend vs CPU on ala_dipeptide)**

Append to `geodesic-gpu/tests/forces_gpu_vs_cpu.rs` (uses the `common` module already imported at the top of the file):

```rust
#[test]
fn ala_dipeptide_full_backend_matches_cpu() {
    use common::{clone_positions, load, params};
    use geodesic_core::ComputeBackend;
    let Some(_ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load("ala_dipeptide");
    let (state2, atoms2, topology2) = load("ala_dipeptide");
    let p = params(state.pos_x.len());
    let n = state.pos_x.len();

    let mut s1 = clone_positions(&state);
    let mut s2 = clone_positions(&state2);

    let mut cpu = geodesic_engine::cpu_backend::CpuBackend::new(atoms, topology, &p);
    let mut gpu = geodesic_gpu::gpu_backend::GpuBackend::try_new(atoms2, topology2, &p).unwrap();

    cpu.build_neighbor_list(&mut s1, &p);
    let fc = cpu.compute_forces(&s1).clone();
    gpu.build_neighbor_list(&mut s2, &p);
    let fg = gpu.compute_forces(&s2).clone();

    for i in 0..n {
        for (cpuv, gpuv) in [(fc.fx[i], fg.fx[i]), (fc.fy[i], fg.fy[i]), (fc.fz[i], fg.fz[i])] {
            let diff = (cpuv - gpuv).abs();
            let bound = 1e-4 * cpuv.abs().max(1.0);
            assert!(diff <= bound, "atom {i}: cpu={cpuv}, gpu={gpuv}, diff={diff}");
        }
    }
}
```

Note: `geodesic-gpu` needs `geodesic-engine` as a dev-dependency for `CpuBackend` in this test — add `geodesic-engine = { workspace = true }` under `[dev-dependencies]` in `geodesic-gpu/Cargo.toml` (it is already a normal dependency, so this is only needed if the normal dep is not visible to tests; normal deps are visible to integration tests, so no change is required — verify at build time).

- [ ] **Step 2: Run — verify it fails**

Run: `cargo test -p geodesic-gpu --test forces_gpu_vs_cpu ala_dipeptide_full_backend_matches_cpu`
Expected: FAIL.

- [ ] **Step 3: Implement `gpu_backend.rs`**

`geodesic-gpu/src/gpu_backend.rs`:

```rust
use crate::device::{self, GpuContext};
use crate::kernel::{NonbondedInput, NonbondedKernel};
use crate::neighbor_csr::build_csr;
use geodesic_core::{
    AtomData, BackendError, BondedTopology, ComputeBackend, ConvergenceError, ForceBuffer,
    NeighborList, SimParams, SimState,
};
use geodesic_engine::force::bonded;
use geodesic_engine::{integrator, neighbor};

pub struct GpuBackend {
    ctx: GpuContext,
    kernel: NonbondedKernel,
    atoms: AtomData,
    topology: BondedTopology,
    neighbor_list: NeighborList,
    offsets: Vec<u32>,
    neighbors: Vec<u32>,
    box_size: [f64; 3],
    max_constr_iter: u32,
    constr_tol: f64,
    reduced: ForceBuffer,
    potential_energy: f64,
}

impl GpuBackend {
    pub fn try_new(
        atoms: AtomData,
        topology: BondedTopology,
        params: &SimParams,
    ) -> Result<Self, BackendError> {
        let ctx = device::try_new()?;
        let kernel = NonbondedKernel::new(&ctx)?;
        let n = atoms.mass.len();
        let neighbor_list = NeighborList {
            pair_i: Vec::new(),
            pair_j: Vec::new(),
            ref_x: vec![0.0; n],
            ref_y: vec![0.0; n],
            ref_z: vec![0.0; n],
            r_cutoff: params.r_cutoff,
            r_skin: params.r_skin,
            r_switch: params.r_switch,
        };
        Ok(Self {
            ctx,
            kernel,
            atoms,
            topology,
            neighbor_list,
            offsets: vec![0; n + 1],
            neighbors: Vec::new(),
            box_size: params.box_size,
            max_constr_iter: params.max_constr_iter,
            constr_tol: params.constr_tol,
            reduced: ForceBuffer { fx: vec![0.0; n], fy: vec![0.0; n], fz: vec![0.0; n] },
            potential_energy: 0.0,
        })
    }
}

impl ComputeBackend for GpuBackend {
    fn build_neighbor_list(&mut self, state: &mut SimState, params: &SimParams) {
        self.neighbor_list = neighbor::build(state, params, &self.topology);
        let n = state.pos_x.len();
        let (off, nbr) = build_csr(&self.neighbor_list.pair_i, &self.neighbor_list.pair_j, n);
        self.offsets = off;
        self.neighbors = nbr;
    }

    fn compute_forces(&mut self, state: &SimState) -> &ForceBuffer {
        let n = state.pos_x.len();
        self.reduced.fx.iter_mut().for_each(|x| *x = 0.0);
        self.reduced.fy.iter_mut().for_each(|x| *x = 0.0);
        self.reduced.fz.iter_mut().for_each(|x| *x = 0.0);
        let mut potential = 0.0;
        potential += bonded::compute_bond_forces(state, &self.topology, &mut self.reduced.fx, &mut self.reduced.fy, &mut self.reduced.fz);
        potential += bonded::compute_angle_forces(state, &self.topology, &mut self.reduced.fx, &mut self.reduced.fy, &mut self.reduced.fz);
        potential += bonded::compute_dihedral_forces(state, &self.topology, &mut self.reduced.fx, &mut self.reduced.fy, &mut self.reduced.fz);

        let input = NonbondedInput {
            pos_x: &state.pos_x,
            pos_y: &state.pos_y,
            pos_z: &state.pos_z,
            sigma: &self.atoms.sigma,
            epsilon: &self.atoms.epsilon,
            offsets: &self.offsets,
            neighbors: &self.neighbors,
            r_cutoff: self.neighbor_list.r_cutoff,
            r_switch: self.neighbor_list.r_switch,
            box_size: self.box_size,
        };
        let (gpu_f, nb_energy) = self.kernel.evaluate(&self.ctx, &input);
        for i in 0..n {
            self.reduced.fx[i] += gpu_f[i][0] as f64;
            self.reduced.fy[i] += gpu_f[i][1] as f64;
            self.reduced.fz[i] += gpu_f[i][2] as f64;
        }
        potential += nb_energy as f64;
        self.potential_energy = potential;
        &self.reduced
    }

    fn geodesic_drift(&mut self, state: &mut SimState, dt: f64) -> Result<(), ConvergenceError> {
        integrator::drift_and_constrain(state, &self.topology, &self.atoms, dt, self.max_constr_iter, self.constr_tol)
    }

    fn reduce_forces(&self) -> ForceBuffer {
        self.reduced.clone()
    }

    fn potential_energy(&self) -> f64 { self.potential_energy }
    fn atoms(&self) -> &AtomData { &self.atoms }
    fn topology(&self) -> &BondedTopology { &self.topology }
    fn needs_rebuild(&self, state: &SimState) -> bool {
        neighbor::needs_rebuild(state, &self.neighbor_list)
    }
    fn n_threads(&self) -> usize { 1 }
}
```

Note: this requires `geodesic_engine::force::bonded`, `neighbor`, and `integrator` to be reachable (they are `pub`). If `force` is re-exported differently, use the actual path (`geodesic_engine::force::bonded`).

- [ ] **Step 4: Wire the test body, run — verify pass (or skip)**

Run: `cargo test -p geodesic-gpu --test forces_gpu_vs_cpu -- --nocapture`
Expected: all three cases PASS where an adapter exists (ala_dipeptide validates exclusions through CSR); clean skip otherwise.

- [ ] **Step 5: clippy + commit**

```bash
cargo clippy -p geodesic-gpu --all-targets -- -D warnings
git add geodesic-gpu/src/gpu_backend.rs geodesic-gpu/src/lib.rs geodesic-gpu/tests/forces_gpu_vs_cpu.rs
git commit -m "feat(gpu): GpuBackend (GPU non-bonded + CPU bonded/constraint), matches CPU on ala_dipeptide"
git push
```

---

### Task 8: GPU determinism test

**Files:**
- Test: `geodesic-gpu/tests/gpu_determinism.rs`

**Interfaces:**
- Consumes: `NonbondedKernel::evaluate`.

- [ ] **Step 1: Write the test**

`geodesic-gpu/tests/gpu_determinism.rs`:

```rust
mod common;
use common::{cpu_nonbonded_reference, load, params};
use geodesic_gpu::device;
use geodesic_gpu::kernel::{NonbondedInput, NonbondedKernel};
use geodesic_gpu::neighbor_csr::build_csr;

#[test]
fn two_gpu_evaluations_are_bit_identical() {
    let Some(ctx) = device::context_or_skip() else { return };
    let (state, atoms, topology) = load("water_box_4");
    let p = params(state.pos_x.len());
    let (s, list, _fx, _fy, _fz) = cpu_nonbonded_reference(&state, &atoms, &topology, &p);
    let n = s.pos_x.len();
    let (offsets, neighbors) = build_csr(&list.pair_i, &list.pair_j, n);
    let kernel = NonbondedKernel::new(&ctx).unwrap();
    let input = NonbondedInput {
        pos_x: &s.pos_x,
        pos_y: &s.pos_y,
        pos_z: &s.pos_z,
        sigma: &atoms.sigma,
        epsilon: &atoms.epsilon,
        offsets: &offsets,
        neighbors: &neighbors,
        r_cutoff: list.r_cutoff,
        r_switch: list.r_switch,
        box_size: p.box_size,
    };
    let (a, ea) = kernel.evaluate(&ctx, &input);
    let (b, eb) = kernel.evaluate(&ctx, &input);
    assert_eq!(a, b);
    assert_eq!(ea.to_bits(), eb.to_bits());
}
```

- [ ] **Step 2: Run — verify it fails, then passes**

Run: `cargo test -p geodesic-gpu --test gpu_determinism`
Expected: FAIL to compile only if `common` isn't present (it is, from Task 6). With an adapter, PASS (bit-identical `Vec<[f32;3]>` and identical energy bits); clean skip otherwise.

- [ ] **Step 3: Run again to confirm**

Run: `cargo test -p geodesic-gpu --test gpu_determinism -- --nocapture`
Expected: PASS or clean skip.

- [ ] **Step 4: clippy + commit**

```bash
cargo clippy -p geodesic-gpu --all-targets -- -D warnings
git add geodesic-gpu/tests/gpu_determinism.rs
git commit -m "test(gpu): two evaluations on the same adapter are bit-identical (f32)"
git push
```

---

### Task 9: Wire `GpuBackend` into the binary behind the `gpu` feature

**Files:**
- Modify: `geodesic/src/lib.rs:186-192` (config-error site) and `:233` (backend construction)

**Interfaces:**
- Consumes: `geodesic_gpu::gpu_backend::GpuBackend`.

- [ ] **Step 1: Write a failing gated test**

Add to `geodesic/tests/` a new `gpu_run.rs` (compiled only with the feature). This mirrors `geodesic/tests/determinism.rs:10-65`'s config setup, changing only `backend = "gpu"`:

```rust
#![cfg(feature = "gpu")]

use std::path::{Path, PathBuf};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("geodesic-engine").join("tests").join("fixtures")
}
fn slashed(p: &Path) -> String {
    p.to_string_lossy().replace('\\', "/")
}

#[test]
fn gpu_backend_runs_ala_dipeptide_or_skips() {
    let fixtures = fixtures_dir();
    let prmtop = slashed(&fixtures.join("ala_dipeptide.prmtop"));
    let inpcrd = slashed(&fixtures.join("ala_dipeptide.inpcrd"));
    let out_dir = std::env::temp_dir();
    let dcd = slashed(&out_dir.join("geodesic_gpu_run.dcd"));
    let csv = slashed(&out_dir.join("geodesic_gpu_run.csv"));
    let cfg_path = out_dir.join("geodesic_gpu_run.toml");
    let config = format!(
        r#"
[run]
n_steps        = 20
frame_interval = 5
backend        = "gpu"
n_threads      = 1

[system]
prmtop   = "{prmtop}"
inpcrd   = "{inpcrd}"
box_size = [1000.0, 1000.0, 1000.0]
periodic = false

[integrator]
dt           = 0.001
total_energy = 5.12

[nonbonded]
r_cutoff = 12.0
r_skin   = 14.0
r_switch = 10.0

[constraints]
max_iter  = 100
tolerance = 1.0e-8

[output]
trajectory = "{dcd}"
energy_log = "{csv}"
"#
    );
    std::fs::write(&cfg_path, config).unwrap();

    match geodesic::run_from_config_file(&cfg_path) {
        Ok(summary) => {
            let bytes = std::fs::read(&summary.trajectory).unwrap();
            assert!(!bytes.is_empty(), "GPU run produced an empty DCD");
        }
        Err(geodesic_core::SimError::Backend(geodesic_core::BackendError::NoAdapter)) => {
            eprintln!("skipping GPU run test: no adapter available");
        }
        Err(e) => panic!("GPU run failed: {e}"),
    }
}
```

Note: confirm the run summary's trajectory-path field name against `determinism.rs` (`summary.trajectory`) and that `run_from_config_file` returns `Result<_, geodesic_core::SimError>`. `geodesic-core` must be a dev-dependency of `geodesic` for the error match — add `geodesic-core = { workspace = true }` to `[dev-dependencies]` in `geodesic/Cargo.toml` if it is not already a normal dependency (it is a normal dependency, visible to tests, so likely no change needed).

- [ ] **Step 2: Run — verify it fails**

Run: `cargo test -p geodesic --features gpu --test gpu_run`
Expected: FAIL.

- [ ] **Step 3: Gate the backend selection**

In `geodesic/src/lib.rs`, replace the hard config error and the `CpuBackend` construction with feature-aware selection. Current code errors for any non-Cpu backend (line ~188); change to:

```rust
    let mut backend: Box<dyn ComputeBackend> = match config.backend {
        Backend::Cpu => Box::new(CpuBackend::new(atoms, topology, &params)),
        #[cfg(feature = "gpu")]
        Backend::Gpu => Box::new(geodesic_gpu::gpu_backend::GpuBackend::try_new(atoms, topology, &params)?),
        #[cfg(not(feature = "gpu"))]
        Backend::Gpu => {
            return Err(ConfigError::InvalidValue {
                key: "run.backend".to_string(),
                value: "gpu".to_string(),
                reason: "this build was compiled without the `gpu` feature; rebuild with --features gpu".to_string(),
            }.into())
        }
        Backend::Hybrid => {
            return Err(ConfigError::InvalidValue {
                key: "run.backend".to_string(),
                value: "hybrid".to_string(),
                reason: "the hybrid backend lands in v0.6".to_string(),
            }.into())
        }
    };
```

`GpuBackend::try_new` returns `BackendError`, which converts into `SimError` via the existing `#[from]`. Confirm `run_from_config_file` returns `Result<_, SimError>`; if it returns a narrower error, add the `From`/`?` plumbing.

- [ ] **Step 4: Run — verify pass (or skip)**

Run: `cargo test -p geodesic --features gpu --test gpu_run -- --nocapture`
Expected: PASS where an adapter exists (DCD written); clean skip on `NoAdapter`.
Run: `cargo test -p geodesic` (no feature) — the existing "gpu backend rejected" behavior still holds via the `cfg(not)` arm.

- [ ] **Step 5: clippy both configs + commit**

```bash
cargo clippy -p geodesic --all-targets -- -D warnings
cargo clippy -p geodesic --all-targets --features gpu -- -D warnings
git add geodesic/src/lib.rs geodesic/tests/gpu_run.rs
git commit -m "feat(gpu): select GpuBackend for run.backend=gpu under the gpu feature"
git push
```

---

### Task 10: CI, ROADMAP, and handoff

**Files:**
- Modify: `.github/workflows/ci.yml`
- Modify: `ROADMAP.md:177-187` (v0.5 checkboxes/exit criteria)
- Modify: `memory.md`

**Interfaces:** none (project bookkeeping).

- [ ] **Step 1: Add a Linux job and gpu-build steps**

In `.github/workflows/ci.yml`, add an `ubuntu-latest` job mirroring the existing windows job, and to both jobs add after the default build/test:

```yaml
      - name: build gpu feature
        run: cargo build -p geodesic --features gpu
      - name: clippy gpu feature
        run: cargo clippy -p geodesic-gpu --all-targets -- -D warnings
      - name: gpu tests (skip if no adapter)
        run: cargo test -p geodesic-gpu -- --nocapture
```

The GPU tests self-skip when no adapter (WARP on the Windows runner may or may not initialize; the tests tolerate absence). Do not make `--features gpu` build the only build — keep the default CPU build/test steps unchanged so M1 coverage is untouched.

- [ ] **Step 2: Run the workspace suite locally (both configs)**

Run: `cargo test --workspace`
Run: `cargo build --features gpu` (from the `geodesic` crate context) and `cargo test -p geodesic-gpu`
Expected: green (GPU tests pass or skip).

- [ ] **Step 3: Tick ROADMAP v0.5 and note the exit criteria met**

In `ROADMAP.md`, mark the v0.5 items done and record: `cargo build --features gpu` green on Windows+Linux; GPU/CPU force agreement at f32 tol on lj_pair/water_box_4/ala_dipeptide; same-adapter determinism bit-identical. Note the three SAD deviations point to ADR 0001-0004.

- [ ] **Step 4: Update `memory.md` handoff**

Set status to "v0.5 shipped (GPU backend, M2)"; next priorities = v0.6 hybrid (position/velocity residency, min-image in bonded terms for periodic systems, GPU-side constraint reduce). Record: GPU is f32 (ADR 0002); the CPU pair set is uploaded as a CSR gather list (ADR 0003); geodesic-gpu depends on geodesic-engine (ADR 0004); platform is Windows+Linux (ADR 0001).

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/ci.yml ROADMAP.md memory.md
git commit -m "ci(gpu): build+clippy --features gpu on windows+linux; close out v0.5"
git push
```

---

## Notes for the implementer

- **wgpu API drift:** the code targets wgpu 22 per SAD §14. Several wgpu APIs used here (`request_adapter` return type, `entry_point: Option<&str>`, `create_pipeline` `cache`/`compilation_options` fields, `device.poll` signature) changed across recent versions. If the resolved version rejects a call, adapt to that version's signature — the algorithm is unchanged. Record the actual version in ADR 0001.
- **Never loosen the f32 tolerance to force a green.** If ala_dipeptide forces don't agree at ~1e-4 relative, the fault is in the CSR expansion, the exclusion set, the `min_image`/switch math, or an upload layout mismatch — fix the root (project CLAUDE.md).
- **Fixture loading:** copy the exact loader calls from `geodesic-engine/tests/fixture_gradient_check.rs` so the CPU reference matches bit-for-bit what the engine computes.
- **Release:** after Task 10 is green, v0.5 meets its ROADMAP exit criteria; publishing the `v0.5` GitHub release (per project CLAUDE.md standing authorization) is a separate manual step, written in Felipe's voice — not part of this plan's automated steps.

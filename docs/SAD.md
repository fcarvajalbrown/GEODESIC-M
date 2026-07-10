## Software Architecture Document: GEODESIC-M

### 1. System Overview
**GEODESIC-M** is a high-precision simulation environment designed to explore protein conformational manifolds through classical mechanics. Unlike probabilistic LLM-based predictors, this system provides a deterministic, physics-based mapping of the high-dimensional energy landscape of a protein sequence. The system's distinguishing claim is that conformational transitions are geodesics on a Riemannian manifold defined by the Jacobi metric (§2.0), not Newtonian paths through flat space.

---

### 2. Mathematical Foundation

The system operates in the classical mechanics regime of the **Karplus-Levitt-Warshel** framework, extended by Riemannian geometry and topological data analysis. The central claim of GEODESIC-M is that protein conformational transitions are not Newtonian paths through flat $\mathbb{R}^{3N}$ — they are **geodesics on a curved manifold** whose geometry is determined by the energy landscape itself.

---

#### 2.0 The Conformational Manifold (Jacobi–Maupertuis Principle)

On the energy hypersurface $H = E$, the **Jacobi metric** on configuration space $\mathcal{M} \subset \mathbb{R}^{3N}$ is:

$$g^J_{ij}(x) = 2\bigl(E - V(x)\bigr)\, m_i\, \delta_{ij}$$

By the **Maupertuis principle**, classical trajectories with total energy $E$ are exactly the **geodesics** of $(\mathcal{M},\, g^J)$. The curvature of $g^J$ encodes the topology of the energy landscape: regions of high $V$ contract the metric (barriers become narrow necks), and regions of low $V$ expand it (minima become wide basins). The "folding funnel" is a geometric object — a region of positive Ricci scalar curvature funneling geodesics toward the native state.

> This is the operational meaning of "manifold" in GEODESIC-M: the metric, not just the potential, *is* the physics.

Empirical validation: Diepeveen et al. (*PNAS* 2024, arXiv:2308.07818) demonstrated that geodesics on such manifolds recover MD trajectory statistics for protein conformational transitions at a fraction of the simulation cost.

---

#### 2.1 Total Potential Energy

$$V_{\text{total}} = V_{\text{bonded}} + V_{\text{non-bonded}}$$

#### 2.1.1 Bonded Interactions

| Term | Equation | Parameters |
| :--- | :--- | :--- |
| Bond stretching | $V_b = k_b(r - r_0)^2$ | $k_b$: force constant, $r_0$: equilibrium length |
| Angle bending | $V_\theta = k_\theta(\theta - \theta_0)^2$ | $k_\theta$: force constant, $\theta_0$: equilibrium angle |
| Dihedral torsion | $V_\phi = k_\phi\bigl[1 + \cos(n\phi - \delta)\bigr]$ | $k_\phi$: barrier height, $n$: multiplicity, $\delta$: phase |

Parameters are sourced from a standard force field parameter file (AMBER or CHARMM format — resolved in §I/O section).

#### 2.1.2 Non-Bonded Interactions

**Lennard-Jones (van der Waals):**

$$V_{LJ}(r) = 4\epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^6 \right]$$

Truncated at cutoff $r_c$ with a smooth switching function over $[r_{\text{sw}},\, r_c]$ to avoid force discontinuities at the boundary.

---

#### 2.2 Equations of Motion

Newton’s second law in the ambient space:

$$F_i = -\nabla_i V_{\text{total}} = m_i \ddot{r}_i$$

When holonomic constraints are active (fixed bond lengths), the dynamics are restricted to the constraint submanifold $\mathcal{C} \subset \mathbb{R}^{3N}$. Integration must respect the geometry of $\mathcal{C}$, not merely reproject into it after each step.

---

#### 2.3 Geodesic BAOAB Integration

Standard Velocity Verlet linearizes constraints and reprojects (SHAKE/RATTLE) — it approximates the geodesic. GEODESIC-M uses the **Geodesic BAOAB** integrator (Leimkuhler & Matthews, *Proc. R. Soc. A* 2016), which computes the **A step as a true geodesic segment on $\mathcal{C}$** via the exponential map, not a linearized drift. For the NVE (deterministic) ensemble, the O (thermostat) step is omitted, yielding the symplectic **BAB** splitting:

$$\underbrace{v \leftarrow v + \tfrac{\Delta t}{2}\,M^{-1}F}_{\textbf{B: force kick}} \;\longrightarrow\; \underbrace{r \leftarrow \exp_r\!\bigl(v\,\Delta t\bigr)\big|_{\mathcal{C}}}_{\textbf{A: geodesic drift on } \mathcal{C}} \;\longrightarrow\; \underbrace{v \leftarrow v + \tfrac{\Delta t}{2}\,M^{-1}F}_{\textbf{B: force kick}}$$

The exponential map $\exp_r$ advances the position along the geodesic of $\mathcal{C}$ initiating at $r$ with velocity $v$, computed iteratively via the constraint Jacobian. Practical consequence: this scheme supports timesteps of **8–9 fs** (vs. 2 fs for Velocity Verlet) with no loss in energy conservation, because it never leaves the constraint manifold and does not accumulate reprojection drift.

> **Key sequencing constraint:** The force at $t + \Delta t$ must be fully evaluated — including neighbor list validity check — before the second B half-step. This is identical to standard Verlet and does not introduce new pipeline dependencies.

---

#### 2.4 Periodic Boundary Conditions (PBC)

For a simulation box of side $L$, the minimum image convention gives the effective pairwise displacement:

$$r_{ij}^* = r_{ij} - L \cdot \operatorname{round}\!\left(\frac{r_{ij}}{L}\right)$$

All pairwise distances use $r_{ij}^*$. Atoms leaving the box are wrapped to $[0,\, L)$ before each neighbor list rebuild.

---

#### 2.5 Verlet Neighbor Lists

Naive pairwise evaluation is $O(N^2)$. A **Verlet neighbor list** with cutoff $r_c$ and skin $r_s > r_c$ reduces amortized cost to $O(N)$:

- List rebuilt when any atom displaces more than $\tfrac{r_s - r_c}{2}$ since the last build.
- Only pairs within the list are evaluated; pairs with $r > r_c$ contribute zero force.

---

#### 2.6 Persistent Sheaf Laplacian (PSL) — Flexibility Analysis

Standard normal mode analysis (GNM/ANM) uses a scalar spring constant per contact. GEODESIC-M replaces this with the **Persistent Sheaf Laplacian** (Hayes et al., *J. Phys. Chem. B* 2025), a sheaf-theoretic generalization of the Hodge Laplacian that captures correlated multi-scale flexibility.

A sheaf $\mathcal{F}$ over the protein’s atomic contact complex assigns a vector space $\mathcal{F}(\sigma)$ of local observables to each simplex $\sigma$ (atom, bond, contact triangle), with restriction maps encoding inter-simplex consistency. The **sheaf Laplacian** $\mathcal{L}^k_{\mathcal{F}}$ is:

$$\mathcal{L}^k_{\mathcal{F}} = \bigl(\delta^k_{\mathcal{F}}\bigr)^T \delta^k_{\mathcal{F}} + \delta^{k-1}_{\mathcal{F}}\bigl(\delta^{k-1}_{\mathcal{F}}\bigr)^T$$

where $\delta^k_{\mathcal{F}}$ is the coboundary map twisted by the sheaf restriction maps. Running this over a filtration of growing contact radius $r$ gives **persistent** spectral invariants of the flexibility landscape. Validated result: PSL achieves 32% higher B-factor prediction accuracy than GNM on a 364-protein benchmark.

---

#### 2.7 Zigzag Persistence — Trajectory Topology

Standard persistent homology requires a monotone filtration and cannot track topological features that appear, disappear, and re-appear. **Zigzag persistence** (Carlsson & de Silva) relaxes this to bidirectional inclusion sequences:

$$\emptyset \leftarrow K_0 \rightarrow K_1 \leftarrow K_2 \rightarrow K_3 \leftarrow \cdots$$

Applied to GEODESIC-M trajectory snapshots, each $K_t$ is the atomic contact complex at frame $t$. Zigzag persistence produces a **barcode** over time: a rigorous, coordinate-free record of when topological features ($H_0$ connectivity, $H_1$ loops, $H_2$ voids) are born and die across the simulation. This replaces RMSD-based trajectory clustering with a topologically invariant signature of the folding/unfolding event.

---

### 3. Component Architecture

| Component | Crate | Responsibility | Architectural Pattern |
| :--- | :--- | :--- | :--- |
| **I/O Layer** | `geodesic-io` | Parses `prmtop` + `inpcrd` + `config.toml`; writes DCD trajectory, CSV energy log, JSON barcode, PDB snapshots | Parser pipeline |
| **Force Field Engine** | `geodesic-engine::force` | Evaluates bonded (bonds, angles, dihedrals) and non-bonded (LJ + switching) forces; SoA layout for SIMD | Data-Oriented Design |
| **Geodesic BAOAB Integrator** | `geodesic-engine::integrator` | BAB time-stepping (4–9 fs) with true geodesic drift on constraint manifold $\mathcal{C}$ via exponential map | Symplectic Numerical Solver |
| **Compute Backend** | `ComputeBackend` trait | Dispatches force evaluation and geodesic drift to CPU (Rayon + SIMD), GPU (wgpu), or hybrid; selected at startup | Strategy Pattern |
| **Topology Pipeline** | `geodesic-topo` | PSL flexibility analysis (§2.6) + Zigzag persistence barcode (§2.7) from completed trajectory | Post-processing Pipeline |
| **GUI Renderer** | `geodesic-gui` | wgpu 3D atom/bond viewer; ring buffer consumer; data export panel (M4+) | Event-driven Observer |

---

### 4. Data Flow

**Startup:**
1. `config.toml` → `SimParams` (timestep, cutoffs, backend selection, output paths)
2. `protein.prmtop` → `AtomData` (masses, LJ parameters) + `BondedTopology` (bonds, angles, dihedrals, constraints)
3. `protein.inpcrd` → initial `SimState` (positions, velocities)

**Integration loop** (repeated for each step):
1. Check if neighbor list needs rebuild; if so, rebuild Verlet pair list (§2.5)
2. Evaluate bonded + non-bonded forces via `ComputeBackend::compute_forces()` → `SimState::force_{x,y,z}`
3. NaN guard: scan force arrays; abort with `NumericalError` on first NaN (§12.3)
4. Geodesic BAB: B half-kick → geodesic drift on $\mathcal{C}$ (constraint solve) → B half-kick (§2.3)
5. Energy drift check; warn or abort per config (§12.4)
6. Every `frame_interval` steps: write `TrajectoryFrame` to ring buffer (non-blocking)

**Concurrent output** (runs alongside integration):
- Ring buffer → DCD writer appends frame to `output.dcd`
- Ring buffer → GUI renderer updates 3D view (M4+)
- Energy values → `energy.csv` appended line-by-line

**Post-processing** (after integration completes):
1. Completed trajectory → PSL contact complex filtration → eigenvalue solve → per-atom flexibility scores
2. Trajectory snapshots → Zigzag persistence → `barcode.json`
3. Flexibility scores backfilled into a PDB snapshot (B-factor column)

---

### 5. Non-Functional Requirements (NFRs)

*   **Determinism:** Bit-for-bit identical output given the same binary, hardware, initial coordinates, and timestep. Floating-point operations must follow IEEE 754 strict ordering; no non-deterministic parallelism (e.g., unordered floating-point reductions).
*   **Auditability:** Every force contribution must be attributable to a specific atom pair and interaction type, queryable by step index. No implicit accumulation or fused operations that obscure intermediate values.
*   **Performance:** Near-linear thread scaling up to the available core count (CPU backend). Force calculation cost reduced from $O(N^2)$ to $O(N)$ amortized via Verlet neighbor lists (§2.5). Concrete ns/day throughput targets benchmarked and tracked per §13.9.

---

### 6. Implementation Constraints
*   **Language:** **Rust** — mandatory for memory safety, zero-cost abstractions, and deterministic floating-point behavior without a garbage collector.
*   **Hardware:** CPU backend uses SIMD (AVX2/AVX-512) via SoA memory layout; GPU backend uses wgpu compute shaders (M2+). Neighbor lists reduce the non-bonded force calculation from $O(N^2)$ to $O(N)$ amortized, making CPU competitive at research-scale atom counts. See §7 for the full backend selection and milestone plan.

---

### 7. Concurrency & Parallelism Architecture

The compute layer is abstracted behind a `ComputeBackend` trait. The user selects the backend at startup; the simulation engine is backend-agnostic. A future GUI layer connects to the simulation via a trajectory channel and shares the wgpu device with the GPU backend, avoiding redundant GPU context creation.

**Milestone plan:**

| Milestone | Backend | GUI |
| :--- | :--- | :--- |
| **M1 (current)** | CPU only | None — file output only |
| **M2** | GPU (wgpu compute shaders) | None |
| **M3** | Hybrid CPU+GPU | None |
| **M4** | Any backend | 3D viewer + live streaming |
| **M5** | Any backend | Data export UI |

---

#### 7.1 Backend Abstraction

All compute work is dispatched through a single trait:

```rust
trait ComputeBackend: Send {
    fn build_neighbor_list(&mut self, state: &SimState, params: &SimParams);
    fn compute_forces(&mut self, state: &SimState) -> &ForceBuffer;
    fn geodesic_drift(&mut self, state: &mut SimState, dt: f64) -> Result<(), ConvergenceError>;
    fn reduce_forces(&self) -> ForceBuffer;
}
```

The simulation loop holds a `Box<dyn ComputeBackend>` and never inspects its concrete type. Selecting CPU, GPU, or hybrid is a startup decision, not scattered throughout the codebase.

---

#### 7.2 Milestone 1 — CPU Backend

**Thread pool:** A fixed Rayon thread pool of size $T$ = physical core count (not hyperthreads), created once at startup. $T$ is recorded in the simulation log — it is a reproducibility parameter.

**Static force decomposition (determinism):** Rayon's work-stealing scheduler produces non-deterministic floating-point reduction order, violating the Determinism NFR. Resolution:

- The $N \times N$ pair interaction space is partitioned into $T$ fixed horizontal strips: thread $k$ evaluates all pairs $(i,j)$ where $i \in \bigl[k \cdot \lfloor N/T \rfloor,\ (k+1) \cdot \lfloor N/T \rfloor\bigr)$.
- Each thread accumulates forces into a **private per-thread buffer** $\mathbf{f}^{(k)} \in \mathbb{R}^{N \times 3}$.
- The master thread reduces in thread-index order: $\mathbf{F} = \sum_{k=0}^{T-1} \mathbf{f}^{(k)}$ — fixed sequential sum, reproducible by construction.

**Memory cost:** $O(N \cdot T \cdot 24)$ bytes. At $N = 10^4$, $T = 16$: ~23 MB. At $N = 10^5$, $T = 16$: ~230 MB.

**SIMD:** SoA memory layout (mandated — see §8) allows the inner $j$-loop to auto-vectorize under LLVM AVX2/AVX-512. Hot loop functions carry `#[target_feature(enable = "avx2")]` with scalar fallback via `std::arch::is_x86_feature_detected!`.

**Constraint solver:** Within each Geodesic BAOAB iteration, constraint forces are dispatched via `rayon::par_iter`. The convergence reduce (global max of $|\lambda_i|$) is performed sequentially over per-thread partial maxima in thread-index order — deterministic. Non-convergence after $I_{\max}$ iterations is a hard `ConvergenceError`, never silent degradation.

---

#### 7.3 Milestone 2 — GPU Backend (wgpu)

**Why wgpu over CUDA:** wgpu targets Vulkan/Metal/DX12 — cross-platform, and crucially the same API used for 3D rendering in the future GUI. The GPU backend and GUI renderer **share a single `wgpu::Device`**, eliminating the cost of a second GPU context and enabling zero-copy frame streaming from simulation to renderer.

**Force calculation:** The $O(N^2)$ non-bonded force loop is mapped to a compute shader dispatched as $\lceil N/64 \rceil$ workgroups of 64 threads. Each workgroup tile loads a block of positions into shared memory (tile size 64) and sweeps all tiles — the standard GPU tiled force evaluation pattern.

**Determinism on GPU:** Warp-level reduction order is fixed by using a tree reduction with a defined evaluation order, not `atomicAdd`. This makes results reproducible on the same GPU model and driver version. Cross-GPU reproducibility is not guaranteed and is documented as a known limitation.

**Constraint solver on GPU:** The convergence reduce (max $|\lambda_i|$) is performed as a GPU tree reduction — no CPU round-trip per iteration. Result is read back once per step after convergence.

---

#### 7.4 Milestone 3 — Hybrid Backend

CPU and GPU divide responsibilities by workload type:

| Workload | Backend | Reason |
| :--- | :--- | :--- |
| Non-bonded force calculation | GPU | $O(N^2)$, embarrassingly parallel |
| Bonded forces (bonds, angles, dihedrals) | CPU | Low $N$, dependency graph limits GPU occupancy |
| Geodesic BAOAB constraint solve | CPU | Iterative, convergence-dependent; CPU overhead acceptable |
| Neighbor list rebuild | CPU | Irregular memory access, hard to vectorize on GPU |
| PSL + Zigzag pipeline | CPU | Sparse graph — poor GPU fit |

Positions and velocities live on GPU during force computation; copied to CPU before the constraint solve step.

---

#### 7.5 GUI Integration Architecture

The GUI runs in a dedicated OS thread, completely decoupled from the simulation loop. Communication is one-directional:

```
Simulation thread  ──[ring buffer: TrajectoryFrame]──►  GUI thread
```

- The ring buffer holds a configurable number of frames (default 256). If the GUI falls behind, the oldest frames are overwritten — the simulation never blocks.
- `TrajectoryFrame` contains: atomic positions, per-atom PSL flexibility scores, current energy, step index, and wall-clock time.
- The GUI renders 3D atomic positions using wgpu (shared device with the GPU backend where available). Bonds are rendered as cylinders; atoms as instanced spheres colored by element or flexibility score.
- **Data export** (available from the GUI): PDB snapshot, DCD trajectory segment, CSV of energy/RMSD over time, and the Zigzag persistence barcode as JSON.
- `OPENBLAS_NUM_THREADS` is pinned to 1 at process startup, regardless of backend, to keep the PSL eigenvalue solve deterministic.

---

#### 7.6 Topology Pipeline (PSL + Zigzag)

Post-processing — not on the integration hot path. Runs after trajectory collection, always on CPU.

| Stage | Parallelism | Notes |
| :--- | :--- | :--- |
| PSL contact complex construction | Rayon `par_iter` over filtration steps | Embarrassingly parallel per radius $r$ |
| PSL eigenvalue solve | Single-threaded LAPACK (`ndarray-linalg`) | `OPENBLAS_NUM_THREADS=1` for determinism |
| Zigzag persistence | Sequential | Boundary matrix reduction; Ripser (C++ FFI) or native Rust in v1 |

---

### 8. Data Structures

**Layout rule:** SoA (Structure of Arrays) for any data touched in the integration hot loop. AoS (Array of Structs) is acceptable only for cold data read once per frame (metadata, parameters). This is not a style preference — it is required for auto-vectorization (§7.2) and wgpu buffer alignment (§7.3).

**Precision rule:** `f64` for all simulation arithmetic. `f32` for GUI trajectory frames only — the renderer does not need sub-angstrom precision, and halving the buffer width halves GPU upload cost.

---

#### 8.1 `SimState` — Mutable Integration State

The core mutable object passed through every step of the BAB integrator. Laid out as SoA so the inner force loop streams contiguous memory.

```rust
struct SimState {
    // Positions (Å) — SoA, f64
    pos_x: Vec<f64>,
    pos_y: Vec<f64>,
    pos_z: Vec<f64>,

    // Velocities (Å/ps) — SoA, f64
    vel_x: Vec<f64>,
    vel_y: Vec<f64>,
    vel_z: Vec<f64>,

    // Net forces (kcal/mol·Å) — SoA, f64; overwritten each step
    force_x: Vec<f64>,
    force_y: Vec<f64>,
    force_z: Vec<f64>,

    // Scalar energies for logging
    potential_energy: f64,
    kinetic_energy:   f64,

    step: u64,
}
```

All three position arrays are length $N$; same for velocities and forces. $N$ is fixed for the lifetime of a run.

---

#### 8.2 `AtomData` — Static Per-Atom Properties

Read-only after initialization. Not on the hot loop — AoS is acceptable, but SoA used for LJ parameters to allow vectorized combination-rule evaluation.

```rust
struct AtomData {
    // LJ parameters per atom — SoA, used in non-bonded inner loop
    epsilon: Vec<f64>,   // kcal/mol
    sigma:   Vec<f64>,   // Å
    mass:    Vec<f64>,   // amu
    charge:  Vec<f64>,   // elementary charge (reserved; not used in v1)

    // Metadata — AoS, used only for I/O and GUI coloring
    meta: Vec<AtomMeta>,
}

struct AtomMeta {
    element:    Element,      // enum: H, C, N, O, S, ...
    residue_id: u32,
    atom_name:  [u8; 4],      // PDB atom name field
    chain_id:   u8,
}
```

`Element` is a `#[repr(u8)]` enum. `AtomMeta` is 8 bytes — fits one cache line per 8 atoms.

---

#### 8.3 `BondedTopology` — Force Field Connectivity

SoA layout throughout: the bonded force loop iterates over all bonds in sequence, reading `i[n]`, `j[n]`, `k[n]`, `r0[n]` — contiguous access in each array.

```rust
struct BondedTopology {
    // Bond stretching — one entry per bond
    bond_i:  Vec<u32>,
    bond_j:  Vec<u32>,
    bond_k:  Vec<f64>,   // force constant (kcal/mol·Å²)
    bond_r0: Vec<f64>,   // equilibrium length (Å)

    // Angle bending — one entry per angle
    angle_i:   Vec<u32>,
    angle_j:   Vec<u32>,  // central atom
    angle_k:   Vec<u32>,
    angle_kth: Vec<f64>,  // force constant (kcal/mol·rad²)
    angle_th0: Vec<f64>,  // equilibrium angle (rad)

    // Dihedral torsion — one entry per dihedral
    dihed_i:     Vec<u32>,
    dihed_j:     Vec<u32>,
    dihed_k:     Vec<u32>,
    dihed_l:     Vec<u32>,
    dihed_kphi:  Vec<f64>,  // barrier height (kcal/mol)
    dihed_n:     Vec<u32>,  // multiplicity
    dihed_delta: Vec<f64>,  // phase (rad)

    // Holonomic constraints for Geodesic BAOAB A-step
    // (typically bond lengths involving hydrogen)
    constr_i:   Vec<u32>,
    constr_j:   Vec<u32>,
    constr_dsq: Vec<f64>,  // target |r_i - r_j|² (Å²)
}
```

---

#### 8.4 `NeighborList` — Verlet Pair List

Rebuilt when any atom displaces more than $(r_s - r_c)/2$ since the last build (checked via squared displacement to avoid a square root).

```rust
struct NeighborList {
    // Flat list of all pairs (i, j) with i < j within r_skin
    pair_i: Vec<u32>,
    pair_j: Vec<u32>,

    // Positions at last rebuild — used for displacement check
    ref_x: Vec<f64>,
    ref_y: Vec<f64>,
    ref_z: Vec<f64>,

    r_cutoff: f64,   // r_c — force goes to zero beyond this
    r_skin:   f64,   // r_s — list radius; r_s > r_c
    r_switch: f64,   // r_sw — switching function onset
}
```

`pair_i` and `pair_j` are SoA: the non-bonded inner loop loads `pair_i[n]` and `pair_j[n]` to gather positions, then scatters forces — contiguous index access improves prefetcher behavior.

---

#### 8.5 `ForceBuffer` — Per-Thread Accumulator

One buffer per thread, each of length $N$. Allocated once at startup, zeroed at the start of each step.

```rust
struct ForceBuffer {
    fx: Vec<f64>,   // length N
    fy: Vec<f64>,
    fz: Vec<f64>,
}

// Simulation engine holds:
thread_force_buffers: Vec<ForceBuffer>,  // length T
```

After all threads complete, the master thread reduces into `SimState::force_{x,y,z}` in thread-index order.

---

#### 8.6 `SimParams` — Immutable Run Configuration

Created once from the input file; never mutated during a run. Shared via `Arc<SimParams>` across threads.

```rust
struct SimParams {
    n_atoms:    usize,
    n_steps:    u64,
    dt:         f64,    // timestep (ps); typically 0.004 ps = 4 fs with Geodesic BAOAB
    box_size:   [f64; 3],  // simulation box (Å); cubic assumed in v1

    r_cutoff:   f64,
    r_skin:     f64,
    r_switch:   f64,

    max_constr_iter: u32,   // I_max for Geodesic A-step
    constr_tol:      f64,   // convergence threshold for |λ_i|

    frame_interval:  u32,   // steps between trajectory snapshots
    n_threads:       usize, // T; recorded in log for reproducibility

    total_energy:    f64,   // E; defines the Jacobi metric (§2.0)
}
```

`total_energy` is set from the initial state and held constant (NVE ensemble). It enters the Jacobi metric $g^J_{ij} = 2(E - V)\,m_i\,\delta_{ij}$ at every step.

---

#### 8.7 `TrajectoryFrame` — Snapshot for Output and GUI

Written to the ring buffer at every `frame_interval` steps. Positions are downcast to `f32` at write time — the simulation continues in `f64`.

```rust
struct TrajectoryFrame {
    step:   u64,
    time_ps: f64,

    // f32 positions for GPU upload to renderer — half the bandwidth of f64
    pos_x: Vec<f32>,
    pos_y: Vec<f32>,
    pos_z: Vec<f32>,

    // Per-atom PSL flexibility score (computed post-hoc; 0.0 until PSL runs)
    flexibility: Vec<f32>,

    potential_energy: f64,
    kinetic_energy:   f64,
}
```

The ring buffer is a `Vec<TrajectoryFrame>` of fixed capacity, indexed with a wrapping atomic counter. The simulation thread writes; the GUI thread and file-output thread read. No allocation after startup.

---

### 9. Crate & Workspace Layout

The repository is a **Cargo workspace**. Each crate has a single, bounded responsibility. The key invariant is that `geodesic-core` depends on nothing internal — all other crates depend on it, never on each other, except where explicitly listed below.

Heavy optional dependencies (wgpu, ndarray-linalg, egui, Ripser FFI) are isolated to their own crates and gated behind Cargo features so that an M1 CPU-only build pulls in none of them.

---

#### 9.1 Directory Structure

```
geodesic-m/
├── Cargo.toml              ← workspace root; no src/
├── Cargo.lock
├── docs/
│   └── SAD.md
├── geodesic-core/          ← shared types, traits, errors — no heavy deps
├── geodesic-engine/        ← force field + integrator + CPU backend
├── geodesic-gpu/           ← GPU backend (feature = "gpu")
├── geodesic-topo/          ← PSL + Zigzag pipeline (feature = "topo")
├── geodesic-io/            ← PDB/DCD/JSON file I/O
├── geodesic-gui/           ← wgpu renderer + export UI (feature = "gui")
└── geodesic/               ← binary: CLI, backend selection, orchestration
```

---

#### 9.2 Crate Responsibilities

| Crate | Type | Responsibility | Key deps |
| :--- | :--- | :--- | :--- |
| `geodesic-core` | lib | `SimState`, `AtomData`, `SimParams`, `BondedTopology`, `NeighborList`, `ForceBuffer`, `TrajectoryFrame`, `ComputeBackend` trait, all error types | none |
| `geodesic-engine` | lib | Force field evaluation (bonded + LJ), Geodesic BAOAB integrator, constraint solver, Verlet list, `CpuBackend` impl | `geodesic-core`, `rayon` |
| `geodesic-gpu` | lib | `GpuBackend` impl, wgpu device management, WGSL force compute shader | `geodesic-core`, `wgpu` |
| `geodesic-topo` | lib | Contact complex construction, PSL eigenvalue pipeline, Zigzag persistence (Ripser FFI) | `geodesic-core`, `ndarray`, `ndarray-linalg` |
| `geodesic-io` | lib | PDB/XYZ parser, DCD trajectory writer, JSON barcode serializer, CSV energy log | `geodesic-core` |
| `geodesic-gui` | lib | wgpu 3D atom/bond renderer, ring buffer consumer, export panel | `geodesic-core`, `wgpu`, `egui`, `winit` |
| `geodesic` | bin | Argument parsing, backend selection, run orchestration | all above |

**v0.1 CLI scope (M1 only):** Two subcommands ship in the first release:
- `geodesic energy <protein.prmtop> <protein.inpcrd>` — evaluates and prints `total_energy` (required to fill `config.toml` before a run)
- `geodesic run <config.toml>` — runs the simulation, writes `output.dcd` + `energy.csv`

All other subcommands (`snapshot`, `topo`, `analyze`) are added in later milestones alongside the features that need them.

---

#### 9.3 Dependency Graph

```
geodesic (bin)
├── geodesic-core
├── geodesic-engine  →  geodesic-core
├── geodesic-io      →  geodesic-core
├── geodesic-gpu     →  geodesic-core          [feature = "gpu"]
├── geodesic-topo    →  geodesic-core          [feature = "topo"]
└── geodesic-gui     →  geodesic-core          [feature = "gui"]
```

No crate other than the binary depends on more than one sibling. The graph is a strict DAG with `geodesic-core` as the root.

---

#### 9.4 Feature Flags

Defined on the `geodesic` binary crate:

| Feature | Enables | Adds deps |
| :--- | :--- | :--- |
| *(none — default)* | M1 CPU-only headless build | `rayon` only |
| `gpu` | `geodesic-gpu`, GPU backend selectable at runtime | `wgpu` |
| `topo` | `geodesic-topo`, PSL + Zigzag post-processing | `ndarray-linalg`, Ripser FFI |
| `gui` | `geodesic-gui`, 3D viewer + export UI | `wgpu`, `egui`, `winit` |

`gui` implies a wgpu dependency regardless of whether `gpu` is also enabled — the renderer always needs a graphics API. When both `gpu` and `gui` are enabled, the binary passes a single `Arc<wgpu::Device>` to both, sharing the GPU context.

---

#### 9.5 Module Map for `geodesic-engine` (M1 Scope)

The engine crate is the largest and most performance-sensitive. Its internal module structure:

```
geodesic-engine/src/
├── lib.rs
├── cpu_backend.rs      ← CpuBackend: ComputeBackend — entry point for all CPU dispatch
├── force/
│   ├── mod.rs
│   ├── nonbonded.rs    ← LJ inner loop (SoA, AVX2 hot path)
│   └── bonded.rs       ← bonds, angles, dihedrals
├── integrator.rs       ← Geodesic BAB outer loop
├── constraint.rs       ← iterative Lagrangian solver, ConvergenceError
└── neighbor.rs         ← Verlet list build + displacement check
```

`force/nonbonded.rs` is the only file that carries `#[target_feature]` annotations. All other modules are architecture-independent.

---

### 10. I/O Formats

All file I/O is handled exclusively by `geodesic-io`. The simulation engine and topology pipeline have no file system access — they operate on in-memory types from `geodesic-core`. This boundary is enforced by the crate graph (§9.3).

---

#### 10.1 Input Pipeline

A complete run requires three inputs:

```
config.toml          →  SimParams
protein.prmtop       →  AtomData + BondedTopology + LJ parameters
protein.inpcrd       →  initial SimState (positions + optional velocities)
```

**Why AMBER `.prmtop` + `.inpcrd` as the primary format:**
CHARMM and GROMACS require matching a topology file, a parameter file, and a coordinate file — three sources that must agree on atom naming conventions. AMBER's `.prmtop` is a single self-contained file encoding topology, atom types, LJ parameters, masses, charges, and connectivity. `.inpcrd` contains only positions (and optionally velocities). There is no cross-file atom-type matching to get wrong. This makes it the least error-prone input surface for v1.

A secondary **PDB-only input mode** is supported for test cases and validation runs where bonded topology is not needed (non-bonded-only or coarse-grained models). In this mode, `BondedTopology` is empty and the constraint list is empty — the integrator degenerates to standard Verlet.

| File | Format | Parser location | What it populates |
| :--- | :--- | :--- | :--- |
| `config.toml` | TOML | `geodesic-io::config` | `SimParams` |
| `protein.prmtop` | AMBER prmtop v1 | `geodesic-io::prmtop` | `AtomData`, `BondedTopology` |
| `protein.inpcrd` | AMBER inpcrd | `geodesic-io::inpcrd` | `SimState` (positions, velocities) |
| `protein.pdb` | PDB ATOM/HETATM | `geodesic-io::pdb` | `AtomData` (positions + metadata only) |

Future: CHARMM PSF/PRM and Gromacs GRO/TOP parsers added as optional `geodesic-io` features.

---

#### 10.2 `config.toml` Schema

All simulation parameters are set in a single TOML file. Unknown keys are rejected (no silent ignoring of typos).

```toml
[run]
n_steps         = 5_000_000
frame_interval  = 500           # write one frame every 500 steps
backend         = "cpu"         # "cpu" | "gpu" | "hybrid"
n_threads       = 0             # 0 = auto-detect physical cores

[system]
prmtop   = "protein.prmtop"
inpcrd   = "protein.inpcrd"
box_size = [80.0, 80.0, 80.0]  # Å; cubic only in v1

[integrator]
dt              = 0.004         # ps (4 fs — Geodesic BAOAB default)
total_energy    = -48230.5      # kcal/mol; defines the Jacobi metric

[nonbonded]
r_cutoff = 12.0                 # Å
r_skin   = 14.0                 # Å
r_switch = 10.0                 # Å

[constraints]
max_iter  = 100
tolerance = 1.0e-6              # convergence threshold for |λ_i|

[output]
trajectory = "output.dcd"
energy_log = "energy.csv"
```

`total_energy` is a required field — it defines the Jacobi metric (§2.0) and cannot be inferred without running an initial energy evaluation. For convenience, the binary provides a `geodesic energy protein.prmtop protein.inpcrd` subcommand that computes and prints it before the run.

---

#### 10.3 Trajectory Output — DCD

Trajectory frames are written to a **DCD** (CHARMM/NAMD binary trajectory) file. DCD is chosen over XTC (Gromacs compressed) because:
- It is a simple binary format with no external library dependency (XTC requires libxdr).
- It is natively supported by VMD, MDAnalysis, MDTraj, and OVITO — the standard visualisation and analysis stack.
- Frames are fixed-size and can be appended without rewriting the header.

Each DCD frame contains `f32` positions for all $N$ atoms — matching `TrajectoryFrame::pos_{x,y,z}` directly. No precision conversion is needed at write time.

**Write strategy:** The file-output thread consumes frames from the same ring buffer as the GUI (§7.5) and appends them to the DCD file. The simulation loop is never stalled by disk I/O.

**DCD header fields set at run start:**

| Field | Value |
| :--- | :--- |
| `NFILE` | 0 (updated on close) |
| `ISTART` | 0 |
| `NSAVC` | `frame_interval` |
| `DELTA` | `dt` (in AKMA time units: 1 AKMA = 48.88 fs) |
| `NATOM` | $N$ |

---

#### 10.4 Energy Log — CSV

Written every `frame_interval` steps, appended line-by-line. Never buffered in memory beyond one line — safe to tail during a run.

```
step,time_ps,potential_kcal,kinetic_kcal,total_kcal,temperature_K
0,0.000,-48230.512,1847.331,-46383.181,298.14
500,2.000,-48228.901,1849.102,-46379.799,298.43
...
```

Temperature is derived from kinetic energy: $T = 2 E_k / (3 N k_B)$.

---

#### 10.5 PDB Snapshot — On-Demand Export

A single-frame PDB written on request from the GUI export panel or via the `geodesic snapshot` subcommand. Positions are taken from the most recent `TrajectoryFrame` in the ring buffer, upcast back to `f64`, and written in standard ATOM record format.

B-factor column is populated with the PSL flexibility score if available, or `0.00` otherwise — this allows VMD and PyMOL to color by flexibility directly without post-processing.

---

#### 10.6 Zigzag Persistence Barcode — JSON

Written at the end of the topology pipeline run (§7.6). One file per analysis.

```json
{
  "metadata": {
    "n_atoms": 1247,
    "n_frames": 10000,
    "frame_interval": 500
  },
  "barcode": [
    { "dim": 0, "birth": 0.0,  "death": 142.5 },
    { "dim": 1, "birth": 18.3, "death": 97.1  },
    { "dim": 1, "birth": 44.0, "death": 44.0  }
  ]
}
```

`dim` is the homological dimension ($H_0$, $H_1$, $H_2$). `birth` and `death` are step indices (not time in ps) to keep the barcode independent of the timestep choice. Infinite bars (`death = ∞`) are encoded as `-1.0`. The format is intentionally minimal — downstream analysis in Python/Julia loads it with a single `json.load` call.

---

### 11. Crate Dependencies

Dependencies are pinned by major version in the workspace `Cargo.toml`. Patch versions are managed by `Cargo.lock` and never edited manually.

---

#### 11.1 Per-Crate Dependency Table

| Crate | Dependency | Version | Purpose |
| :--- | :--- | :--- | :--- |
| `geodesic-core` | *(none)* | — | Zero external deps by design |
| `geodesic-engine` | `rayon` | 1 | Static force decomposition + constraint par_iter |
| `geodesic-gpu` | `wgpu` | 22 | Compute shaders + device management |
| `geodesic-gpu` | `bytemuck` | 1 | Safe `&[T]` → `&[u8]` for GPU buffer uploads |
| `geodesic-gpu` | `pollster` | 0.3 | Minimal async executor for wgpu init (avoids tokio) |
| `geodesic-topo` | `ndarray` | 0.16 | Dense array operations for PSL matrix |
| `geodesic-topo` | `ndarray-linalg` | 0.16 | LAPACK eigenvalue solve for PSL spectrum |
| `geodesic-topo` | `sprs` | 0.11 | Sparse matrix for sheaf coboundary maps |
| `geodesic-topo` | `cc` | 1 | Build script to compile Ripser C++ source |
| `geodesic-io` | `serde` | 1 | Derive serialization for JSON barcode |
| `geodesic-io` | `serde_json` | 1 | JSON barcode writer |
| `geodesic-io` | `toml` | 0.8 | `config.toml` parser (strict: unknown keys rejected) |
| `geodesic-gui` | `wgpu` | 22 | 3D atom/bond renderer |
| `geodesic-gui` | `egui` | 0.29 | Immediate-mode UI panels and export controls |
| `geodesic-gui` | `egui-wgpu` | 0.29 | egui ↔ wgpu integration layer |
| `geodesic-gui` | `winit` | 0.30 | Window creation and OS event loop |
| `geodesic-gui` | `bytemuck` | 1 | Vertex buffer casting |
| `geodesic` (bin) | `clap` | 4 | CLI argument parsing, derive API |
| `geodesic` (bin) | `tracing` | 0.1 | Structured logging (spans per step, per constraint iter) |
| `geodesic` (bin) | `tracing-subscriber` | 0.3 | Log formatting and filtering via `RUST_LOG` |

No `tokio` or async runtime anywhere — the simulation is synchronous CPU work; async would add overhead with no benefit.

---

#### 11.2 Workspace `Cargo.toml` (skeleton)

```toml
[workspace]
members = [
    "geodesic-core",
    "geodesic-engine",
    "geodesic-gpu",
    "geodesic-topo",
    "geodesic-io",
    "geodesic-gui",
    "geodesic",
]
resolver = "2"

[workspace.dependencies]
geodesic-core   = { path = "geodesic-core" }
geodesic-engine = { path = "geodesic-engine" }
geodesic-gpu    = { path = "geodesic-gpu",  optional = true }
geodesic-topo   = { path = "geodesic-topo", optional = true }
geodesic-io     = { path = "geodesic-io" }
geodesic-gui    = { path = "geodesic-gui",  optional = true }

rayon             = "1"
wgpu              = "22"
bytemuck          = { version = "1", features = ["derive"] }
pollster          = "0.3"
ndarray           = "0.16"
ndarray-linalg    = { version = "0.16", features = ["openblas-static"] }
sprs              = "0.11"
cc                = "1"
serde             = { version = "1", features = ["derive"] }
serde_json        = "1"
toml              = { version = "0.8", features = ["parse"] }
egui              = "0.29"
egui-wgpu         = "0.29"
winit             = "0.30"
clap              = { version = "4", features = ["derive"] }
tracing           = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
```

Workspace-level dependency declarations mean version constraints are defined once — individual crate `Cargo.toml` files reference these with `{ workspace = true }`, no version duplication.

---

#### 11.3 BLAS / LAPACK Backend

`ndarray-linalg` requires a BLAS/LAPACK implementation. The workspace uses **OpenBLAS with static linking** (`openblas-static` feature) for the following reasons:

- Static linking removes the need for users to install a system BLAS — the binary is self-contained.
- OpenBLAS is the standard open-source BLAS for scientific computing on x86.
- Intel MKL is faster on Intel hardware but requires a separate license agreement and download; it is supported as a future opt-in build profile.

**Runtime constraint (enforced at startup):**
```rust
std::env::set_var("OPENBLAS_NUM_THREADS", "1");
```
Set unconditionally at process start in `geodesic/src/main.rs`, before any crate initializes. This pins the OpenBLAS thread pool to 1, ensuring the PSL eigenvalue solve is deterministic regardless of how Rayon has partitioned other work.

---

#### 11.4 Ripser FFI

Ripser has no published Rust crate. The `geodesic-topo` build script (`build.rs`) compiles the Ripser C++ source directly:

```
crates/geodesic-topo/
├── build.rs              ← cc::Build::new().file("ripser/ripser.cpp").compile("ripser")
├── ripser/
│   └── ripser.cpp        ← vendored Ripser source (MIT licence)
└── src/
    └── zigzag.rs         ← unsafe extern "C" declarations + safe Rust wrapper
```

Ripser is vendored (copied into the repository) rather than fetched at build time. This keeps builds reproducible without a network connection and avoids a git submodule dependency on an external repository.

---

### 12. Error Handling

**Strategy: `Result<T, E>` everywhere in library code. No panics on the hot path. Every error is actionable — it names the step, atom, or constraint that caused it.**

A NaN that goes undetected for even one step spreads to all atom positions via the integrator and makes the trajectory appear to diverge chaotically — hiding the root cause. A constraint that silently fails to converge produces a physically invalid trajectory with no audit trail. Both violate the Auditability NFR. The error handling strategy is designed to catch both as early as possible and stop the run with full diagnostic context.

---

#### 12.1 Error Type Hierarchy

All error types live in `geodesic-core::error` so every crate can reference them without circular deps.

```rust
// geodesic-core/src/error.rs

pub enum SimError {
    Io(IoError),
    Config(ConfigError),
    Numerical(NumericalError),
    Convergence(ConvergenceError),
    Backend(BackendError),
    Topology(TopologyError),
}

pub enum NumericalError {
    NanInForce  { step: u64, atom: usize, component: Axis },
    NanInPos    { step: u64, atom: usize, component: Axis },
    EnergyDrift { step: u64, drift_kcal: f64, threshold_kcal: f64 },
}

pub enum ConvergenceError {
    ConstraintSolver {
        step:           u64,
        constraint_idx: usize,
        atom_i:         usize,
        atom_j:         usize,
        residual:       f64,
        max_iter:       u32,
    },
}

pub enum ConfigError {
    UnknownKey(String),
    InvalidValue   { key: String, value: String, reason: String },
    MissingRequired(String),
    PhysicallyInvalid { description: String },  // e.g. r_switch >= r_cutoff
}

pub enum BackendError {
    DeviceLost,
    ShaderCompilation(String),
    OutOfGpuMemory,
}

pub enum TopologyError {
    EigensolverFailed { reason: String },
    RipserOutOfMemory,
}
```

`SimError` implements `std::error::Error`. All variants carry enough context to produce an actionable message without inspecting call-stack frames.

---

#### 12.2 No Panics in Library Code

`unwrap()` and `expect()` are forbidden in all library crates (`geodesic-core`, `geodesic-engine`, `geodesic-gpu`, `geodesic-topo`, `geodesic-io`, `geodesic-gui`). The only permitted exceptions are:

- **Tests** — `assert!`, `unwrap()` in `#[test]` functions.
- **Truly unreachable branches** — annotated with `unreachable!("reason")`, not `panic!()`, and only where the invariant is provable from the type system.

`geodesic/src/main.rs` (the binary) may `unwrap()` after startup validation is complete, since at that point all inputs have been checked.

---

#### 12.3 NaN Detection

NaN and infinite values in forces are checked **after every force evaluation, before positions are updated**. The check is a parallel scan over the three force SoA arrays — negligible cost compared to force calculation, but it catches divergence at the step it first occurs rather than several steps later when positions are already corrupted.

```
each step:
  1. compute_forces()          → writes SimState::force_{x,y,z}
  2. check_forces_finite()     → returns Err(NumericalError::NanInForce{..}) on first NaN
  3. geodesic_drift()          → updates positions
  4. check_positions_finite()  → returns Err(NumericalError::NanInPos{..}) on first NaN
  5. second B half-step
```

On detection, the run stops immediately. The DCD file is closed (header frame count written), the energy CSV is flushed, and the error is printed with atom index, residue name, and step number. The partially written trajectory up to the failing step is valid and readable.

**Diagnostic hint printed with NaN errors:**

```
ERROR geodesic::engine: NaN in force_x
  step=8421  atom=334  residue=GLY-22  chain=A
  hint: check for clashing atoms in the initial geometry
        (run `geodesic energy` to inspect initial forces before a long run)
```

---

#### 12.4 Energy Drift Monitoring

In the NVE ensemble, total energy $E = E_k + V$ is conserved. Drift indicates numerical instability (timestep too large, constraint solver tolerance too loose). Monitoring is configurable:

```toml
[monitoring]
energy_drift_threshold_kcal = 1.0   # warn if |E(t) - E(0)| exceeds this
energy_drift_action = "warn"        # "warn" | "stop"
```

Default: warn every `frame_interval` steps, never stop. Setting `energy_drift_action = "stop"` converts drift beyond threshold into a hard `NumericalError::EnergyDrift` that stops the run.

---

#### 12.5 Error Propagation in the Main Loop

The simulation loop in `geodesic-engine` is a single function returning `Result<(), SimError>`. Every fallible call uses `?`:

```rust
pub fn run(state: &mut SimState, backend: &mut dyn ComputeBackend,
           params: &SimParams) -> Result<(), SimError>
{
    for step in 0..params.n_steps {
        if neighbor_list_needs_rebuild(state, &params) {
            backend.build_neighbor_list(state, params)?;
        }
        backend.compute_forces(state)?;
        check_forces_finite(state, step)?;          // NaN guard
        backend.geodesic_drift(state, params.dt)?;  // ConvergenceError here
        check_positions_finite(state, step)?;
        check_energy_drift(state, params, step)?;   // EnergyDrift monitor
        emit_frame(state, step, params);            // ring buffer write, infallible
    }
    Ok(())
}
```

`main.rs` calls `run(...)` and on `Err` executes the shutdown sequence — flush CSV, close DCD, print error — then exits with code 1.

---

#### 12.6 GPU Backend Errors

wgpu operations fail asynchronously. `GpuBackend` maps wgpu errors to `BackendError` variants:

- `DeviceLost` → hard stop, no fallback to CPU. The user selected GPU explicitly; silently falling back would produce a different (CPU) result and violate the Determinism NFR.
- `ShaderCompilation` → hard stop at startup before any steps run. The WGSL shader is vendored and should never fail to compile — if it does, it is a bug, not a recoverable condition.
- `OutOfGpuMemory` → hard stop with a message suggesting reducing $N$ or switching to CPU.

---

#### 12.7 Topology Pipeline Errors

PSL and Zigzag run after the simulation completes. Errors here do not retroactively invalidate the trajectory:

- `EigensolverFailed` → logged as error, PSL output skipped, JSON barcode not written. The DCD trajectory is still intact.
- `RipserOutOfMemory` → logged as error, Zigzag barcode skipped. Suggest reducing the number of trajectory frames fed to Ripser via a frame-stride parameter.

---

### 13. Testing Strategy

**Core principle: force correctness and integrator correctness are tested independently via different methods. A passing energy conservation test does not prove force correctness, and a passing finite-difference gradient check does not prove the integrator is symplectic.**

Tests live in each crate. `cargo test` (debug build) runs the full suite in CI. Benchmarks run in release build on a dedicated CI runner and compare against a stored baseline — a regression of more than 10% fails the build.

---

#### 13.1 Test Fixtures

Small systems stored in `crates/geodesic-engine/tests/fixtures/` and `crates/geodesic-io/tests/fixtures/`. All fixture files are committed to the repository.

| Fixture | $N$ | Purpose |
| :--- | :--- | :--- |
| `lj_pair.inpcrd` / `.prmtop` | 2 | Analytical LJ force verification |
| `harmonic_dimer.inpcrd` / `.prmtop` | 2 | Bond force + constraint solver verification |
| `water_box_4.inpcrd` / `.prmtop` | 12 | PBC + angle forces + constraints (4 TIP3P waters) |
| `ala_dipeptide.inpcrd` / `.prmtop` | 22 | Full force field (bonds, angles, dihedrals, LJ) |
| `ala_dipeptide_ref.dcd` | 22 | Golden reference trajectory (100 steps) |

`ala_dipeptide` is the standard small-molecule MD benchmark — widely used in the literature for force field validation, with published reference energies.

---

#### 13.2 Force Correctness — Finite Difference Gradient Check

The canonical test for force implementation correctness. For each atom $i$ and each Cartesian component $\alpha$:

$$F_{i\alpha}^{\text{numeric}} = -\frac{V(\mathbf{r} + \varepsilon\,\hat{e}_{i\alpha}) - V(\mathbf{r} - \varepsilon\,\hat{e}_{i\alpha})}{2\varepsilon}$$

with $\varepsilon = 10^{-4}$ Å. The analytic force from `compute_forces()` is compared against $F^{\text{numeric}}$:

$$\frac{|F_{i\alpha}^{\text{analytic}} - F_{i\alpha}^{\text{numeric}}|}{|F_{i\alpha}^{\text{numeric}}| + 1} < \delta = 10^{-5}$$

This test runs on all four fixture systems and covers LJ, bond, angle, and dihedral forces independently. It catches sign errors, missing terms, wrong power-law exponents, and incorrect chain-rule applications. It does **not** test the integrator.

Location: `crates/geodesic-engine/tests/gradient_check.rs`

---

#### 13.3 Newton's Third Law

For every pair $(i, j)$ evaluated by the force engine, $F_{ij} = -F_{ji}$. Tested by summing all force contributions and asserting:

$$\left| \sum_{i=0}^{N-1} \mathbf{F}_i \right| < N \cdot \varepsilon_{\text{machine}}$$

A non-zero total force indicates an asymmetric accumulation bug — typically a missing Newton's-third-law shortcut in the pair loop. This is a property test that runs on all fixture systems.

Location: `crates/geodesic-engine/tests/newton_third_law.rs`

---

#### 13.4 Integrator Correctness — Energy Conservation

Tests that the Geodesic BAOAB BAB integrator is correctly implemented (not that forces are correct — that is §13.2). Runs on the harmonic dimer fixture with exact forces.

**Test:** Run 100,000 steps ($\Delta t = 0.004$ ps). Assert:

$$\frac{|E(t) - E(0)|}{|E(0)|} < 10^{-4}$$

A standard Velocity Verlet integrator on the same system would conserve energy to ~$10^{-5}$ — the tolerance is relaxed slightly to allow for the geodesic drift convergence tolerance. If this test fails it indicates a sequencing error in the BAB loop (e.g., forces evaluated at the wrong configuration) or a bug in the constraint solver that corrupts the Hamiltonian.

Location: `crates/geodesic-engine/tests/energy_conservation.rs`

---

#### 13.5 Determinism

Two independent runs with identical `SimParams` and initial `SimState` must produce bit-for-bit identical DCD output. Tested by running the `ala_dipeptide` fixture twice (same process, new `SimState` each time) and comparing the output byte arrays.

This test is **the only guard** against non-deterministic parallelism regressions. If a future change accidentally introduces an unordered Rayon reduction, this test catches it.

Location: `crates/geodesic-engine/tests/determinism.rs`

---

#### 13.6 Constraint Solver

Tested independently of the integrator, on the harmonic dimer:

- **Convergence:** Assert solver converges within $I_{\max}$ iterations for valid inputs.
- **Manifold adherence:** After the geodesic drift step, assert $\bigl| |r_i - r_j|^2 - d_0^2 \bigr| < \varepsilon_{\text{constr}}$ for every constrained pair.
- **Non-convergence error:** Set $I_{\max} = 1$ and assert `ConvergenceError` is returned, not a silent wrong result.

Location: `crates/geodesic-engine/tests/constraint_solver.rs`

---

#### 13.7 Golden Reference Trajectory

The stored `ala_dipeptide_ref.dcd` contains 100 frames produced by the first verified-correct build of GEODESIC-M. Subsequent builds must produce byte-identical output for the same input.

This test is **frozen on first passing build** and never regenerated unless the physics model changes (in which case it is regenerated deliberately with a commit message explaining why). Any change to force field parameters, integrator constants, or reduction order that silently changes the trajectory will break this test.

Location: `crates/geodesic-engine/tests/golden_reference.rs`

---

#### 13.8 I/O Round-Trip Tests

In `geodesic-io`:

| Test | Assertion |
| :--- | :--- |
| prmtop parser | Atom count, bond count, LJ parameters match known values for `ala_dipeptide.prmtop` |
| inpcrd parser | Positions round-trip to `f64` within machine epsilon |
| DCD writer | Written frame count matches `NFILE` header field on close |
| TOML config | Unknown keys rejected; missing required keys rejected; physical invalidity (`r_switch >= r_cutoff`) rejected |
| JSON barcode | Serialized barcode deserializes to identical struct; infinite bars encoded as `-1.0` |

---

#### 13.9 Benchmarks

Run with `cargo bench` (release build, `criterion` crate). Stored baselines committed to the repository. A regression of >10% on any benchmark fails CI.

| Benchmark | System | Measures |
| :--- | :--- | :--- |
| `bench_lj_inner_loop` | $N = 10{,}000$ | Non-bonded force throughput (atom-pairs/second) |
| `bench_neighbor_rebuild` | $N = 10{,}000$ | Verlet list build time |
| `bench_constraint_solver` | $N = 1{,}000$ constraints | Iterations/second for geodesic A-step |
| `bench_full_step` | `ala_dipeptide` | Wall time per simulation step end-to-end |

Benchmarks run on a pinned core count ($T = 1$ and $T = 8$) to detect parallelism regressions separately from single-core regressions.

---

#### 13.10 CI Matrix

```
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - cargo test                          # default features (CPU only)
      - cargo test --features topo          # topology pipeline
      - cargo clippy -- -D warnings
      - cargo fmt --check

  bench:
    runs-on: [self-hosted, bench]           # dedicated runner, pinned hardware
    steps:
      - cargo bench -- --baseline main      # compare against main branch baseline
```

GPU tests (`--features gpu`) are not in CI — they require physical GPU hardware. They are run manually before GPU milestone releases.
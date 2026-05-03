## Software Architecture Document: Deterministic Protein Manifold Simulator (DPMS)

### 1. System Overview
The **DPMS** is a high-precision simulation environment designed to explore protein conformational manifolds through classical mechanics. Unlike probabilistic LLM-based predictors, this system provides a deterministic, physics-based mapping of the high-dimensional energy landscape of a protein sequence.

---

### 2. Mathematical Foundation

The system operates in the classical mechanics regime of the **Karplus-Levitt-Warshel** framework. The total potential energy is partitioned into bonded and non-bonded contributions:

$$V_{\text{total}} = V_{\text{bonded}} + V_{\text{non-bonded}}$$

#### 2.1 Bonded Interactions

| Term | Equation | Parameters |
| :--- | :--- | :--- |
| Bond stretching | $V_b = k_b(r - r_0)^2$ | $k_b$: force constant, $r_0$: equilibrium length |
| Angle bending | $V_\theta = k_\theta(\theta - \theta_0)^2$ | $k_\theta$: force constant, $\theta_0$: equilibrium angle |
| Dihedral torsion | $V_\phi = k_\phi[1 + \cos(n\phi - \delta)]$ | $k_\phi$: barrier height, $n$: multiplicity, $\delta$: phase |

Parameters are sourced from a standard force field file (AMBER or CHARMM format — TBD in I/O section).

#### 2.2 Non-Bonded Interactions

**Lennard-Jones (van der Waals):**

$$V_{LJ}(r) = 4\epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^6 \right]$$

Truncated at cutoff $r_c$ with a smooth switching function over $[r_{\text{sw}}, r_c]$ to avoid discontinuous forces at the boundary.

#### 2.3 Equations of Motion

Newton’s second law governs each particle:

$$F_i = -\nabla_i V_{\text{total}} = m_i \ddot{r}_i$$

#### 2.4 Velocity Verlet Integration (2 fs timestep)

$$r_i(t + \Delta t) = r_i(t) + v_i(t)\,\Delta t + \tfrac{1}{2}a_i(t)\,\Delta t^2$$

$$a_i(t + \Delta t) = \frac{F_i(t + \Delta t)}{m_i}$$

$$v_i(t + \Delta t) = v_i(t) + \tfrac{1}{2}\bigl[a_i(t) + a_i(t + \Delta t)\bigr]\Delta t$$

Velocity Verlet is symplectic and time-reversible, preserving the Hamiltonian on long timescales. Forces at $t + \Delta t$ must be fully computed before velocities are updated — this is the key sequencing constraint for the integrator.

#### 2.5 Periodic Boundary Conditions (PBC)

For a simulation box of side $L$, the minimum image convention gives the effective pairwise displacement:

$$r_{ij}^* = r_{ij} - L \cdot \operatorname{round}\!\left(\frac{r_{ij}}{L}\right)$$

All pairwise distances use $r_{ij}^*$. Atoms that leave the box are wrapped back to $[0, L)$.

#### 2.6 Neighbor Lists

Naive pairwise evaluation is $O(N^2)$. A **Verlet neighbor list** reduces average cost to $O(N)$ amortized per step:

- A skin distance $r_s > r_c$ defines the list radius.
- The list is rebuilt when any atom has displaced more than $\tfrac{r_s - r_c}{2}$ since the last build.
- Only pairs within the list are evaluated; pairs with $r > r_c$ contribute zero force.

---

### 3. Component Architecture

| Component | Responsibility | Architectural Pattern |
| :--- | :--- | :--- |
| **Atomic Discretizer** | Converts protein coordinate data into a high-dimensional point cloud. | Data-Oriented Design |
| **Force Field Engine** | Computes bonded and non-bonded atomic interactions. | Parallel Task Processing |
| **ODE Integrator** | Executes time-stepping (approx. **2fs** increments) using Velocity Verlet algorithms. | Numerical Solver |
| **Manifold Projector** | Applies **TDA** or **PCA** to reduce trajectory data into interpretable 3D energy landscapes. | Dimensionality Reduction |

---

### 4. Data Flow
1.  **Input:** Protein sequence and initial spatial coordinates.
2.  **Processing:** Continuous loop of force summation followed by spatial coordinate updates (Integration).
3.  **Storage:** Trajectory snapshots representing movement across the manifold.
4.  **Output:** A topological map of the **Folding Funnel** and identified local energy minima.

---

### 5. Non-Functional Requirements (NFRs)

*   **Determinism:** Bit-for-bit identical output given the same binary, hardware, initial coordinates, and timestep. Floating-point operations must follow IEEE 754 strict ordering; no non-deterministic parallelism (e.g., unordered floating-point reductions).
*   **Auditability:** Every force contribution must be attributable to a specific atom pair and interaction type, queryable by step index. No implicit accumulation or fused operations that obscure intermediate values.
*   **Performance:** Near-linear thread scaling up to the available core count. Concrete throughput targets (ns/day per atom count) to be defined once simulation scale is decided (see open decision: §3 parallelism model).

---

### 6. Implementation Constraints
*   **Language:** **Rust** — mandatory for memory safety, zero-cost abstractions, and deterministic floating-point behavior without a garbage collector.
*   **Hardware:** Optimized for SIMD (Single Instruction, Multiple Data) or GPU-accelerated computing to handle $O(N^2)$ force calculations.
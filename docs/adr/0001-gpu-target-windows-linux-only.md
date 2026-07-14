# 0001 — GPU target is Windows + Linux only

Status: Accepted
Date: 2026-07-14
Deciders: Felipe Carvajal Brown

## Context

The v0.5 GPU backend (`geodesic-gpu`, M2) uses wgpu compute shaders. wgpu can target DX12, Vulkan, Metal, and GL. Supporting Metal means owning a macOS toolchain, macOS CI, and a third graphics path we cannot test on the developer's Windows machine. GEODESIC-M's development and CI target Windows and Linux; no macOS hardware is in scope.

## Decision

The GPU backend targets Windows and Linux only. wgpu is initialized with the DX12 and Vulkan backends enabled; Metal and GL are not requested. macOS is unsupported for `--features gpu`. The default CPU-only build remains cross-platform.

## Consequences

- `--features gpu` builds and runs on Windows (DX12, including the WARP software adapter) and Linux (Vulkan, including lavapipe). It will not produce a GPU adapter on macOS.
- GPU tests are adapter-adaptive: they skip with a logged reason when no DX12/Vulkan adapter is present, so a macOS or headless machine sees a clean skip, not a failure.
- Dependency version: SAD §14 pins `wgpu = "22"`. Resolved and built against **wgpu 22.1.0** (naga 22.1.0, wgpu-hal 22.0.0) on the current Rust toolchain (Windows); no bump was needed.

## Alternatives considered

- **Include Metal/macOS.** Rejected: no macOS hardware or CI in scope; a third untestable graphics path.
- **Single backend (Vulkan only).** Rejected: DX12 is the reliable default on Windows and gives WARP for CI without a Vulkan runtime.

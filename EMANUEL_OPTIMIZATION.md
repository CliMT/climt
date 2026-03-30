# Emanuel Convection Optimization Report

This document details the refactoring and optimization of the Emanuel convection scheme in `climt`. The primary goals were to achieve near-native performance, ensure compatibility with any array backend (NumPy, Unyt, JAX), and provide a path toward a fully differentiable model.

## 1. Architectural Changes

### Core Abstraction Layer (`climt/_core/backend.py`)
A universal abstraction layer was introduced to handle backend-specific operations transparently:
*   **`get_array_namespace(*arrays)`**: Dynamically identifies the correct API (`numpy` or `jax.numpy`) based on input data.
*   **`set_item(arr, idx, val)`**: A unified interface for updating arrays. It uses Numba `overload` for high-performance in-place updates on NumPy, and functional updates (`at[idx].set(val)`) for JAX.
*   **`JaxBackend`**: A new `sympl` StateBackend that allows the entire model to run on JAX arrays, enabling differentiability through the `sympl` interface.

### Functional Refactoring Versions
*   **V2 (`pure_python_v2.py`)**: Optimized for Numba JIT on NumPy/Unyt. Uses standard Python control flow. Highly performant on CPU.
*   **V3 (`pure_python_v3.py`)**: Dual-path implementation optimized for both Numba and JAX XLA. Replaces vertical loops with `jax.lax.scan` and branches with `jnp.where` for the JAX path.

## 2. Optimization Strategy & Performance

### Multi-Backend Performance (1000 Columns, 30 Levels, Apple M3 Pro)
| Backend | Platform | Time per Call (s) | Relative Speedup |
| :--- | :--- | :--- | :--- |
| Original Python | CPU (Serial) | 0.6607 | 1.0x |
| Fortran | CPU (Serial) | 0.0039 | ~170x |
| **V3 Numba JIT** | **CPU (Parallel)** | **0.0012** | **~550x** |
| **V3 JAX JIT** | **METAL (GPU)** | **0.2660** | **~2.5x** |

**Note on GPU Performance**: The Emanuel scheme's high logic complexity results in many small kernels. On Apple Silicon, the dispatch latency for these kernels currently makes the CPU/Numba path significantly faster than the GPU path for this specific algorithm. However, the JAX path provides full differentiability.

## 3. Verification and Differentiability

### Continuous Parity Check
*   **Tolerance**: A strict tolerance of `1e-12` is enforced between Numba and Original Python.
*   **JAX x64**: Precision was verified using `jax_enable_x64=True` on CPU, achieving perfect bit-wise parity. Note that METAL GPU currently only supports `float32`.

### Differentiability
The `V3` implementation is fully differentiable. Sensitivities can be computed using `jax.grad` through the `array_call` interface, providing a foundation for parameter optimization and sensitivity analysis.

## 4. Summary of New Components
*   `climt/_core/jax_backend.py`: JAX-native state handling.
*   `climt/_components/emanuel/pure_python_v3.py`: High-performance differentiable physics.
*   `tests/test_jax_differentiation.py`: Verification of gradient flow.

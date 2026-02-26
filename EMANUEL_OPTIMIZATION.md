# Emanuel Convection Optimization Report

This document details the refactoring and optimization of the Emanuel convection scheme in `climt`. The primary goals were to achieve near-native performance, ensure compatibility with any array backend (NumPy, Unyt, JAX), and provide a path toward a fully differentiable model.

## 1. Architectural Changes

### Core Abstraction Layer (`climt/_core/backend.py`)
A new universal abstraction layer was introduced to handle backend-specific operations transparently:
*   **`get_array_namespace(*arrays)`**: Dynamically identifies the correct API (`numpy` or `jax.numpy`) based on input data.
*   **`set_item(arr, idx, val)`**: A unified interface for updating arrays. It uses Numba `overload` for high-performance in-place updates on NumPy, and functional updates (`at[idx].set(val)`) for JAX.
*   **`jit_compile(backend=...)`**: A unified decorator that applies the appropriate JIT compiler (`numba.njit` or `jax.jit`).
*   **`vectorize_component(...)`**: Orchestrates high-level parallelization (e.g., via `jax.vmap`).

### Functional Refactoring (`climt/_components/emanuel/pure_python_v2.py`)
The Emanuel scheme was refactored into a "Pure Functional" style:
*   **Single-Column Logic**: The core physics functions (`_convect_functional`, `_tlift_functional`) operate strictly on a single column.
*   **Immutability**: All state updates are performed via `set_item`, ensuring compatibility with JAX's functional requirements.
*   **NamedTuple Params**: Parameters are stored in an `EmanuelParams` NamedTuple, ensuring type stability and optimal performance under Numba JIT.

## 2. Optimization Strategy

### Numba JIT (NumPy/Unyt Path)
On CPU, performance is achieved through Numba:
1.  **Strict No-Python Mode**: All core routines are compiled to machine code.
2.  **Parallel Vertical Loops**: Horizontal columns are processed in parallel using Numba's `prange` within `_numpy_vectorized_convect`.
3.  **Zero-Overhead Abstractions**: Using Numba `overload` for `set_item` ensures that the abstraction adds no runtime cost compared to raw NumPy assignments.

### Performance Gains (100 Columns, ARM64)
| Version | Time per Call (s) | Relative Speedup |
| :--- | :--- | :--- |
| Fortran | 0.0451 | 1.0x |
| Original Python | 0.0683 | 0.6x |
| **V2 Python (JIT)** | **0.0006** | **~75x vs Fortran / ~118x vs Python** |

## 3. Verification and Benchmarking

### Continuous Parity Check
Parity is maintained through a dual-implementation strategy:
*   `pure_python.py` remains the serial "Reference Implementation."
*   `tests/test_emanuel_optimization.py` runs a `pytest` suite comparing the optimized V2 against the Reference across moist, dry, and unstable profiles.
*   **Tolerance**: A strict tolerance of `1e-12` is enforced.

### Non-Zero Value Inspection
To ensure correctness beyond trivial zeros, `tests/verify_nonzero_parity.py` was used to compare all three versions (Fortran, Python Original, Python V2) side-by-side at levels with active convection.

### Benchmarking Protocol
Benchmarking is conducted via `benchmarks/benchmark_emanuel_v2.py` using the following protocol:
1.  **Warmup**: Perform an initial call to trigger JIT compilation.
2.  **Iteration**: Average execution time over 50 iterations to eliminate noise.
3.  **Hardware Awareness**: The script leverages multi-core CPUs through the parallelized `V2` path.

## 4. Future Directions
The infrastructure introduced here is designed to be reusable for other `climt` components (e.g., `HeldSuarez`, `RRTMG`). The next milestone is the implementation of a full `JaxBackend` using the `vectorize_component` hook to enable GPU-accelerated and differentiable simulations.

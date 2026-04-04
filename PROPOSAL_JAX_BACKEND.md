# Proposal: JAX Backend for Climt/Sympl (CPU-Optimized)

## Overview
This document outlines the design for a `JaxBackend` in `climt`. The primary goal is to leverage JAX's XLA compilation on CPU to fuse the many small arithmetic operations (unit conversions, dimension alignments, state updates) that occur in a typical `sympl` timestep.

## Core Components

### 1. JaxStateContainer
Similar to `UnytStateContainer`, this object wraps a `jax.numpy.ndarray` with metadata.
- **Attributes**: `data` (jnp.ndarray), `dims` (tuple), `units` (string/object).
- **Immutability**: Every operation (`__add__`, `__mul__`, etc.) must return a *new* `JaxStateContainer`.
- **Numpy Integration**: Implement `__array__` to return `np.asarray(self.data)`. On CPU, this should be a zero-copy operation, allowing seamless passing to Cython/Fortran wrappers.
- **Array Priority**: Set `__array_priority__ = 100.0` to ensure JAX-aware operations are preferred over standard NumPy.

### 2. JaxBackend (StateBackend Implementation)
- **`get_array`**: Should return a `jnp.ndarray`. On CPU, it can return a NumPy view if the target component is known to be a non-JAX wrapper.
- **`create_quantity`**: Converts input data to JAX arrays.
- **Device Locking**: Force CPU usage at initialization to avoid unnecessary GPU/TPU discovery overhead:
  ```python
  import jax
  jax.config.update('jax_platform_name', 'cpu')
  ```

### 3. JaxTimeDelta
A subclass of `timedelta` (similar to `UnytTimeDelta`).
- **`total_seconds()`**: Returns a JAX scalar (`jnp.float32/64`).
- **Purpose**: Allows the `AdamsBashforth` or `Leapfrog` steppers to remain inside a JAX tracer during a `jit` call.

## The "Metadata at the Edge" Strategy
To keep JAX computations fast and JIT-compatible:
1. **Unwrap on Entry**: `JaxBackend.get_array` extracts the raw numerical array.
2. **Compute**: The component performs pure JAX math on "naked" arrays.
3. **Wrap on Exit**: `JaxBackend.create_quantity` re-attaches units and dimensions based on the component's `properties` metadata.
4. **Unit Conversion**: Perform conversion using a lightweight library or internal logic *before* entering the JIT-compiled block.

## Handling Immutability in Sympl
`sympl` often assumes a mutable state dictionary (`state.update()`). 
- **Recommendation**: Create a `JaxState` dictionary wrapper.
- **In-place blocking**: Intercept `__setitem__` with slice notation (e.g., `state['temp'][:] = x`) and raise a `TypeError` explaining that JAX arrays are immutable.
- **Functional Updates**: Encourage the use of the `TendencyStepper` pattern where new state objects are returned rather than modified in place.

## Performance Opportunities (CPU Only)
1. **XLA Kernel Fusion**: Fusing `tendency * dt + state` into a single CPU kernel. This reduces cache misses and eliminates the overhead of multiple Python-to-C calls in NumPy. **Note**: While JAX uses XLA for all operations, kernel fusion and the resulting speedup only occur when code is wrapped in `@jax.jit`.
2. **Zero-Copy Interop**: On CPU, JAX and NumPy can often share the same memory buffer. This allows `climt` to use JAX for the "glue" logic (steppers, unit handling) while still calling heavy Fortran physics (like RRTMG) without significant data movement penalties.
3. **Control Flow**: JAX's `jax.lax.scan` can be used to JIT-compile the entire multi-thousand-iteration simulation loop, effectively removing the Python interpreter from the main execution path.

## JIT Strategy & Component Analysis

To maximize performance, the implementation will target different components based on their computational nature:

| Component Category | JIT Strategy | Expected Benefit |
| :--- | :--- | :--- |
| **"Glue" Logic** (Steppers, Unit/Dim handling) | **Full JIT** | **High**: Fuses many small "death by a thousand cuts" operations into a single kernel. |
| **Pure Python Physics** (`HeldSuarez`, `SimplePhysics`) | **Full JIT** | **High**: Eliminates temporary array allocations and fuses arithmetic. |
| **Fortran/C Wrappers** (`RRTMG`, `Emanuel`) | **Metadata at Edge** | **None**: JAX cannot JIT through compiled binaries; these must remain outside JIT blocks. |
| **Simulation Loop** | **`jax.lax.scan`** | **Critical**: The ultimate goal to achieve near-native execution speed for the entire model. |

## Why JAX over Numba?
While Numba is excellent for optimizing explicit `for` loops in imperative code, JAX is better suited for `climt/sympl` for several reasons:
- **Framework Overhead**: JAX can optimize the "graph" of operations across multiple component calls, which is where `sympl` spends much of its time. Numba struggles to "see" across the framework's object-oriented boundaries.
- **Composition**: JAX's functional paradigm matches the "State in, State out" philosophy of `sympl`.
- **Differentiability**: JAX-native components automatically become differentiable, opening doors for sensitivity analysis and parameter optimization in the future.
- **Whole-Model JIT**: `jax.lax.scan` provides a more robust path to compiling the entire simulation loop compared to Numba's `jitclass`.

## Implementation Roadmap
1. **Skeleton**: `climt/_core/jax_backend.py`.
2. **Operator Support**: Ensure `JaxStateContainer` supports all math required by `sympl` steppers.
3. **Unit Helper**: Implement a non-wrapping unit conversion utility.
4. **Validation**: Run `benchmarks/benchmark_gmd.py` with the JAX backend to compare CPU performance against `Unyt` and `DataArray`.

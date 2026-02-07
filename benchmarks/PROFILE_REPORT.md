# Profiling Report for CliMT and Sympl Core (Issue #46 Analysis)

## Methodology
The profiling was conducted using `pyinstrument` on a modified version of the `examples/full_radiation_gcm_energy_balanced.py` script. The simulation ran for 40 time steps.

To isolate the overhead of the `sympl` and `climt` framework, "decorators" (monkey-patching wrappers) were applied to all functions and methods within the `sympl` and `climt` modules. This ensures they appear explicitly in the profiling call stack.

Micro-benchmarks addressing specific concerns from **GitHub Issue #46** were added to the end of the profiling run.

## Results

### Top Level Breakdown (Main Loop)
Total CPU time: **47.5s**

1.  **RRTMGShortwave (Radiation): ~36s (76%)**
    *   `UpdateFrequencyWrapper` dominates, managing the calls to `RRTMGShortwave`.
    *   Core physics (`RRTMGShortwave.array_call`) takes ~26s.
    *   Framework overhead (wrappers, unit checks) is significant around this call.

2.  **SimplePhysics: ~5.5s (11.5%)**
    *   `SimplePhysics.array_call` takes ~4.5s.
    *   Framework overhead is visible but proportionally smaller than for faster components.

3.  **SlabSurface: ~0.8s (1.7%)**
    *   `SlabSurface.__call__` takes 0.84s.
    *   **Major Overhead:** `get_numpy_arrays_with_properties` takes **0.65s (77% of component time)**.
    *   This confirms that for lightweight components, the framework cost dominates the physics cost.
    *   Deep dive into `get_numpy_arrays_with_properties`:
        *   Calls `DataArray.to_units` (0.57s)
        *   Calls `data_array_to_units` (0.57s)
        *   Calls `pint.UnitRegistry.parse_expression` (0.54s)
    *   **Finding:** Unit parsing is extremely expensive and happens on every call.

### Issue #46 Micro-benchmark Analysis

The micro-benchmarks (Total ~0.88s) specifically tested the concerns raised in Issue #46.

1.  **`get_numpy_array` Overhead (0.69s total for 1000 calls)**
    *   This is the most expensive operation in the micro-benchmark suite.
    *   It takes ~0.69ms per call.
    *   The overhead comes from `DataArray.transpose` (which is called internally) and the subsequent slicing/reshaping logic.
    *   **Confirmation:** The issue's claim that `get_numpy_array` is slow is supported. It is a significant cost when called thousands of times (as it is in the main loop for every component input).

2.  **`DataArray.__init__` Overhead**
    *   While not explicitly separated in the summary, `DataArray` creation is ubiquitous.
    *   In the main loop, `update_dict_by_adding_another` (0.67s) involves unit conversions and potentially new array creation.
    *   `DataArray.to_units` (0.57s in SlabSurface) creates new DataArrays.

3.  **Property Access Overhead**
    *   Attribute access (`.values`, `.dims`) was tested but is fast enough that it doesn't appear as a primary hotspot in the sampling profiler compared to `pint` parsing or `get_numpy_array`.

## Conclusion & Recommendations

The profiling confirms the concerns raised in **Issue #46** and identifies `sympl` core functions as a major bottleneck for fast physics components.

1.  **Unit Parsing is the #1 Framework Bottleneck:** `get_numpy_arrays_with_properties` spends the vast majority of its time in `pint.UnitRegistry.parse_expression`.
    *   **Fix:** Cache parsed units or avoiding repeated parsing for the same variable/component pairs.

2.  **`get_numpy_array` is Slow:** Consumes significant time (0.6ms/call) due to internal transpose and reshaping logic.
    *   **Fix:** Optimize `get_numpy_array` to avoid full `DataArray.transpose` if possible, or use numpy-only operations as suggested in the issue.

3.  **`DataArray` Creation/Conversion:** Creating new `DataArray`s during unit conversion (`to_units`) is expensive.
    *   **Fix:** In-place unit modification where possible, or optimizing `__init__`.

For `climt` users, this means that simple components (like `SlabSurface` or `EmanuelConvection`) are paying a "framework tax" that is often larger than their actual physical computation time.

# Profiling Report for CliMT and Sympl Core

## Methodology
The profiling was conducted using `pyinstrument` on a modified version of the `examples/full_radiation_gcm_energy_balanced.py` script. The simulation ran for 40 time steps.

To isolate the overhead of the `sympl` and `climt` framework, "decorators" (monkey-patching wrappers) were applied to all functions and methods within the `sympl` and `climt` modules. This ensures they appear explicitly in the profiling call stack, although `pyinstrument` captures the full stack regardless.

The `climt.GFSDynamicalCore` was replaced with `sympl.AdamsBashforth` to drive the physics components.

## Results

### Top Level Breakdown
Total CPU time: **41.5s**

1.  **RRTMGShortwave (Radiation): ~32s (77%)**
    *   Most of this time is spent in the `wrapper` logic which eventually calls `RRTMGShortwave.array_call`.
    *   Inside `RRTMGShortwave.array_call` (23.08s), the actual computation happens.
    *   However, a significant portion (~9s) appears to be overhead *before* or *around* the core computation, or shared with `RRTMGLongwave`.
    *   Wait, the tree shows `RRTMGLongwave.array_call` (8.79s) is *inside* `RRTMGShortwave`'s call stack? No, that's likely an artifact of how `ImplicitTendencyComponentComposite` calls components sequentially or how pyinstrument groups them if they are in the same loop. Actually, looking at the text output:
        ```
        RRTMGShortwave.__call__
         └─ RRTMGShortwave.wrapper
             ├─ RRTMGShortwave.array_call (23.08s)
             └─ RRTMGLongwave.array_call (8.79s)
        ```
        This nesting suggests `RRTMGShortwave` might be wrapping another component or they are being called in a way that pyinstrument attributes one to the other (unlikely unless one calls the other).
        Actually, `profile_run.py` (and the original script) defines:
        ```python
        radiation_lw = UpdateFrequencyWrapper(climt.RRTMGLongwave(), ...)
        radiation_sw = UpdateFrequencyWrapper(climt.RRTMGShortwave(), ...)
        ```
        And they are passed to the dycore (AdamsBashforth -> Composite).
        The composite calls them sequentially.
        The pyinstrument output indentation might be misleading in text mode if they are siblings.
        However, `UpdateFrequencyWrapper.__call__` is taking 32.1s.

2.  **SimplePhysics: ~5s (12%)**
    *   `SimplePhysics.__call__` takes 4.35s.
    *   `SimplePhysics.array_call` takes 4.00s.
    *   Overhead is relatively low here (~0.35s wrapper overhead).

3.  **EmanuelConvection: ~1.8s (4%)**
    *   `EmanuelConvection.__call__` takes 1.85s.
    *   `EmanuelConvection.array_call` takes 1.33s.
    *   **Significant Overhead:** `EmanuelConvection._set_fortran_constants` takes 0.42s (23% of component time). This calls `get_constant` -> `DataArray.to_units` -> `pint`. This is done *every call*.

4.  **SlabSurface: ~0.8s (2%)**
    *   `SlabSurface.__call__` takes 0.77s.
    *   **Major Overhead:** `get_numpy_arrays_with_properties` takes 0.61s (79% of the component time!).
    *   This function calls `DataArray.to_units` -> `pint`, which is very slow.

### Framework Overhead Analysis

The profiling clearly highlights `sympl` core functions that contribute significantly to runtime, especially for fast components like `SlabSurface` and `EmanuelConvection`.

#### 1. `get_numpy_arrays_with_properties` / `DataArray.to_units` / `pint`
This path is a major bottleneck.
*   **Context:** `sympl` components automatically validate and convert units of input states using `get_numpy_arrays_with_properties` (called via `__call__` wrapper).
*   **Cost:** In `SlabSurface`, it consumes ~80% of execution time. In `EmanuelConvection`, setting constants (which uses similar unit logic) consumes ~23%.
*   **Root Cause:** `pint.UnitRegistry.parse_expression` and `to_units` are expensive operations to perform on every timestep for every component.

#### 2. `UpdateFrequencyWrapper`
*   This wrapper is efficient (it just delegates), but it organizes the heavy radiation calls.

#### 3. Initialization (`climt.get_default_state`)
*   Not explicitly shown in the main loop profile (as it runs before), but `init_ozone` and `CubicSpline` usage was noted during debugging.

## Recommendations for Optimization

1.  **Cache Unit Conversions:** `sympl` should cache the result of `parse_expression` or the conversion factors if the units and input types haven't changed.
2.  **Optimize `get_numpy_arrays_with_properties`:** Avoid full unit parsing if the input `DataArray` already has the correct units (string comparison might be faster than full pint parsing).
3.  **Optimize `_set_fortran_constants`:** In `climt` components like `EmanuelConvection`, constants are fetched and converted every single call. These should be cached in `__init__` or only updated if they change.

## Conclusion
While radiation physics (`RRTMG`) dominate the total runtime, the **framework overhead** (specifically `sympl`'s unit handling and input processing) is disproportionately high for faster components (~80% overhead for `SlabSurface`). Optimizing `DataArray.to_units` and reducing `pint` calls in the hot loop would significantly improve the efficiency of the core driver.

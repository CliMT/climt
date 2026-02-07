# Profiling Report for CliMT Example

## Methodology
The example script `examples/full_radiation_gcm_energy_balanced.py` was used as a base. However, the script required modifications to run in the current environment:

1.  **Missing `GFSDynamicalCore`**: The component `climt.GFSDynamicalCore` referenced in the example is missing from the current `climt` codebase. It was replaced with `sympl.AdamsBashforth` to step the physics components forward in time (without dynamics).
2.  **Bug Fix in `RRTMGShortwave`**: A `TypeError` was encountered in `climt/_components/rrtmg/sw/component.py` where `solar_cycle_fraction` (a 0-d numpy array) was passed to a Cython function expecting a double. This was patched by using `.item()` to extract the scalar value.
3.  **Reduced Iterations**: The simulation loop was reduced from 216,000 steps to 20 steps for feasible profiling.
4.  **Removed Plotting**: `PlotFunctionMonitor` was removed to avoid GUI dependencies.

## Results
Profiling was performed using `cProfile` for 20 time steps. The radiation components (`RRTMGShortwave` and `RRTMGLongwave`) were configured to update every 6 steps (4 calls total in 20 steps).

| Component | Cumulative Time (s) | Number of Calls | Time Per Call (s) |
| :--- | :--- | :--- | :--- |
| **RRTMGShortwave** | **14.80** | **4** | **3.70** |
| RRTMGLongwave | 5.59 | 4 | 1.40 |
| SimplePhysics | 2.40 | 20 | 0.12 |

## Conclusion
**`climt.RRTMGShortwave` is the slowest part of the code.** It consumes the majority of the execution time, taking approximately 3.7 seconds per call, which is significantly more than `RRTMGLongwave` (1.4s) and `SimplePhysics` (0.12s).

The primary bottleneck is in the `array_call` method of `RRTMGShortwave`, specifically the call to the underlying Fortran/Cython extension `_rrtmg_sw.rrtm_calculate_shortwave_fluxes`.

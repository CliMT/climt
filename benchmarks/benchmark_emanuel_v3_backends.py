import time
from datetime import timedelta

import numpy as np
import sympl
from sympl._core.backend import DataArrayBackend

import climt
from climt import (
    EmanuelConvection,
    EmanuelConvectionPython,
    RRTMGLongwave,
    RRTMGShortwave,
    SimplePhysics,
    SlabSurface,
    UnytBackend,
    get_default_state,
)


def run_benchmark(convection_cls, backend_obj, iterations=1000):
    """
    Runs a radiative-convective equilibrium simulation for a given backend and convection scheme.
    Based on examples/gmd_radiative_convective_python_emanuelv3.py
    """
    sympl.set_backend(backend_obj)

    timestep = timedelta(minutes=5)
    if isinstance(backend_obj, UnytBackend):
        from climt import UnytTimeDelta

        timestep = UnytTimeDelta(timestep)

    # Initialize components
    convection = convection_cls()
    radiation_sw = RRTMGShortwave()
    radiation_lw = RRTMGLongwave()
    slab = SlabSurface()
    simple_physics = SimplePhysics()

    # Get default state
    state = get_default_state(
        [simple_physics, convection, radiation_lw, radiation_sw, slab]
    )

    # Helper to set values across backends
    def set_values(obj, val):
        if hasattr(obj, "values"):
            obj.values[:] = val
        elif hasattr(obj, "data"):
            obj.data[:] = val
        else:
            obj[:] = val

    # Set initial conditions as in the example
    set_values(state["surface_albedo_for_direct_shortwave"], 0.5)
    set_values(state["surface_albedo_for_direct_near_infrared"], 0.5)
    set_values(state["surface_albedo_for_diffuse_shortwave"], 0.5)
    set_values(state["zenith_angle"], np.pi / 2.5)
    set_values(state["surface_temperature"], 300.0)
    set_values(state["ocean_mixed_layer_thickness"], 5)
    set_values(state["area_type"], "sea")

    time_stepper = sympl.AdamsBashforth([convection, radiation_lw, radiation_sw, slab])

    # --- Warmup Phase ---
    # We run a full iteration of the loop logic to trigger any JIT compilation (Numba/RRTMG)
    convection.current_time_step = timestep
    diagnostics, state = time_stepper(state, timestep)
    state.update(diagnostics)
    diagnostics, new_state = simple_physics(state, timestep)
    state.update(diagnostics)
    state.update(new_state)
    state["time"] += timestep
    set_values(state["eastward_wind"], 3.0)

    # --- Timed Phase ---
    start_time = time.perf_counter()

    for i in range(iterations):
        # Update convection timestep (required by Emanuel scheme in climt)
        convection.current_time_step = timestep

        # Step the main components
        diagnostics, state = time_stepper(state, timestep)
        state.update(diagnostics)

        # Step simple physics
        diagnostics, new_state = simple_physics(state, timestep)
        state.update(diagnostics)
        state.update(new_state)

        # Update time and boundary conditions
        state["time"] += timestep
        set_values(state["eastward_wind"], 3.0)

    end_time = time.perf_counter()
    return end_time - start_time


if __name__ == "__main__":
    iterations = 1000
    print(f"Benchmarking Emanuel Convection backends ({iterations} iterations)...")
    print("=" * 80)
    print(f"{'Configuration':<45} | {'Duration (s)':<15} | {'ms/step':<10}")
    print("-" * 80)

    results = {}

    configs = [
        ("Fortran + DataArray", EmanuelConvection, DataArrayBackend()),
        ("Fortran + Unyt", EmanuelConvection, UnytBackend()),
        ("V3 Python + DataArray", EmanuelConvectionPython, DataArrayBackend()),
        ("V3 Python + Unyt", EmanuelConvectionPython, UnytBackend()),
    ]

    for label, conv_cls, backend in configs:
        try:
            duration = run_benchmark(conv_cls, backend, iterations=iterations)
            ms_per_step = (duration / iterations) * 1000
            print(f"{label:<45} | {duration:<15.4f} | {ms_per_step:<10.2f}")
            results[label] = duration
        except Exception as e:
            print(f"{label:<45} | {'FAILED':<15} | N/A")
            print(f"Error: {e}")

    print("=" * 80)
    if "Fortran + DataArray" in results and "Fortran + Unyt" in results:
        f_da = results["Fortran + DataArray"]
        f_unyt = results["Fortran + Unyt"]
        print(f"Fortran Backend Speedup (Unyt/DA):   {f_da / f_unyt:.2f}x")

    if "V3 Python + DataArray" in results and "V3 Python + Unyt" in results:
        v3_da = results["V3 Python + DataArray"]
        v3_unyt = results["V3 Python + Unyt"]
        print(f"V3 Python Backend Speedup (Unyt/DA): {v3_da / v3_unyt:.2f}x")

    if "Fortran + Unyt" in results and "V3 Python + Unyt" in results:
        fortran_unyt = results["Fortran + Unyt"]
        v3_unyt = results["V3 Python + Unyt"]
        print(f"V3 (Unyt) vs Fortran (Unyt):         {fortran_unyt / v3_unyt:.2f}x")
    print("=" * 80)

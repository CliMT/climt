import time
from datetime import timedelta

import numpy as np
import sympl
from sympl._core.backend import DataArrayBackend

import climt
from climt import (
    EmanuelConvection,
    RRTMGLongwave,
    RRTMGShortwave,
    SimplePhysics,
    SlabSurface,
    UnytBackend,
    get_default_state,
)


def set_values(obj, val):
    if hasattr(obj, "values"):  # DataArray
        obj.values[:] = val
    elif hasattr(obj, "data"):  # UnytStateContainer
        # Ensure we don't overwrite with unitless if units expected, or unyt handles it?
        # unyt_array[:] = val works if val is compatible.
        # If val is unitless float, unyt assumes current units?
        obj.data[:] = val


def get_values_as_numpy(obj):
    if hasattr(obj, "values"):  # DataArray
        return np.array(obj.values)
    elif hasattr(obj, "data"):  # UnytStateContainer
        return np.array(obj.data)  # Strips units
    return np.array(obj)


def run_gmd_simulation(backend_obj, iterations=1000):
    print(f"Initializing for {type(backend_obj).__name__}...")
    sympl.set_backend(backend_obj)

    timestep = timedelta(minutes=5)
    if isinstance(backend_obj, UnytBackend):
        from climt import UnytTimeDelta

        timestep = UnytTimeDelta(timestep)

    convection = EmanuelConvection()
    radiation_sw = RRTMGShortwave()
    radiation_lw = RRTMGLongwave()
    slab = SlabSurface()
    simple_physics = SimplePhysics()

    # Note: get_default_state must be called AFTER set_backend to get correct state objects
    state = get_default_state(
        [simple_physics, convection, radiation_lw, radiation_sw, slab]
    )

    # Set initial conditions
    set_values(state["air_temperature"], 270)
    set_values(state["surface_albedo_for_direct_shortwave"], 0.5)
    set_values(state["surface_albedo_for_direct_near_infrared"], 0.5)
    set_values(state["surface_albedo_for_diffuse_shortwave"], 0.5)
    set_values(state["zenith_angle"], np.pi / 2.5)
    set_values(state["surface_temperature"], 300.0)
    set_values(state["ocean_mixed_layer_thickness"], 5)
    set_values(
        state["area_type"], "sea"
    )  # 'sea' string works for both? unyt backend stores string array as object/numpy

    time_stepper = sympl.AdamsBashforth([convection, radiation_lw, radiation_sw, slab])

    output_data = []

    start_time = time.perf_counter()

    for i in range(iterations):
        convection.current_time_step = timestep

        # Stepper call
        diagnostics, state = time_stepper(state, timestep)
        state.update(diagnostics)

        # SimplePhysics (Stepper-like) call
        diagnostics, new_state = simple_physics(state, timestep)
        state.update(diagnostics)
        state.update(new_state)

        if (i + 1) % 20 == 0:
            # Capture outputs
            t_surf = get_values_as_numpy(state["surface_temperature"]).flatten()[0]
            sh_flux = get_values_as_numpy(
                state["surface_upward_sensible_heat_flux"]
            ).flatten()[0]
            lh_flux = get_values_as_numpy(
                state["surface_upward_latent_heat_flux"]
            ).flatten()[0]
            output_data.append((i, t_surf, sh_flux, lh_flux))

        state["time"] += timestep
        set_values(state["eastward_wind"], 3.0)

    end_time = time.perf_counter()
    duration = end_time - start_time

    return duration, output_data


if __name__ == "__main__":
    # Run DataArrayBackend
    print("\n--- Running DataArrayBackend ---")
    try:
        dur_da, res_da = run_gmd_simulation(DataArrayBackend(), iterations=1000)
        print(f"DataArrayBackend Duration: {dur_da:.4f} s")
    except Exception as e:
        print(f"DataArrayBackend Failed: {e}")
        import traceback

        traceback.print_exc()
        dur_da, res_da = None, []

    # Run UnytBackend
    print("\n--- Running UnytBackend ---")
    try:
        dur_unyt, res_unyt = run_gmd_simulation(UnytBackend(), iterations=1000)
        print(f"UnytBackend Duration: {dur_unyt:.4f} s")
    except Exception as e:
        print(f"UnytBackend Failed: {e}")
        import traceback

        traceback.print_exc()
        dur_unyt, res_unyt = None, []

    # Compare
    if dur_da and dur_unyt:
        print(f"\nSpeedup: {dur_da / dur_unyt:.2f}x")

        # Compare results
        print("\nVerifying outputs match...")
        all_match = True
        if len(res_da) != len(res_unyt):
            print(f"Length mismatch: {len(res_da)} vs {len(res_unyt)}")
            all_match = False
        else:
            for i, (row_da, row_unyt) in enumerate(zip(res_da, res_unyt)):
                # Row format: (iter, t_surf, sh_flux, lh_flux)
                if row_da[0] != row_unyt[0]:
                    print(f"Iteration mismatch at index {i}")
                    all_match = False
                    break

                # Compare floats
                vals_da = np.array(row_da[1:])
                vals_unyt = np.array(row_unyt[1:])
                if not np.allclose(vals_da, vals_unyt, rtol=1e-5):
                    print(f"Mismatch at iter {row_da[0]}:")
                    print(f"  DA:   {vals_da}")
                    print(f"  Unyt: {vals_unyt}")
                    all_match = False
                    break

        if all_match:
            print("SUCCESS: Outputs match!")
        else:
            print("FAILURE: Outputs differ.")

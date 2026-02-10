import time
from datetime import timedelta

import numpy as np
import sympl
from sympl._core.backend import DataArrayBackend
from gfs_dynamical_core import GFSDynamicalCore

import climt
from climt import (
    HeldSuarez,
    UnytBackend,
    get_default_state,
    get_grid,
)


def set_values(obj, val):
    if hasattr(obj, "values"):  # DataArray
        obj.values[:] = val
    elif hasattr(obj, "data"):  # UnytStateContainer
        obj.data[:] = val


def get_values_as_numpy(obj):
    if hasattr(obj, "values"):  # DataArray
        return np.array(obj.values)
    elif hasattr(obj, "data"):  # UnytStateContainer
        return np.array(obj.data)  # Strips units
    return np.array(obj)


def run_held_suarez_simulation(backend_obj, iterations=100):
    print(f"Initializing for {type(backend_obj).__name__}...")
    sympl.set_backend(backend_obj)

    timestep = timedelta(seconds=600)
    if isinstance(backend_obj, UnytBackend):
        from climt import UnytTimeDelta
        timestep = UnytTimeDelta(timestep)

    grid = get_grid(nx=64, ny=32) # Smaller grid for faster benchmark
    held_suarez = HeldSuarez()
    dycore = GFSDynamicalCore([held_suarez])

    # Note: get_default_state must be called AFTER set_backend to get correct state objects
    state = get_default_state([dycore], grid_state=grid)

    # Set initial conditions
    # We use a fixed seed for reproducibility when comparing backends
    rng = np.random.default_rng(42)
    initial_wind = rng.standard_normal(get_values_as_numpy(state["eastward_wind"]).shape)
    set_values(state["eastward_wind"], initial_wind)

    output_data = []

    start_time = time.perf_counter()

    for i in range(iterations):
        # GFSDynamicalCore call
        diagnostics, new_state = dycore(state, timestep)
        state.update(diagnostics)
        state.update(new_state)

        if (i + 1) % 10 == 0:
            # Capture outputs: Mean zonal wind and mean temperature
            u_mean = np.mean(get_values_as_numpy(state["eastward_wind"]))
            t_mean = np.mean(get_values_as_numpy(state["air_temperature"]))
            u_max = np.max(get_values_as_numpy(state["eastward_wind"]))
            output_data.append((i, u_mean, t_mean, u_max))

        state["time"] += timestep

    end_time = time.perf_counter()
    duration = end_time - start_time

    return duration, output_data


if __name__ == "__main__":
    iterations = 100
    
    # Run DataArrayBackend
    print("\n--- Running DataArrayBackend ---")
    try:
        dur_da, res_da = run_held_suarez_simulation(DataArrayBackend(), iterations=iterations)
        print(f"DataArrayBackend Duration: {dur_da:.4f} s")
    except Exception as e:
        print(f"DataArrayBackend Failed: {e}")
        import traceback
        traceback.print_exc()
        dur_da, res_da = None, []

    # Run UnytBackend
    print("\n--- Running UnytBackend ---")
    try:
        dur_unyt, res_unyt = run_held_suarez_simulation(UnytBackend(), iterations=iterations)
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
                # Row format: (iter, u_mean, t_mean, u_max)
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
                    # break # Check all mismatches? No, break is fine for now

        if all_match:
            print("SUCCESS: Outputs match!")
        else:
            print("FAILURE: Outputs differ.")

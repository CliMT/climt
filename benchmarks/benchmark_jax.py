import time
from datetime import timedelta
import numpy as np
import sympl
from sympl._core.backend import DataArrayBackend

import climt
from climt import (
    HeldSuarez,
    UnytBackend,
    JaxBackend,
    get_default_state,
    get_grid,
)
from gfs_dynamical_core import GFSDynamicalCore

def run_simulation(backend_obj, iterations=50):
    print(f"\nInitializing for {type(backend_obj).__name__}...")
    sympl.set_backend(backend_obj)

    timestep = timedelta(seconds=600)
    if isinstance(backend_obj, (climt.UnytBackend, climt.JaxBackend)):
        if isinstance(backend_obj, climt.UnytBackend):
            from climt import UnytTimeDelta as TimeDelta
        else:
            from climt import JaxTimeDelta as TimeDelta
        timestep = TimeDelta(timestep)

    grid = get_grid(nx=64, ny=32)
    held_suarez = HeldSuarez()
    dycore = GFSDynamicalCore([held_suarez])

    state = get_default_state([dycore], grid_state=grid)

    # Use fixed seed for reproducibility
    rng = np.random.default_rng(42)
    u_shape = state["eastward_wind"].shape
    initial_u = rng.standard_normal(u_shape)
    
    # Set values. For JAX, we might need a different approach if it's strictly immutable,
    # but JaxStateContainer.values returns a numpy view, so this might work if we're lucky
    # or if we implement create_quantity properly.
    if hasattr(state["eastward_wind"], "values"):
        try:
            state["eastward_wind"].values[:] = initial_u
        except (TypeError, ValueError):
            # If immutable, we replace the quantity
            state["eastward_wind"] = sympl.get_backend().create_quantity(
                initial_u, "eastward_wind", "m/s", state["eastward_wind"].dims
            )
    else:
        state["eastward_wind"].data[:] = initial_u

    output_data = []
    start_time = time.perf_counter()

    for i in range(iterations):
        diagnostics, new_state = dycore(state, timestep)
        state.update(diagnostics)
        state.update(new_state)

        if (i + 1) % 10 == 0:
            u_max = np.max(np.asarray(state["eastward_wind"]))
            t_mean = np.mean(np.asarray(state["air_temperature"]))
            output_data.append((i, u_max, t_mean))

        state["time"] += timestep

    end_time = time.perf_counter()
    return end_time - start_time, output_data

if __name__ == "__main__":
    iterations = 1000
    
    # DataArray
    dur_da, res_da = run_simulation(DataArrayBackend(), iterations)
    print(f"DataArrayBackend: {dur_da:.4f}s")

    # Unyt
    dur_unyt, res_unyt = run_simulation(UnytBackend(), iterations)
    print(f"UnytBackend:      {dur_unyt:.4f}s")

    # Jax
    try:
        dur_jax, res_jax = run_simulation(JaxBackend(), iterations)
        print(f"JaxBackend:       {dur_jax:.4f}s")
    except Exception as e:
        print(f"JaxBackend Failed: {e}")
        import traceback
        traceback.print_exc()
        dur_jax, res_jax = None, []

    # Comparison
    if res_da and res_jax:
        print("\nVerifying JAX vs DataArray...")
        match = True
        for r_da, r_jax in zip(res_da, res_jax):
            if not np.allclose(r_da[1:], r_jax[1:], rtol=1e-5):
                print(f"Mismatch at iter {r_da[0]}: DA={r_da[1:]}, JAX={r_jax[1:]}")
                match = False
                break
        if match:
            print("SUCCESS: JAX outputs match DataArray!")

# -*- coding: utf-8 -*-
import climt
from sympl import TimeDifferencingWrapper
from gfs_dynamical_core import GFSDynamicalCore
import numpy as np
from datetime import timedelta
import time
import os

# Use EmanuelConvectionPython to test Numba optimization
from climt import EmanuelConvectionPython, GrayLongwaveRadiation, SlabSurface, SimplePhysics

def run_benchmark(steps=20):
    climt.set_constants_from_dict({
        'stellar_irradiance': {'value': 200, 'units': 'W m^-2'}})

    model_time_step = timedelta(seconds=600)

    # Create components
    convection = EmanuelConvectionPython()
    simple_physics = TimeDifferencingWrapper(SimplePhysics())
    radiation = GrayLongwaveRadiation()
    slab_surface = SlabSurface()

    dycore = GFSDynamicalCore(
        [simple_physics, radiation, slab_surface,
         convection], number_of_damped_levels=5
    )
    # Smaller grid for faster benchmark if needed, 
    # but let's keep the original size to see real impact
    grid = climt.get_grid(nx=128, ny=62)

    # Create model state
    my_state = climt.get_default_state([dycore], grid_state=grid)

    # Set initial/boundary conditions
    latitudes = my_state['latitude'].values
    sw_flux_equator = 300
    sw_flux_pole = 0

    sw_flux_profile = sw_flux_equator - (
        (sw_flux_equator - sw_flux_pole)*(np.sin(np.radians(latitudes))**2))

    my_state['downwelling_shortwave_flux_in_air'].values[:] = sw_flux_profile[np.newaxis, :]
    my_state['surface_temperature'].values[:] = 290.
    my_state['ocean_mixed_layer_thickness'].values[:] = 5
    my_state['eastward_wind'].values[:] = np.random.randn(*my_state['eastward_wind'].shape)

    print(f"Starting benchmark for {steps} steps...")
    
    # Warmup step (to exclude JIT compilation time from main loop if desired, 
    # but usually we want to see it or just run enough steps)
    dycore(my_state, model_time_step)
    
    start_time = time.time()
    for i in range(steps):
        diag, my_state = dycore(my_state, model_time_step)
        my_state.update(diag)
        my_state['time'] += model_time_step
    
    end_time = time.time()
    total_time = end_time - start_time
    avg_time = total_time / steps
    
    print(f"Benchmark finished.")
    print(f"Total time for {steps} steps: {total_time:.2f}s")
    print(f"Average time per step: {avg_time:.4f}s")
    
    return total_time

if __name__ == "__main__":
    numba_disabled = os.environ.get("NUMBA_DISABLE_JIT", "0") == "1"
    mode = "WITHOUT Numba" if numba_disabled else "WITH Numba"
    print(f"Running in mode: {mode}")
    run_benchmark(steps=20)

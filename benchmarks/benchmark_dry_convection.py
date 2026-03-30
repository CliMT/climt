import time
import numpy as np
import climt
from climt import DryConvectiveAdjustment, get_grid, get_default_state
import sympl
from datetime import timedelta

def benchmark_dry_convection(ncol=1024, nlev=30, iterations=10):
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    dc = DryConvectiveAdjustment()
    
    # NumPy path
    sympl.set_backend(sympl.DataArrayBackend())
    state = get_default_state([dc], grid_state=grid)
    
    # Create an unstable profile
    # Temperature decreasing too fast with height (super-adiabatic)
    unstable_level = 5
    state['air_temperature'].values[:unstable_level, :, :] += 20.0
    state['specific_humidity'].values[:unstable_level, :, :] = 0.02
    
    timestep = timedelta(minutes=10)
    
    # Run once
    dc(state, timestep)
    
    print(f"Benchmarking DryConvectiveAdjustment with {ncol} columns, {iterations} iterations")

    start = time.perf_counter()
    for _ in range(iterations):
        dc(state, timestep)
    end = time.perf_counter()
    print(f"Dry Convection NumPy: {end - start:.4f}s")

    # JAX path
    try:
        from climt import JaxBackend
        import jax
        jax.config.update("jax_enable_x64", True)
        jax.config.update('jax_platform_name', 'cpu')
        sympl.set_backend(JaxBackend())
        state_jax = get_default_state([dc], grid_state=grid)
        
        # Unstable profile for JAX
        state_jax['air_temperature'].data = state_jax['air_temperature'].data.at[:unstable_level].add(20.0)
        state_jax['specific_humidity'].data = state_jax['specific_humidity'].data.at[:unstable_level].set(0.02)

        # Run once to JIT
        dc(state_jax, timestep)
        
        start = time.perf_counter()
        for _ in range(iterations):
            dc(state_jax, timestep)
        end = time.perf_counter()
        print(f"Dry Convection JAX:   {end - start:.4f}s")
    except (ImportError, Exception) as e:
        print(f"JAX path failed or not available: {e}")

if __name__ == "__main__":
    benchmark_dry_convection()

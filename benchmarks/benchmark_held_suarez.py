import time
import numpy as np
import climt
from climt import HeldSuarez, get_grid, get_default_state
import sympl

def benchmark_held_suarez(ncol=8192, nlev=30, iterations=100):
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    hs = HeldSuarez()
    
    # NumPy/Numba path
    sympl.set_backend(sympl.DataArrayBackend())
    state = get_default_state([hs], grid_state=grid)
    
    # Set some non-zero winds and varied temperature
    state['eastward_wind'].values[:] = 10.0
    state['air_temperature'].values[:] = 300.0
    
    # Run once to JIT
    hs(state)
    
    print(f"Benchmarking Held-Suarez with {ncol} columns, {iterations} iterations")

    start = time.perf_counter()
    for _ in range(iterations):
        hs(state)
    end = time.perf_counter()
    print(f"Held-Suarez Numba: {end - start:.4f}s")

    # JAX path
    try:
        from climt import JaxBackend
        import jax
        jax.config.update("jax_enable_x64", True)
        jax.config.update('jax_platform_name', 'cpu')
        sympl.set_backend(JaxBackend())
        state_jax = get_default_state([hs], grid_state=grid)
        
        # Run once to JIT
        hs(state_jax)
        
        start = time.perf_counter()
        for _ in range(iterations):
            hs(state_jax)
        end = time.perf_counter()
        print(f"Held-Suarez JAX:   {end - start:.4f}s")
    except ImportError:
        print("JAX not available")

if __name__ == "__main__":
    benchmark_held_suarez()

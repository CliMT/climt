import time
import numpy as np
import climt
from climt import GridScaleCondensation, get_grid, get_default_state
import sympl
from datetime import timedelta

def benchmark_grid_scale_condensation(ncol=8192, nlev=30, iterations=100):
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    gsc = GridScaleCondensation()
    
    # NumPy path
    sympl.set_backend(sympl.DataArrayBackend())
    state = get_default_state([gsc], grid_state=grid)
    
    # Create saturated profile
    state['air_temperature'].values[:] = 280.0
    state['specific_humidity'].values[:] = 0.05 # Very moist
    
    timestep = timedelta(minutes=10)
    
    # Run once
    gsc(state, timestep)
    
    print(f"Benchmarking GridScaleCondensation with {ncol} columns, {iterations} iterations")

    start = time.perf_counter()
    for _ in range(iterations):
        gsc(state, timestep)
    end = time.perf_counter()
    print(f"Grid Scale Condensation NumPy: {end - start:.4f}s")

    # JAX path
    try:
        from climt import JaxBackend
        import jax
        import jax.numpy as jnp
        jax.config.update("jax_enable_x64", True)
        jax.config.update('jax_platform_name', 'cpu')
        sympl.set_backend(JaxBackend())
        state_jax = get_default_state([gsc], grid_state=grid)
        
        # Saturated profile for JAX
        state_jax['air_temperature'].data = jnp.full(state_jax['air_temperature'].shape, 280.0)
        state_jax['specific_humidity'].data = jnp.full(state_jax['specific_humidity'].shape, 0.05)

        # Run once to JIT
        gsc(state_jax, timestep)
        
        start = time.perf_counter()
        for _ in range(iterations):
            gsc(state_jax, timestep)
        end = time.perf_counter()
        print(f"Grid Scale Condensation JAX:   {end - start:.4f}s")
    except (ImportError, Exception) as e:
        import traceback
        # traceback.print_exc()
        print(f"JAX path failed or not available: {e}")

if __name__ == "__main__":
    benchmark_grid_scale_condensation()

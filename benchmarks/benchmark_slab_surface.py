import time
import numpy as np
import climt
from climt import SlabSurface, get_grid, get_default_state
import sympl
from datetime import timedelta

def benchmark_slab_surface(ncol=8192, iterations=100):
    grid = get_grid(nx=ncol, ny=1, nz=20)
    slab = SlabSurface()
    
    # NumPy path
    sympl.set_backend(sympl.DataArrayBackend())
    state = get_default_state([slab], grid_state=grid)
    
    # Run once
    slab(state)
    
    print(f"Benchmarking SlabSurface with {ncol} columns, {iterations} iterations")

    start = time.perf_counter()
    for _ in range(iterations):
        slab(state)
    end = time.perf_counter()
    print(f"SlabSurface NumPy: {end - start:.4f}s")

    # JAX path
    try:
        from climt import JaxBackend
        import jax
        import jax.numpy as jnp
        jax.config.update("jax_enable_x64", True)
        jax.config.update('jax_platform_name', 'cpu')
        sympl.set_backend(JaxBackend())
        state_jax = get_default_state([slab], grid_state=grid)
        
        # Run once to JIT
        slab(state_jax)
        
        start = time.perf_counter()
        for _ in range(iterations):
            slab(state_jax)
        end = time.perf_counter()
        print(f"SlabSurface JAX:   {end - start:.4f}s")
    except (ImportError, Exception) as e:
        print(f"JAX path failed or not available: {e}")

if __name__ == "__main__":
    benchmark_slab_surface()

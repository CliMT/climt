import time
import numpy as np
import climt
from climt import Instellation, get_grid, get_default_state
import sympl
from datetime import datetime

def benchmark_instellation(ncol=8192, iterations=100):
    grid = get_grid(nx=ncol, ny=1, nz=10)
    inst = Instellation()
    
    # NumPy path
    sympl.set_backend(sympl.DataArrayBackend())
    state = get_default_state([inst], grid_state=grid)
    state['time'] = datetime(2026, 6, 21, 12, 0)
    
    # Run once
    inst(state)
    
    print(f"Benchmarking Instellation with {ncol} columns, {iterations} iterations")

    start = time.perf_counter()
    for _ in range(iterations):
        inst(state)
    end = time.perf_counter()
    print(f"Instellation NumPy: {end - start:.4f}s")

    # JAX path
    try:
        from climt import JaxBackend
        import jax
        import jax.numpy as jnp
        jax.config.update("jax_enable_x64", True)
        jax.config.update('jax_platform_name', 'cpu')
        sympl.set_backend(JaxBackend())
        state_jax = get_default_state([inst], grid_state=grid)
        state_jax['time'] = datetime(2026, 6, 21, 12, 0)

        # Run once to JIT
        inst(state_jax)
        
        start = time.perf_counter()
        for _ in range(iterations):
            inst(state_jax)
        end = time.perf_counter()
        print(f"Instellation JAX:   {end - start:.4f}s")
    except (ImportError, Exception) as e:
        print(f"JAX path failed or not available: {e}")

if __name__ == "__main__":
    benchmark_instellation()

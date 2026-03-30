import time
import numpy as np
import climt
from climt import GrayLongwaveRadiation, Frierson06LongwaveOpticalDepth, get_grid, get_default_state
import sympl

def benchmark_gray_radiation(ncol=8192, nlev=30, iterations=50):
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    tau_comp = Frierson06LongwaveOpticalDepth()
    lw_comp = GrayLongwaveRadiation()
    
    # NumPy/Numba path
    sympl.set_backend(sympl.DataArrayBackend())
    state = get_default_state([tau_comp, lw_comp], grid_state=grid)
    
    # Run once to JIT
    tau_comp(state)
    lw_comp(state)
    
    print(f"Benchmarking with {ncol} columns, {iterations} iterations")

    start = time.perf_counter()
    for _ in range(iterations):
        tau_comp(state)
    end = time.perf_counter()
    print(f"Frierson Numba: {end - start:.4f}s")

    start = time.perf_counter()
    for _ in range(iterations):
        lw_comp(state)
    end = time.perf_counter()
    print(f"Gray Radiation Numba: {end - start:.4f}s")

    # JAX path
    try:
        from climt import JaxBackend
        import jax
        jax.config.update("jax_enable_x64", True)
        jax.config.update('jax_platform_name', 'cpu')
        sympl.set_backend(JaxBackend())
        state_jax = get_default_state([tau_comp, lw_comp], grid_state=grid)
        
        # Run once to JIT
        tau_comp(state_jax)
        lw_comp(state_jax)
        
        start = time.perf_counter()
        for _ in range(iterations):
            tau_comp(state_jax)
        end = time.perf_counter()
        print(f"Frierson JAX:   {end - start:.4f}s")

        start = time.perf_counter()
        for _ in range(iterations):
            lw_comp(state_jax)
        end = time.perf_counter()
        print(f"Gray Rad JAX:   {end - start:.4f}s")
    except ImportError:
        print("JAX not available")

if __name__ == "__main__":
    benchmark_gray_radiation()

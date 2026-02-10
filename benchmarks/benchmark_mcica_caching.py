import time
from datetime import timedelta
import numpy as np
import sympl
from sympl._core.backend import DataArrayBackend
import climt
from climt import RRTMGLongwave, RRTMGShortwave, get_default_state, get_grid

def run_benchmark(component_cls, mcica=True, iterations=20):
    print(f"Benchmarking {component_cls.__name__} (mcica={mcica}, iterations={iterations})...")
    
    sympl.set_backend(DataArrayBackend())
    
    # Larger grid to make allocation overhead more visible
    grid = get_grid(nx=64, ny=64, nz=30)
    
    component = component_cls(mcica=mcica)
    state = get_default_state([component], grid_state=grid)
    
    timestep = timedelta(minutes=20)
    
    # Warm up
    component(state)
    
    start_time = time.perf_counter()
    for _ in range(iterations):
        component(state)
        state['time'] += timestep
    end_time = time.perf_counter()
    
    duration = end_time - start_time
    print(f"  Duration: {duration:.4f} s")
    return duration

if __name__ == "__main__":
    # Test LW
    t1 = run_benchmark(RRTMGLongwave, mcica=True, iterations=20)
    
    # Test SW
    t2 = run_benchmark(RRTMGShortwave, mcica=True, iterations=20)
    
    print("\nBenchmark complete.")

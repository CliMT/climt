import time
from datetime import timedelta
import numpy as np
import sympl
from sympl._core.backend import DataArrayBackend
import climt
from climt import EmanuelConvection, EmanuelConvectionPython, get_default_state, get_grid

def run_benchmark():
    print("Benchmarking Emanuel Convection: Fortran vs Pure Python...")
    
    sympl.set_backend(DataArrayBackend())
    grid = get_grid(nx=32, ny=32, nz=30)
    
    # Initialize components
    conv_fortran = EmanuelConvection()
    conv_python = EmanuelConvectionPython()
    
    # Initial state
    state = get_default_state([conv_fortran], grid_state=grid)
    # Ensure specific humidity is reasonable for convection
    state['specific_humidity'].values[:] = 0.01
    state['air_temperature'].values[:] = 290.0
    
    timestep = timedelta(minutes=20)
    
    # Run Fortran version
    start = time.perf_counter()
    t_fort, d_fort = conv_fortran(state, timestep)
    dur_fort = time.perf_counter() - start
    print(f"  Fortran Duration: {dur_fort:.4f} s")
    
    # Run Python version
    start = time.perf_counter()
    t_py, d_py = conv_python(state, timestep)
    dur_py = time.perf_counter() - start
    print(f"  Python Duration:  {dur_py:.4f} s")
    
    # Compare outputs
    print("\nVerifying Outputs...")
    t_vars = ["air_temperature", "specific_humidity"]
    for var in t_vars:
        diff = np.abs(t_fort[var].values - t_py[var].values)
        print(f"  {var} max diff: {np.max(diff):.2e}")

    print(f"\nSpeedup (Fortran/Python): {dur_fort/dur_py:.2f}x")

if __name__ == "__main__":
    run_benchmark()

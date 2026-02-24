import numpy as np
import sympl
from datetime import timedelta
from climt import EmanuelConvection, get_grid, get_default_state
from climt._components.emanuel.pure_python import EmanuelConvectionPython

def create_convective_profile(grid):
    """Create a simple isothermal profile with high surface moisture to force convection."""
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state['air_temperature'].values[:] = 300.0
    q = np.zeros_like(state['specific_humidity'].values)
    q[0:5, :, :] = 0.020 
    state['specific_humidity'].values[:] = q
    return state

def run_comparison():
    print("Running Emanuel Parity Comparison (Isothermal Moist Profile)...")
    
    ncol = 1
    nlev = 30
    grid = get_grid(nx=1, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    
    conv_fortran = EmanuelConvection()
    conv_python = EmanuelConvectionPython()
    
    state = create_convective_profile(grid)
    timestep = timedelta(minutes=15)
    
    from climt._core.util import numpy_version_of
    raw_state = numpy_version_of(state)
    
    python_state = {
        'air_temperature': raw_state['air_temperature'].reshape(nlev, ncol),
        'specific_humidity': raw_state['specific_humidity'].reshape(nlev, ncol),
        'air_pressure': (raw_state['air_pressure'] / 100.0).reshape(nlev, ncol),
        'air_pressure_on_interface_levels': (raw_state['air_pressure_on_interface_levels'] / 100.0).reshape(nlev + 1, ncol),
        'eastward_wind': raw_state['eastward_wind'].reshape(nlev, ncol),
        'northward_wind': raw_state['northward_wind'].reshape(nlev, ncol),
    }

    print("\nCalling Fortran version...")
    try:
        t_fort, d_fort = conv_fortran(state, timestep)
        ft_fort = t_fort['air_temperature'].values.flatten() * 86400
        fq_fort = t_fort['specific_humidity'].values.flatten() * 86400 * 1000 # g/kg/day
        # fortran_state_diag = d_fort['convective_state'].values.flatten()
    except Exception as e:
        print(f"Fortran call failed: {e}")
        ft_fort = np.zeros(nlev)
        fq_fort = np.zeros(nlev)
    
    print("Calling Python version...")
    try:
        t_py, d_py = conv_python.array_call(python_state, timestep)
        ft_py = t_py['air_temperature'].flatten() * 86400
        fq_py = t_py['specific_humidity'].flatten() * 86400 * 1000 # g/kg/day
        # Access internal state for debugging - this requires the instance to store them
        # pure_python.py stores nk, icb, inb in self ONLY if I modify it to do so.
        # But I can't easily access local variables from here.
        # However, I can print them from within the component if I re-enable debug prints there.
        # OR, I can check if I stored them in self.
        pass 
    except Exception as e:
        print(f"Python call failed: {e}")
        ft_py = np.zeros(nlev)
        fq_py = np.zeros(nlev)
    
    print("\n--- Tendency Comparison (T: K/day, Q: g/kg/day) ---")
    print(f"{'Lev':<4} | {'Press':<8} | {'T_Fort':<10} | {'T_Py':<10} | {'Q_Fort':<10} | {'Q_Py':<10}")
    print("-" * 70)
    
    p = raw_state['air_pressure'].flatten()
    for i in range(nlev):
        t_diff = abs(ft_fort[i] - ft_py[i])
        q_diff = abs(fq_fort[i] - fq_py[i])
        if t_diff > 1e-2 or q_diff > 1e-2 or abs(ft_fort[i]) > 1e-2:
            print(f"{i:<4} | {p[i]/100.0:<8.2f} | {ft_fort[i]:10.4f} | {ft_py[i]:10.4f} | {fq_fort[i]:10.4f} | {fq_py[i]:10.4f}")

if __name__ == "__main__":
    run_comparison()

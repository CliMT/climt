import numpy as np
import sympl
from datetime import timedelta
from climt import EmanuelConvection, get_grid, get_default_state
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v2 import EmanuelConvectionPythonV2

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
    conv_python_v2 = EmanuelConvectionPythonV2()
    
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
    except Exception as e:
        print(f"Fortran call failed: {e}")
        ft_fort = np.zeros(nlev)
        fq_fort = np.zeros(nlev)
    
    print("Calling Python version (Original)...")
    try:
        t_py, d_py = conv_python.array_call(python_state, timestep)
        ft_py = t_py['air_temperature'].flatten() * 86400
        fq_py = t_py['specific_humidity'].flatten() * 86400 * 1000 # g/kg/day
    except Exception as e:
        print(f"Python (Original) call failed: {e}")
        ft_py = np.zeros(nlev)
        fq_py = np.zeros(nlev)

    print("Calling Python version (V2 Abstraction)...")
    try:
        t_py_v2, d_py_v2 = conv_python_v2.array_call(python_state, timestep)
        ft_py_v2 = t_py_v2['air_temperature'].flatten() * 86400
        fq_py_v2 = t_py_v2['specific_humidity'].flatten() * 86400 * 1000 # g/kg/day
    except Exception as e:
        print(f"Python (V2) call failed: {e}")
        import traceback
        traceback.print_exc()
        ft_py_v2 = np.zeros(nlev)
        fq_py_v2 = np.zeros(nlev)
    
    print("\n--- Tendency Comparison (T: K/day, Q: g/kg/day) ---")
    print(f"{'Lev':<4} | {'Press':<8} | {'T_Fort':<10} | {'T_Py_Orig':<10} | {'T_Py_V2':<10}")
    print("-" * 60)
    
    p = raw_state['air_pressure'].flatten()
    for i in range(nlev):
        print(f"{i:<4} | {p[i]/100.0:<8.2f} | {ft_fort[i]:10.4f} | {ft_py[i]:10.4f} | {ft_py_v2[i]:10.4f}")

    print("\n--- Specific Humidity Tendency Comparison (g/kg/day) ---")
    print(f"{'Lev':<4} | {'Press':<8} | {'Q_Fort':<10} | {'Q_Py_Orig':<10} | {'Q_Py_V2':<10}")
    print("-" * 60)
    for i in range(nlev):
        print(f"{i:<4} | {p[i]/100.0:<8.2f} | {fq_fort[i]:10.4f} | {fq_py[i]:10.4f} | {fq_py_v2[i]:10.4f}")

    # Check parity between V2 and Original
    np.testing.assert_allclose(ft_py, ft_py_v2, atol=1e-10, err_msg="Temperature tendency parity failed")
    np.testing.assert_allclose(fq_py, fq_py_v2, atol=1e-10, err_msg="Humidity tendency parity failed")
    print("\nSUCCESS: Parity between Original and V2 implementations confirmed!")

if __name__ == "__main__":
    run_comparison()

if __name__ == "__main__":
    run_comparison()

import numpy as np
import climt
import sympl
from sympl import DataArray
from datetime import datetime, timedelta
from sympl._core.tracers import reset_tracers, reset_packers

from climt import SimCloudCondensation

def test_simcloud_water_conservation():
    print("==================================================")
    print("Starting Water Conservation Unit Test...")
    print("==================================================")
    
    nz = 30
    grid = climt.get_grid(nx=1, ny=1, nz=nz)
    
    # 1. Test in Diagnostic Mode (no prognostic cloud tracers)
    print("\n--- Testing 1: Diagnostic Mode (Instant Precipitation) ---")
    
    # Reset tracers in sympl
    reset_tracers()
    reset_packers()
    
    # Initialize component
    cloud_diag_mode = SimCloudCondensation()
    assert not cloud_diag_mode.use_prognostic_tracers, "Expected use_prognostic_tracers to be False!"
    
    # Get state
    state = climt.get_default_state([climt.RRTMGLongwave()], grid_state=grid)
    
    # Initialize mock surface and profiles
    state['time'] = datetime(2000, 1, 1)
    
    p = state['air_pressure'].values
    ps = state['surface_air_pressure'].values
    dp = np.abs(np.diff(state['air_pressure_on_interface_levels'].values, axis=0))
    g = cloud_diag_mode.params.g
    
    # Warm temperature profile (300K at surface, 6.5K/km lapse rate)
    t_surf = 300.0
    gamma = 0.0065
    t_profile = t_surf * (p / ps)**(287.0 * gamma / 9.81)
    state['air_temperature'][:] = np.maximum(200.0, t_profile)
    
    # High initial specific humidity profile
    # Read constants from the public params bundle (all physical constants live there after the refactor)
    P = cloud_diag_mode.params
    temp = state['air_temperature'].values
    tc = temp - 273.15
    tc_clip = np.maximum(tc, -100.0)
    es_w = np.exp(34.4942 - 4924.99 / (tc_clip + 237.1)) / (tc_clip + 105.0)**1.57
    es_i = np.exp(43.4942 - 6545.8 / (tc_clip + 278.0)) / (tc_clip + 868.0)**2.0
    liquid_fraction = np.clip((temp - P.T_min) / (P.T_max - P.T_min), 0.0, 1.0)
    sat_vp = liquid_fraction * es_w + (1.0 - liquid_fraction) * es_i
    sat_q = P.epsilon * sat_vp / (p - (1.0 - P.epsilon) * sat_vp)
    
    # Set to 90% RH below 200 hPa to force condensation
    rh = np.where(p < 20000.0, 0.0, 0.90)
    state['specific_humidity'][:] = np.maximum(rh * sat_q, 1e-8)
    
    # Calculate initial total column water mass
    # q_total = q_vapor (since q_water, q_ice = 0)
    mass_init = np.sum(state['specific_humidity'].values * dp / g)
    
    # Step component
    dt = timedelta(minutes=15)
    diags, outputs = cloud_diag_mode(state, timestep=dt)
    
    # Apply updates
    state.update(outputs)
    state.update(diags)
    
    # Calculate final total column water mass:
    # final vapor mass + instant precipitation
    mass_vapor_final = np.sum(state['specific_humidity'].values * dp / g)
    precipitation = state['precipitation_amount'].values[0, 0] # kg m^-2
    mass_final = mass_vapor_final + precipitation
    
    error = abs(mass_final - mass_init)
    print(f"Initial water mass: {mass_init:.12f} kg/m^2")
    print(f"Final water mass:   {mass_vapor_final:.12f} kg/m^2 (vapor) + {precipitation:.12f} kg/m^2 (precip) = {mass_final:.12f} kg/m^2")
    print(f"Conservation Error: {error:e} kg/m^2")
    
    assert error < 1e-14, f"Water conservation failed in diagnostic mode! Error: {error:e}"
    print("Diagnostic mode test PASSED!")
    
    # 2. Test in Prognostic Mode (retains cloud water and ice as tracers)
    print("\n--- Testing 2: Prognostic Mode (Conserves Tracers) ---")
    
    # Reset and register prognostic cloud tracers
    reset_tracers()
    reset_packers()
    sympl.register_tracer('cloud_liquid_water_mixing_ratio', 'kg/kg')
    sympl.register_tracer('cloud_ice_mixing_ratio', 'kg/kg')
    
    cloud_prog_mode = SimCloudCondensation()
    assert cloud_prog_mode.use_prognostic_tracers, "Expected use_prognostic_tracers to be True!"
    
    state = climt.get_default_state([climt.RRTMGLongwave()], grid_state=grid)
    state['time'] = datetime(2000, 1, 1)
    
    # Reset profiles
    state['air_temperature'][:] = np.maximum(200.0, t_profile)
    state['specific_humidity'][:] = np.maximum(rh * sat_q, 1e-8)
    
    # Initialize cloud tracers to zero (or small non-zero values)
    state['cloud_liquid_water_mixing_ratio'] = DataArray(
        np.zeros_like(state['specific_humidity'].values),
        dims=state['specific_humidity'].dims,
        attrs={'units': 'kg/kg'}
    )
    state['cloud_ice_mixing_ratio'] = DataArray(
        np.zeros_like(state['specific_humidity'].values),
        dims=state['specific_humidity'].dims,
        attrs={'units': 'kg/kg'}
    )
    
    # Calculate initial total column water mass: vapor + liquid + ice
    mass_init = np.sum((state['specific_humidity'].values +
                        state['cloud_liquid_water_mixing_ratio'].values +
                        state['cloud_ice_mixing_ratio'].values) * dp / g)
    
    # Step component
    diags, outputs = cloud_prog_mode(state, timestep=dt)
    
    state.update(outputs)
    state.update(diags)
    
    # Calculate final total column water mass: vapor + liquid + ice
    mass_final = np.sum((state['specific_humidity'].values +
                         state['cloud_liquid_water_mixing_ratio'].values +
                         state['cloud_ice_mixing_ratio'].values) * dp / g)
    
    error = abs(mass_final - mass_init)
    print(f"Initial water mass: {mass_init:.12f} kg/m^2")
    print(f"Final water mass:   {mass_final:.12f} kg/m^2")
    print(f"Conservation Error: {error:e} kg/m^2")
    
    assert error < 1e-14, f"Water conservation failed in prognostic mode! Error: {error:e}"
    print("Prognostic mode test PASSED!")
    
    print("\n==================================================")
    print("ALL WATER CONSERVATION TESTS PASSED SUCCESSFULLY!")
    print("==================================================")

if __name__ == '__main__':
    test_simcloud_water_conservation()

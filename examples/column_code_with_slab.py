from sympl import (
    DataArray, PlotFunctionMonitor, RelaxationPrognostic,
    AdamsBashforth)
from climt import SimplePhysics, get_default_state
import numpy as np
from datetime import timedelta

from climt import EmanuelConvection, RRTMGShortwave, RRTMGLongwave, SlabSurface


def get_interface_pressures(p, ps):
    """Given 3D pressure on model mid levels (cell centers) and the 2D surface
    pressure, return the 3D pressure on model full levels (cell interfaces).
    If the z-dimension of p is length K, the returned p_full will have a
    z-dimension of length K+1."""
    interface_pressures = np.zeros(
        (p.shape[0], p.shape[1], p.shape[2]+1), dtype=np.double)
    interface_pressures[:, :, 1:-1] = 0.5*(p[:, :, 1:] + p[:, :, :-1])
    interface_pressures[:, :, 0] = ps[:, :]
    return interface_pressures

def plot_function(fig, state):
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(
        state['convective_heating_rate'].to_units('degK day^-1').values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.plot(
        state['longwave_heating_rate'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.plot(
        state['shortwave_heating_rate'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.axes.invert_yaxis()
    # print(state['eastward_wind'].values.flatten())
    # print(state['air_pressure'].values.flatten())
    # ax.set_yscale('log')
    # ax.set_ylim(1e5, 100.)
    ax = fig.add_subplot(2, 2, 2)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.axes.invert_yaxis()
    # ax.set_ylim(1e5, 100.)

    ax = fig.add_subplot(2, 2, 3)
    ax.plot(
        state['eastward_wind'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.axes.invert_yaxis()
    # ax.set_ylim(1e5, 100.)
    ax = fig.add_subplot(2, 2, 4)
    net_flux = (state['upwelling_longwave_flux_in_air'] +
                state['upwelling_shortwave_flux_in_air'] -
                state['downwelling_longwave_flux_in_air'] -
                state['downwelling_shortwave_flux_in_air'])
    ax.plot(
        net_flux.values.flatten(),
        state['air_pressure_on_interface_levels'].to_units('mbar').values.flatten(), '-o')
    ax.axes.invert_yaxis()

monitor = PlotFunctionMonitor(plot_function)

convection = EmanuelConvection()
radiation_sw = RRTMGShortwave()
radiation_lw = RRTMGLongwave()
slab = SlabSurface()
simple_physics = SimplePhysics(use_external_surface_specific_humidity=False)

timestep = timedelta(minutes=5)

state = get_default_state([simple_physics, convection,
                           radiation_lw, radiation_sw, slab])

state['air_temperature'].values[:] = 250
state['surface_albedo_for_direct_shortwave'].values[:] = 0.5
state['surface_albedo_for_direct_near_infrared'].values[:] = 0.5
state['surface_albedo_for_diffuse_shortwave'].values[:] = 0.5

state['mole_fraction_of_ozone_in_air'].values[0, 0, :] = np.load('ozone_profile.npy')

state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'].loc[dict(mid_levels=18)] = 0.02
state['cloud_area_fraction_in_atmosphere_layer'].loc[dict(mid_levels=18)] = 1.

equilibrium_value = DataArray(
        np.ones((1, 1, len(state['air_pressure'])))*10.,
        dims=('x', 'y', 'mid_levels'),
        attrs={'units': 'm s^-1'})

tau = DataArray(
    np.array(2.), dims = [], attrs={'units': 'hour'})

relaxation = RelaxationPrognostic('eastward_wind', equilibrium_value, tau)
time_stepper = AdamsBashforth([relaxation, convection, radiation_lw, radiation_sw, slab])

import time
for i in range(60000):
    #print(i)
    # print('before:', state['atmosphere_convective_mass_flux'].values)
    print(state['surface_temperature'].values)
    convection.current_time_step = timestep
    diagnostics, state = time_stepper(state, timestep)
    # print('after: ', diagnostics['atmosphere_convective_mass_flux'].values)
    # print(diagnostics['convective_state'].values)
    state.update(diagnostics)
    new_state, diagnostics = simple_physics(state, timestep)
    # print(diagnostics['stratiform_precipitation_rate'].to_units('mm day^-1').values)
    # print('SH Flux:', diagnostics['surface_upward_sensible_heat_flux'].values.item())
    # print('LH Flux:', diagnostics['surface_upward_latent_heat_flux'].values.item())
    state.update(diagnostics)
    if i%3 == 0:
        monitor.store(state)
    state.update(new_state)
    state['eastward_wind'].values[:] = 6.
    # time.sleep(0.4  )

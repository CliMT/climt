from sympl import (
    DataArray, PlotFunctionMonitor,
    AdamsBashforth
)
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt

from climt import (
    EmanuelConvection, RRTMGShortwave, RRTMGLongwave, SlabSurface,
    DryConvectiveAdjustment, SimplePhysics, get_default_state
)


def plot_function(fig, state):
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(
        state['air_temperature_tendency_from_convection'].to_units('degK day^-1').values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.set_title('Conv. heating rate')
    ax.set_xlabel('K/day')
    ax.set_ylabel('millibar')
    ax.grid()

    ax.axes.invert_yaxis()
    ax = fig.add_subplot(2, 2, 2)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o')
    ax.set_title('Air temperature')
    ax.axes.invert_yaxis()
    ax.set_xlabel('K')
    ax.grid()

    ax = fig.add_subplot(2, 2, 3)
    ax.plot(
        state['air_temperature_tendency_from_longwave'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o',
        label='LW')
    ax.plot(
        state['air_temperature_tendency_from_shortwave'].values.flatten(),
        state['air_pressure'].to_units('mbar').values.flatten(), '-o',
        label='SW')
    ax.set_title('LW and SW Heating rates')
    ax.legend()
    ax.axes.invert_yaxis()
    ax.set_xlabel('K/day')
    ax.grid()
    ax.set_ylabel('millibar')

    ax = fig.add_subplot(2, 2, 4)
    net_flux = (state['upwelling_longwave_flux_in_air'] +
                state['upwelling_shortwave_flux_in_air'] -
                state['downwelling_longwave_flux_in_air'] -
                state['downwelling_shortwave_flux_in_air'])
    ax.plot(
        net_flux.values.flatten(),
        state['air_pressure_on_interface_levels'].to_units('mbar').values.flatten(), '-o')
    ax.set_title('Net Flux')
    ax.axes.invert_yaxis()
    ax.set_xlabel('W/m^2')
    ax.grid()
    plt.tight_layout()


monitor = PlotFunctionMonitor(plot_function)
timestep = timedelta(minutes=10)

convection = EmanuelConvection()
radiation_sw = RRTMGShortwave()
radiation_lw = RRTMGLongwave()
slab = SlabSurface()
simple_physics = SimplePhysics()
dry_convection = DryConvectiveAdjustment()

state = get_default_state(
    [simple_physics, convection, dry_convection,
     radiation_lw, radiation_sw, slab]
)

state['air_temperature'].values[:] = 290
state['surface_albedo_for_direct_shortwave'].values[:] = 0.4
state['surface_albedo_for_direct_near_infrared'].values[:] = 0.4
state['surface_albedo_for_diffuse_shortwave'].values[:] = 0.4

# Uncomment the following two lines to see how clouds change the radiative balance!

# state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 0.03
# state['cloud_area_fraction_in_atmosphere_layer'].loc[dict(mid_levels=slice(4, 8))] = 1.

state['zenith_angle'].values[:] = np.pi/2.5
state['surface_temperature'].values[:] = 300.
state['ocean_mixed_layer_thickness'].values[:] = 5
state['area_type'].values[:] = 'sea'

equilibrium_value = DataArray(
    np.ones(len(state['air_pressure']))*10.,
    dims=('mid_levels'),
    attrs={'units': 'm s^-1'})

tau = DataArray(
    np.array(2.), dims=[], attrs={'units': 'hour'})

time_stepper = AdamsBashforth([convection, radiation_lw, radiation_sw, slab])

for i in range(60000):
    diagnostics, state = time_stepper(state, timestep)
    state.update(diagnostics)

    diagnostics, new_state = simple_physics(state, timestep)
    state.update(diagnostics)
    # state.update(new_state)

    # diagnostics, new_state = dry_convection(state, timestep)
    # state.update(diagnostics)
    if (i+1) % 100 == 0:
        print('Surface Temperature: ', state['surface_temperature'].values.item())
        monitor.store(state)
    state.update(new_state)
    state['eastward_wind'].values[:] = 3.

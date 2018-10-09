from sympl import (
    PlotFunctionMonitor, AdamsBashforth, NetCDFMonitor
)
from climt import SimplePhysics, get_default_state
import numpy as np
from datetime import timedelta

from climt import EmanuelConvection, RRTMGShortwave, RRTMGLongwave, SlabSurface
import matplotlib.pyplot as plt


def plot_function(fig, state):
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(
        state['air_temperature_tendency_from_convection'].to_units(
            'degK day^-1').values.flatten(),
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
        state['air_pressure_on_interface_levels'].to_units(
            'mbar').values.flatten(), '-o')
    ax.set_title('Net Flux')
    ax.axes.invert_yaxis()
    ax.set_xlabel('W/m^2')
    ax.grid()
    plt.tight_layout()


monitor = PlotFunctionMonitor(plot_function)

timestep = timedelta(minutes=5)

convection = EmanuelConvection()
radiation_sw = RRTMGShortwave()
radiation_lw = RRTMGLongwave()
slab = SlabSurface()
simple_physics = SimplePhysics()

store_quantities = ['air_temperature',
                    'air_pressure',
                    'specific_humidity',
                    'air_pressure_on_interface_levels',
                    'air_temperature_tendency_from_convection',
                    'air_temperature_tendency_from_longwave',
                    'air_temperature_tendency_from_shortwave']
netcdf_monitor = NetCDFMonitor('rad_conv_eq.nc',
                               store_names=store_quantities,
                               write_on_store=True)
convection.current_time_step = timestep


state = get_default_state([simple_physics, convection,
                           radiation_lw, radiation_sw, slab])

state['air_temperature'].values[:] = 270
state['surface_albedo_for_direct_shortwave'].values[:] = 0.5
state['surface_albedo_for_direct_near_infrared'].values[:] = 0.5
state['surface_albedo_for_diffuse_shortwave'].values[:] = 0.5

state['zenith_angle'].values[:] = np.pi/2.5
state['surface_temperature'].values[:] = 300.
state['ocean_mixed_layer_thickness'].values[:] = 5
state['area_type'].values[:] = 'sea'

time_stepper = AdamsBashforth([convection, radiation_lw, radiation_sw, slab])

for i in range(20000):
    convection.current_time_step = timestep
    diagnostics, state = time_stepper(state, timestep)
    state.update(diagnostics)
    diagnostics, new_state = simple_physics(state, timestep)
    state.update(diagnostics)
    if (i+1) % 20 == 0:
        monitor.store(state)
        netcdf_monitor.store(state)
        print(i, state['surface_temperature'].values)
        print(state['surface_upward_sensible_heat_flux'])
        print(state['surface_upward_latent_heat_flux'])

    state.update(new_state)
    state['time'] += timestep
    state['eastward_wind'].values[:] = 3.

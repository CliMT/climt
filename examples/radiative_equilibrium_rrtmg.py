from sympl import (
    AdamsBashforth, PlotFunctionMonitor)
from climt import RRTMGShortwave, RRTMGLongwave, get_default_state
import numpy as np
from datetime import timedelta


def plot_function(fig, state):
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(
        state['shortwave_heating_rate'].values.flatten(),
        state['air_pressure'].values.flatten()/100, '-o', label='SW')
    ax.axes.invert_yaxis()
    ax.plot(
        state['longwave_heating_rate'].values.flatten(),
        state['air_pressure'].values.flatten()/100, '-o', label='LW')
    ax.axes.invert_yaxis()
    ax.set_title('Heating Rates')
    ax.grid()
    ax.set_xlabel('K/day')
    ax.set_ylabel('millibar')
    ax.legend()

    # ax.set_yscale('log')
    ax.set_ylim(1e3, 10.)
    ax = fig.add_subplot(1, 2, 2)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].values.flatten()/100, '-o')
    ax.axes.invert_yaxis()

    # ax.set_yscale('log')
    ax.set_ylim(1e3, 10.)
    ax.set_title('Temperature')
    ax.grid()
    ax.set_xlabel('K')
    ax.set_yticklabels([])


monitor = PlotFunctionMonitor(plot_function)
rad_sw = RRTMGShortwave()
rad_lw = RRTMGLongwave()
time_stepper = AdamsBashforth([rad_sw, rad_lw])
timestep = timedelta(hours=3)

mid_levels = {'label': 'mid_level',
              'values': np.arange(60),
              'units': ''}

int_levels = {'label': 'interface_level',
              'values': np.arange(61),
              'units': ''}
state = get_default_state([rad_sw, rad_lw], mid_levels=mid_levels, interface_levels=int_levels)

tp_profiles = np.load('thermodynamic_profiles.npz')
mol_profiles = np.load('molecule_profiles.npz')

state['air_pressure'].values[0, 0, :] = tp_profiles['air_pressure']
state['air_temperature'].values[0, 0, :] = tp_profiles['air_temperature']
state['air_pressure_on_interface_levels'].values[0, 0, :] = tp_profiles['interface_pressures']

state['specific_humidity'].values[0, 0, :] = mol_profiles['specific_humidity']*1e-3
state['mole_fraction_of_carbon_dioxide_in_air'].values[0, 0, :] = mol_profiles['carbon_dioxide']
state['mole_fraction_of_ozone_in_air'].values[0, 0, :] = mol_profiles['ozone']


for i in range(100000):

    # print(i)
    diagnostics, new_state = time_stepper(state, timestep)
    state.update(diagnostics)
    if i % 2 == 0:
        monitor.store(state)
    state = new_state

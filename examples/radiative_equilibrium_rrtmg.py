from sympl import (
    DataArray, AdamsBashforth, PlotFunctionMonitor)
from climt import RRTMLongwave, get_default_state
import numpy as np
from datetime import timedelta
import sympl

def get_interface_values(
    mid_level_values,
    surface_value,
    mid_level_pressure,
    interface_level_pressure):
    """Given 3D values of a quantity on model mid levels (cell centers) and the 2D surface
    value, return the 3D values of that quantity on model full levels (cell interfaces).
    If the z-dimension of p is length K, the returned p_full will have a
    z-dimension of length K+1.
    
    Routine borrowed from CESM (radiation.F90 in rrtmg folder)
    """

    interface_values = np.zeros(
        (mid_level_values.shape[0], mid_level_values.shape[1], mid_level_values.shape[2]+1),
        dtype=np.double, order='F')

    log_mid_p = np.log(mid_level_pressure)
    log_int_p = np.log(interface_level_pressure)

    interp_weight = (log_int_p[:,:,1:-1] - log_mid_p[:,:,1::])/(log_mid_p[:,:,:-1] -
                                                                log_mid_p[:,:,1::])

    interface_values[:,:,1:-1] = \
        mid_level_values[:,:,1::] - interp_weight*(mid_level_values[:,:,1::] -
                                                   mid_level_values[:,:,0:-1])

    interface_values[:, :, 0] = surface_value[:, :]
    #interface_values[:,:,-1] = 2*interface_values[:,:,-2] - interface_values[:,:,-3]
    interface_values[:,:,-1] = mid_level_values[:,:,-1]
    return interface_values


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    #print(state['air_temperature'].values.flatten())
    #print(state['air_pressure'].values.flatten())
    #ax.set_yscale('log')
    ax.set_ylim(1e5, 1.)

monitor = PlotFunctionMonitor(plot_function)
radiation = RRTMLongwave()
time_stepper = AdamsBashforth([radiation])
timestep = timedelta(hours=4)

z_information = {'label' : 'vertical_level',
                 'values' : np.arange(60),
                 'units' : ''}

state = get_default_state([radiation], z=z_information)
state['air_temperature_on_interface_levels'].values = get_interface_values(
            state['air_temperature'].values, state['surface_temperature'].values,
            state['air_pressure'].values, state['air_pressure_on_interface_levels'].values)

print(state['specific_humidity'])
state['specific_humidity'].values[:,:,:17] = 10.
print(state['specific_humidity'])
#format: p_mid, t_mid, p_int, t_int
tp_profiles = np.loadtxt('tp_profiles').transpose()
#format: h2o, co2, o3
mol_profiles = np.loadtxt('mol_profiles').transpose()

state['air_pressure'].values[0,0,:] = tp_profiles[0,:]*100
state['air_temperature'].values[0,0,:] = tp_profiles[1,:]
state['air_pressure_on_interface_levels'].values[0,0,1:] = tp_profiles[2,:]*100
state['air_temperature_on_interface_levels'].values[0,0,1:] = tp_profiles[3,:]
state['air_pressure_on_interface_levels'].values[0,0,0] = 101300.
state['air_temperature_on_interface_levels'].values[0,0,0] = 300.

state['specific_humidity'].values[0,0,:] = mol_profiles[0,:]*(18.02/28.964)*1000
state['carbon_dioxide_mixing_ratio'].values[0,0,:] = mol_profiles[1,:]
state['ozone_mixing_ratio'].values[0,0,:] = mol_profiles[2,:]


for i in range(8000):
    
    #print(i)
    diagnostics, new_state = time_stepper.__call__(state, timestep)
    state.update(diagnostics)
    new_state['air_temperature_on_interface_levels'].values = get_interface_values(
            new_state['air_temperature'].values, new_state['surface_temperature'].values,
            new_state['air_pressure'].values, new_state['air_pressure_on_interface_levels'].values)
    if i % 200 == 0:
        monitor.store(state)
    state = new_state

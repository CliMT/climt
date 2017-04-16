from sympl import (
    DataArray, PlotFunctionMonitor)
from climt import SimplePhysics
import numpy as np
from datetime import timedelta

pressure_axis = np.array(
    [1e5, 9.5e4, 9.2e4, 9e4, 8.8e4, 8.5e4, 8.2e4, 8e4, 7e4, 5e4, 3e4, 1e4, 8e3, 4e3, 1e3, 7e2, 4e2, 1e2,
     7., 4., 1.])


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


state = {
    'air_temperature': DataArray(
        np.ones((1, 1, len(pressure_axis)))*250.,
        dims=('x', 'y', 'mid_levels'),
        attrs={'units': 'degK'}),
    'northward_wind': DataArray(
        np.ones((1, 1, len(pressure_axis)))*2.,
        dims=('x', 'y', 'mid_levels'),
        attrs={'units': 'm s^-1'}),
    'eastward_wind': DataArray(
        np.ones((1, 1, len(pressure_axis)))*5.,
        dims=('x', 'y', 'mid_levels'),
        attrs={'units': 'm s^-1'}),
    'specific_humidity': DataArray(
        np.zeros((1, 1, len(pressure_axis))),
        dims=('x', 'y', 'mid_levels'),
        attrs={'units': 'kg kg^-1'})
}
constant_state = {
    'surface_temperature': DataArray(
        np.ones((1, 1))*300., dims=('x', 'y'), attrs={'units': 'degK'}),
    'surface_specific_humidity': DataArray(
        np.ones((1, 1)) * 0.01, dims=('x', 'y'), attrs={'units': 'g/g'}),
    'surface_air_pressure': DataArray(
        np.ones((1, 1))*1e5, dims=('x', 'y'), attrs={'units': 'Pa'}),
    'air_pressure': DataArray(
        pressure_axis[None, None, :], dims=('x', 'y', 'mid_levels'),
        attrs={'units': 'Pa'}),
    'latitude': DataArray(
        np.zeros((1,)), dims=('y',), attrs={'units': 'degrees_north'}),
}
interface_pressures = get_interface_pressures(
    constant_state['air_pressure'].values,
    constant_state['surface_air_pressure'].values)
interface_sigma = (
    interface_pressures/constant_state['surface_air_pressure'].values[:, :, None])
constant_state['air_pressure_on_interface_levels'] = DataArray(
    interface_pressures, dims=('x', 'y', 'interface_levels'),
    attrs={'units': 'Pa'})
constant_state['sigma_on_interface_levels'] = DataArray(
    interface_sigma, dims=('x', 'y', 'interface_levels'),
    attrs={'units': ''})

state.update(constant_state)


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['specific_humidity'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    # print(state['eastward_wind'].values.flatten())
    # print(state['air_pressure'].values.flatten())
    # ax.set_yscale('log')
    ax.set_ylim(1e5, 100.)


monitor = PlotFunctionMonitor(plot_function)
simple_physics = SimplePhysics(use_external_surface_specific_humidity=False)
timestep = timedelta(hours=1)

for i in range(600):
    #print(i)

    new_state, diagnostics = simple_physics(state, timestep)
    print('SH Flux:', diagnostics['surface_upward_sensible_heat_flux'].values.item())
    print('LH Flux:', diagnostics['surface_upward_latent_heat_flux'].values.item())
    state.update(diagnostics)
    monitor.store(state)
    state.update(new_state)

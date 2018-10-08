import climt
from sympl import (
    PlotFunctionMonitor,
    TimeDifferencingWrapper,
)
import numpy as np
from datetime import timedelta


def plot_function(fig, state):

    ax = fig.add_subplot(2, 2, 1)
    state['surface_temperature'].plot.contourf(
        ax=ax, levels=16, robust=True)
    ax.set_title('Surface Temperature')

    ax = fig.add_subplot(2, 2, 3)
    state['eastward_wind'].mean(dim='lon').plot.contourf(
        ax=ax, levels=16, robust=True)
    ax.set_title('Zonal Wind')

    ax = fig.add_subplot(2, 2, 2)
    state['air_temperature_tendency_from_convection'].transpose().mean(
        dim='lon').plot.contourf(
        ax=ax, levels=16, robust=True)
    ax.set_title('Conv. Heating Rate')

    ax = fig.add_subplot(2, 2, 4)
    state['air_temperature'].mean(dim='lon').plot.contourf(
        ax=ax, levels=16)
    ax.set_title('Temperature')

    fig.tight_layout()


# Create plotting object
monitor = PlotFunctionMonitor(plot_function)

climt.set_constants_from_dict({
    'stellar_irradiance': {'value': 200, 'units': 'W m^-2'}})

model_time_step = timedelta(seconds=600)

# Create components
convection = climt.EmanuelConvection()
simple_physics = TimeDifferencingWrapper(climt.SimplePhysics())

radiation = climt.GrayLongwaveRadiation()
slab_surface = climt.SlabSurface()

dycore = climt.GFSDynamicalCore(
    [simple_physics, radiation, slab_surface,
     convection], number_of_damped_levels=5
)
grid = climt.get_grid(nx=128, ny=62)

# Create model state
my_state = climt.get_default_state([dycore], grid_state=grid)

# Set initial/boundary conditions
latitudes = my_state['latitude'].values
longitudes = my_state['longitude'].values

surface_shape = latitudes.shape

# Set initial/boundary conditions
sw_flux_equator = 300
sw_flux_pole = 0

sw_flux_profile = sw_flux_equator - (
    (sw_flux_equator - sw_flux_pole)*(np.sin(np.radians(latitudes))**2))

my_state['downwelling_shortwave_flux_in_air'].values[:] = sw_flux_profile[np.newaxis, :]
my_state['surface_temperature'].values[:] = 290.
my_state['ocean_mixed_layer_thickness'].values[:] = 5

my_state['eastward_wind'].values[:] = np.random.randn(
    *my_state['eastward_wind'].shape)


for i in range(1500*24*6):
    diag, my_state = dycore(my_state, model_time_step)
    my_state.update(diag)
    my_state['time'] += model_time_step

    if i % 6 == 0:
        monitor.store(my_state)

    print(my_state['time'])

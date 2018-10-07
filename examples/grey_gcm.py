import climt
from sympl import (
    PlotFunctionMonitor,
    TimeDifferencingWrapper,
    DataArray
)
import numpy as np
from datetime import timedelta


def plot_function(fig, state):

    ax = fig.add_subplot(2, 2, 1)
    state['specific_humidity'].mean(
        dim='lon').plot.contourf(
            ax=ax, levels=16, robust=True)
    ax.set_title('Specific Humidity')

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

dycore = climt.GFSDynamicalCore(
    [simple_physics, radiation,
     convection], number_of_damped_levels=5
)
grid = climt.get_grid(nx=128, ny=62)

# Create model state
my_state = climt.get_default_state([dycore], grid_state=grid)

# Set initial/boundary conditions
latitudes = my_state['latitude'].values
longitudes = my_state['longitude'].values
surface_shape = latitudes.shape

temperature_equator = 300
temperature_pole = 240

temperature_profile = temperature_equator - (
    (temperature_equator - temperature_pole)*(
        np.sin(np.radians(latitudes))**2))

my_state['surface_temperature'] = DataArray(
    temperature_profile*np.ones(surface_shape),
    dims=['lat', 'lon'], attrs={'units': 'degK'})
my_state['eastward_wind'].values[:] = np.random.randn(
    *my_state['eastward_wind'].shape)

for i in range(1500*24*6):
    diag, my_state = dycore(my_state, model_time_step)
    my_state.update(diag)
    my_state['time'] += model_time_step

    if i % 6 == 0:
        monitor.store(my_state)

    print(my_state['time'])

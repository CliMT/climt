import climt
from sympl import PlotFunctionMonitor
import numpy as np
import matplotlib.pyplot as plt


def plot_function(fig, state):

    fig.set_size_inches(8, 8)
    ax = fig.add_subplot(2, 2, 1)
    state['surface_temperature'].transpose().plot.contourf(ax=ax, levels=16)
    ax.set_title('Surface Temperature')

    ax = fig.add_subplot(2, 2, 2)
    state['convective_heating_rate'].mean(dim='longitude').transpose().plot.contourf(
        ax=ax, levels=16)
    ax.set_title('Convective Heating')

    ax = fig.add_subplot(2, 2, 3)
    state['eastward_wind'].mean(dim='longitude').transpose().plot.contourf(
        ax=ax, levels=16)
    ax.set_title('Zonal Wind')

    ax = fig.add_subplot(2, 2, 4)
    state['air_temperature'].mean(dim='longitude').transpose().plot.contourf(
        ax=ax, levels=16)
    ax.set_title('Air Temperature')

    plt.suptitle('Time: '+str(state['time']))
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])


# Create plotting object
monitor = PlotFunctionMonitor(plot_function)

# Create components
dycore = climt.GFSDynamicalCore(number_of_longitudes=128,
                                number_of_latitudes=64,
                                number_of_damped_levels=5)

model_time_step = dycore._time_step

convection = climt.EmanuelConvection(convective_momentum_transfer_coefficient=1)
simple_physics = climt.SimplePhysics()

simple_physics = simple_physics.prognostic_version()
simple_physics.current_time_step = model_time_step
convection.current_time_step = model_time_step

radiation = climt.GrayLongwaveRadiation()
slab_surface = climt.SlabSurface()

# Create model state
my_state = climt.get_default_state([dycore, radiation, convection, simple_physics, slab_surface],
                                   x=dycore.grid_definition['x'],
                                   y=dycore.grid_definition['y'],
                                   mid_levels=dycore.grid_definition['mid_levels'],
                                   interface_levels=dycore.grid_definition['interface_levels'])

# Set initial/boundary conditions
temperature_equator = 300
temperature_pole = 240
latitudes = my_state['latitude'].values

temperature_profile = temperature_equator - (
    (temperature_equator - temperature_pole)*(np.sin(np.radians(latitudes))**2))

my_state['surface_temperature'].values[:] = temperature_profile[np.newaxis, :]

my_state['eastward_wind'].values[:] = np.random.randn(*my_state['eastward_wind'].shape)
my_state['air_temperature'].values[:] = 240
my_state['ocean_mixed_layer_thickness'].values[:] = 30

dycore.prognostics = [simple_physics, slab_surface, radiation, convection]

for i in range(50000):
    output, diag = dycore(my_state)
    my_state.update(output)
    my_state.update(diag)

    if i % 20 == 0:
        monitor.store(my_state)
        print('max. zonal wind: ', np.amax(my_state['eastward_wind'].values))
        print('max. surf temp: ',
              my_state['surface_temperature'].max(keep_attrs=True))

    my_state['time'] += dycore._time_step

import climt
from sympl import (
    PlotFunctionMonitor, NetCDFMonitor,
    TimeDifferencingWrapper, UpdateFrequencyWrapper,
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


fields_to_store = ['air_temperature', 'air_pressure', 'eastward_wind',
                   'northward_wind', 'air_pressure_on_interface_levels',
                   'surface_pressure', 'upwelling_longwave_flux_in_air',
                   'specific_humidity', 'surface_temperature',
                   'latitude', 'longitude',
                   'convective_heating_rate']
# Create plotting object
monitor = PlotFunctionMonitor(plot_function)
netcdf_monitor = NetCDFMonitor('gcm_without_seasonal_cycle.nc',
                               write_on_store=True,
                               store_names=fields_to_store)

climt.set_constants_from_dict({
    'stellar_irradiance': {'value': 200, 'units': 'W m^-2'}})

model_time_step = timedelta(seconds=600)
# Create components


convection = climt.EmanuelConvection()
simple_physics = TimeDifferencingWrapper(climt.SimplePhysics())

constant_duration = 6

radiation_lw = UpdateFrequencyWrapper(
    climt.RRTMGLongwave(), constant_duration*model_time_step)

radiation_sw = UpdateFrequencyWrapper(
    climt.RRTMGShortwave(), constant_duration*model_time_step)

slab_surface = climt.SlabSurface()

dycore = climt.GFSDynamicalCore(
    [simple_physics, slab_surface, radiation_sw,
     radiation_lw, convection], number_of_damped_levels=5
)
grid = climt.get_grid(nx=128, ny=62)

# Create model state
my_state = climt.get_default_state([dycore], grid_state=grid)

# Set initial/boundary conditions
latitudes = my_state['latitude'].values
longitudes = my_state['longitude'].values

zenith_angle = np.radians(latitudes)
surface_shape = [len(longitudes), len(latitudes)]

my_state['zenith_angle'].values = zenith_angle
my_state['eastward_wind'].values[:] = np.random.randn(
    *my_state['eastward_wind'].shape)
my_state['ocean_mixed_layer_thickness'].values[:] = 50

surf_temp_profile = 290 - (40*np.sin(zenith_angle)**2)
my_state['surface_temperature'].values = surf_temp_profile

for i in range(1500*24*6):
    diag, my_state = dycore(my_state, model_time_step)
    my_state.update(diag)
    my_state['time'] += model_time_step

    if i % (6*24) == 0:
        netcdf_monitor.store(my_state)
        monitor.store(my_state)
        print('max. zonal wind: ',
              np.amax(my_state['eastward_wind'].values))
        print('max. humidity: ',
              np.amax(my_state['specific_humidity'].values))
        print('max. surf temp: ',
              my_state['surface_temperature'].max(keep_attrs=True).values)

    print(my_state['time'])

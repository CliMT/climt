from pyinstrument import Profiler
import sys
import numpy as np
from datetime import timedelta
import types
import inspect

# Import sympl and climt
import sympl
import climt

# ----------------------------------------------------------------------
# Monkey-patching / Decoration Logic
# ----------------------------------------------------------------------

def profile_decorator(func):
    """
    A simple decorator that doesn't do much logic itself but ensures
    the function has a distinct frame if we were using a deterministic profiler,
    or just as a placeholder to say "we decorated this".
    """
    if hasattr(func, '_is_profiled'):
        return func

    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    wrapper._is_profiled = True
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper

def decorate_module_functions(module, predicate=None):
    """
    Recursively (or flatly) decorate functions in a module.
    Only decorates objects that are defined within the 'sympl' or 'climt' packages.
    """
    for name, obj in inspect.getmembers(module):
        # Determine the module of the object
        obj_module = getattr(obj, '__module__', None)

        # Only proceed if the object belongs to sympl or climt
        if obj_module and (obj_module.startswith('sympl') or obj_module.startswith('climt')):

            if inspect.isfunction(obj) or inspect.ismethod(obj):
                if predicate and not predicate(name, obj):
                    continue
                # Try to decorate
                try:
                    setattr(module, name, profile_decorator(obj))
                except (TypeError, AttributeError):
                    pass
            elif inspect.isclass(obj):
                # Decorate methods of the class
                for method_name, method in inspect.getmembers(obj):
                    if inspect.isfunction(method) or inspect.ismethod(method):
                        method_module = getattr(method, '__module__', None)
                        # Check module of the method too, just in case
                        if method_module and (method_module.startswith('sympl') or method_module.startswith('climt')):
                            if predicate and not predicate(method_name, method):
                                continue
                            try:
                                setattr(obj, method_name, profile_decorator(method))
                            except (TypeError, AttributeError):
                                pass

# Apply decorators to core modules as requested
# sympl core
decorate_module_functions(sympl)
# If sympl exposes submodules in its namespace or if we imported them:
import sympl._core.base_components
import sympl._core.util
import sympl._core.init_np_arrays
import sympl._core.restore_dataarray
import sympl._core.dataarray
import sympl._core.state

decorate_module_functions(sympl._core.base_components)
decorate_module_functions(sympl._core.util)
decorate_module_functions(sympl._core.init_np_arrays)
decorate_module_functions(sympl._core.restore_dataarray)
decorate_module_functions(sympl._core.dataarray.DataArray) # Decorate methods of DataArray class
decorate_module_functions(sympl._core.state)

# climt core
import climt._core.initialization
import climt._core.util
import climt._core.constants

decorate_module_functions(climt._core.initialization)
decorate_module_functions(climt._core.util)
decorate_module_functions(climt._core.constants)

# ----------------------------------------------------------------------
# Benchmark Script
# ----------------------------------------------------------------------

from sympl import (
    NetCDFMonitor,
    TimeDifferencingWrapper, UpdateFrequencyWrapper,
    AdamsBashforth
)

fields_to_store = ['air_temperature', 'air_pressure', 'eastward_wind',
                   'northward_wind', 'air_pressure_on_interface_levels',
                   'surface_pressure', 'upwelling_longwave_flux_in_air',
                   'specific_humidity', 'surface_temperature',
                   'latitude', 'longitude',
                   'convective_heating_rate']

# monitor = PlotFunctionMonitor(plot_function)
netcdf_monitor = NetCDFMonitor('gcm_without_seasonal_cycle.nc',
                               write_on_store=True,
                               store_names=fields_to_store)

climt.set_constants_from_dict({
    'stellar_irradiance': {'value': 200, 'units': 'W m^-2'}})

model_time_step = timedelta(seconds=600)

convection = climt.EmanuelConvection()
simple_physics = TimeDifferencingWrapper(climt.SimplePhysics())

constant_duration = 6

radiation_lw = UpdateFrequencyWrapper(
    climt.RRTMGLongwave(), constant_duration*model_time_step)

radiation_sw = UpdateFrequencyWrapper(
    climt.RRTMGShortwave(), constant_duration*model_time_step)

slab_surface = climt.SlabSurface()

dycore = AdamsBashforth(
    [simple_physics, slab_surface, radiation_sw,
     radiation_lw, convection]
)
grid = climt.get_grid(nx=128, ny=62)

my_state = climt.get_default_state([dycore], grid_state=grid)

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

# Fix for solar_cycle_fraction issue
if 'solar_cycle_fraction' in my_state:
    pass

# ----------------------------------------------------------------------
# Profiling Execution
# ----------------------------------------------------------------------

profiler = Profiler(interval=0.001) # 1ms interval
profiler.start()

print("Starting simulation loop for 40 steps...")
for i in range(40):
    diag, my_state = dycore(my_state, model_time_step)
    my_state.update(diag)
    my_state['time'] += model_time_step

    # if i % 5 == 0:
    #     netcdf_monitor.store(my_state)

profiler.stop()

print("Simulation complete. Generating report...")
profiler.write_html('benchmarks/pyinstrument_report.html')
print("Report saved to benchmarks/pyinstrument_report.html")

# Also print text report
print(profiler.output_text(unicode=True, color=True))

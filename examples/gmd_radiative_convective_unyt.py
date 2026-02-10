from sympl import (
    PlotFunctionMonitor, AdamsBashforth, NetCDFMonitor, set_backend
)
from climt import (
    SimplePhysics, get_default_state, UnytBackend, UnytTimeDelta,
    EmanuelConvection, RRTMGShortwave, RRTMGLongwave, SlabSurface
)
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt


def get_data(obj):
    """Helper to get numpy array from either DataArray or UnytStateContainer."""
    if hasattr(obj, 'values'):  # DataArray
        return obj.values
    elif hasattr(obj, 'data'):  # UnytStateContainer
        return np.asarray(obj.data)
    return np.asarray(obj)


def to_units(obj, units):
    """Helper to convert units for either DataArray or UnytStateContainer."""
    if hasattr(obj, 'to_units'):  # DataArray
        return obj.to_units(units)
    elif hasattr(obj, 'data'):  # UnytStateContainer
        # Unyt handles mbar as hPa or mb
        target_units = units.replace('mbar', 'hPa')
        # Sanitization logic similar to UnytBackend
        sanitized_units = target_units.replace("^", "**").replace(" ", "*")
        return obj.data.to(sanitized_units)
    return obj


def plot_function(fig, state):
    ax = fig.add_subplot(2, 2, 1)
    
    heating_rate = to_units(state['air_temperature_tendency_from_convection'], 'degK day^-1')
    pressure = to_units(state['air_pressure'], 'mbar')
    
    ax.plot(get_data(heating_rate).flatten(), get_data(pressure).flatten(), '-o')
    ax.set_title('Conv. heating rate')
    ax.set_xlabel('K/day')
    ax.set_ylabel('millibar')
    ax.grid()

    ax.axes.invert_yaxis()
    ax = fig.add_subplot(2, 2, 2)
    
    temp = state['air_temperature']
    ax.plot(get_data(temp).flatten(), get_data(pressure).flatten(), '-o')
    ax.set_title('Air temperature')
    ax.axes.invert_yaxis()
    ax.set_xlabel('K')
    ax.grid()

    ax = fig.add_subplot(2, 2, 3)
    
    lw_tendency = state['air_temperature_tendency_from_longwave']
    sw_tendency = state['air_temperature_tendency_from_shortwave']
    
    ax.plot(get_data(lw_tendency).flatten(), get_data(pressure).flatten(), '-o', label='LW')
    ax.plot(get_data(sw_tendency).flatten(), get_data(pressure).flatten(), '-o', label='SW')
    
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
    
    pressure_iface = to_units(state['air_pressure_on_interface_levels'], 'mbar')
    
    ax.plot(get_data(net_flux).flatten(), get_data(pressure_iface).flatten(), '-o')
    ax.set_title('Net Flux')
    ax.axes.invert_yaxis()
    ax.set_xlabel('W/m^2')
    ax.grid()
    plt.tight_layout()


def to_dataarray_state(state):
    """Helper to convert a state with UnytStateContainers to DataArrays."""
    from sympl import DataArray
    da_state = {}
    for name, value in state.items():
        if hasattr(value, 'data') and hasattr(value, 'dims'):  # UnytStateContainer
            # Handle numeric data with units
            if hasattr(value.data, 'units'):
                units = str(value.data.units)
                # Normalize units for sympl (pint) if needed
                units = units.replace('degK', 'K')
            else:
                units = ''  # or 'dimensionless' for strings/objects
            
            data = np.asarray(value.data)
            da_state[name] = DataArray(data, dims=value.dims, attrs={'units': units})
        else:
            da_state[name] = value
    return da_state


# Set the backend to Unyt
set_backend(UnytBackend())

monitor = PlotFunctionMonitor(plot_function)

# Use UnytTimeDelta for unit-aware time stepping
timestep = UnytTimeDelta(minutes=5)

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

# Note: NetCDFMonitor might still have issues with UnytBackend if it
# strictly expects xarray objects. Users might need to convert state 
# to numpy or xarray before passing to NetCDFMonitor if it fails.
try:
    netcdf_monitor = NetCDFMonitor('rad_conv_eq_unyt.nc',
                                   store_names=store_quantities,
                                   write_on_store=True)
except Exception as e:
    print(f"Warning: NetCDFMonitor initialization failed: {e}")
    netcdf_monitor = None

convection.current_time_step = timestep

state = get_default_state([simple_physics, convection,
                           radiation_lw, radiation_sw, slab])

def set_values(obj, val):
    if hasattr(obj, 'values'):
        obj.values[:] = val
    elif hasattr(obj, 'data'):
        obj.data[:] = val
    else:
        obj[:] = val

set_values(state['air_temperature'], 270)
set_values(state['surface_albedo_for_direct_shortwave'], 0.5)
set_values(state['surface_albedo_for_direct_near_infrared'], 0.5)
set_values(state['surface_albedo_for_diffuse_shortwave'], 0.5)

set_values(state['zenith_angle'], np.pi/2.5)
set_values(state['surface_temperature'], 300.)
set_values(state['ocean_mixed_layer_thickness'], 5)
set_values(state['area_type'], 'sea')

time_stepper = AdamsBashforth([convection, radiation_lw, radiation_sw, slab])

for i in range(20000):
    convection.current_time_step = timestep
    diagnostics, state = time_stepper(state, timestep)
    state.update(diagnostics)
    diagnostics, new_state = simple_physics(state, timestep)
    state.update(diagnostics)
    if (i+1) % 20 == 0:
        monitor.store(state)
        if netcdf_monitor:
            netcdf_monitor.store(to_dataarray_state(state))
        print(i+1, get_data(state['surface_temperature']).flatten()[0])

    state.update(new_state)
    state['time'] += timestep
    set_values(state['eastward_wind'], 3.)

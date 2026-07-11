import numpy as np
import matplotlib.pyplot as plt
from sympl import PlotFunctionMonitor, AdamsBashforth, set_backend
from climt import (
    SimplePhysics, get_default_state, UnytBackend, UnytTimeDelta,
    EmanuelConvection, SocratesShortwave, SocratesLongwave, SlabSurface
)


def get_data(obj):
    if hasattr(obj, 'values'):
        return obj.values
    elif hasattr(obj, 'data'):
        return np.asarray(obj.data)
    return np.asarray(obj)


def to_units(obj, units):
    if hasattr(obj, 'to_units'):
        return obj.to_units(units)
    elif hasattr(obj, 'data'):
        target_units = units.replace('mbar', 'hPa')
        sanitized_units = target_units.replace("^", "**").replace(" ", "*")
        return obj.data.to(sanitized_units)
    return obj


def plot_function(fig, state):
    ax = fig.add_subplot(2, 2, 1)
    heating_rate = to_units(
        state['air_temperature_tendency_from_convection'], 'degK day^-1'
    )
    pressure = to_units(state['air_pressure'], 'mbar')
    ax.plot(
        get_data(heating_rate).flatten(), get_data(pressure).flatten(), '-o'
    )
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
    ax.plot(
        get_data(lw_tendency).flatten(), get_data(pressure).flatten(),
        '-o', label='LW'
    )
    ax.plot(
        get_data(sw_tendency).flatten(), get_data(pressure).flatten(),
        '-o', label='SW'
    )
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
    pressure_iface = to_units(
        state['air_pressure_on_interface_levels'], 'mbar'
    )
    ax.plot(
        get_data(net_flux).flatten(), get_data(pressure_iface).flatten(), '-o'
    )
    ax.set_title('Net Flux')
    ax.axes.invert_yaxis()
    ax.set_xlabel('W/m^2')
    ax.grid()
    plt.tight_layout()


# Set backend to Unyt
set_backend(UnytBackend())
monitor = PlotFunctionMonitor(plot_function)
timestep = UnytTimeDelta(minutes=5)

convection = EmanuelConvection()
radiation_sw = SocratesShortwave(spectral_file='sp_sw_ga9')
radiation_lw = SocratesLongwave(spectral_file='sp_lw_ga9')
slab = SlabSurface()
simple_physics = SimplePhysics()

convection.current_time_step = timestep
state = get_default_state(
    [simple_physics, convection, radiation_lw, radiation_sw, slab]
)


def set_values(obj, val):
    if hasattr(obj, 'values'):
        obj.values[:] = val
    elif hasattr(obj, 'data'):
        obj.data[:] = val
    else:
        obj[:] = val


# Set initial/boundary values
set_values(state['air_temperature'], 270)
set_values(state['zenith_angle'], np.pi / 2.5)
set_values(state['surface_temperature'], 300.)
set_values(state['ocean_mixed_layer_thickness'], 5)
set_values(state['area_type'], 'sea')

time_stepper = AdamsBashforth([convection, radiation_lw, radiation_sw, slab])

print("Running Radiative Convective Equilibrium (RCE) example...")
# Run time loop
for i in range(2000):
    convection.current_time_step = timestep
    diagnostics, state = time_stepper(state, timestep)
    state.update(diagnostics)
    diagnostics, new_state = simple_physics(state, timestep)
    state.update(diagnostics)
    if (i + 1) % 100 == 0:
        monitor.store(state)
        temp_val = get_data(state['surface_temperature']).flatten()[0]
        print(f"Step {i+1:04d}: Surface Temp = {temp_val:.2f} K")
    state.update(new_state)
    state['time'] += timestep
    set_values(state['eastward_wind'], 3.)

print("Simulation finished. Close the plot window to exit.")
plt.show()

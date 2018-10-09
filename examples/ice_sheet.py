import climt
from datetime import timedelta
from sympl import PlotFunctionMonitor


ice = climt.IceSheet(maximum_snow_ice_height=30.)

state = climt.get_default_state([ice])
state['area_type'].values[:] = 'sea_ice'
state['sea_ice_thickness'].values[:] = 5.
state['surface_snow_thickness'].values[:] = 4.
state['surface_temperature'].values[:] = 260.
state['surface_upward_sensible_heat_flux'].values[:] = -0.5


def plot_function(fig, model_state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(model_state['snow_and_ice_temperature'].values.squeeze(),
            model_state['height_on_ice_interface_levels'].values.squeeze(), 'o-')
    ax.set_title('Evolution of Temperature in Ice-Snow layers')
    ice_thickness = state['sea_ice_thickness'].values.item()
    ax.axhline(ice_thickness, linestyle='--', label='ice-snow interface')
    ax.set_ylabel('model levels')
    ax.set_xlabel('Kelvin')
    ax.set_xlim(270, 273.5)
    ax.legend()


monitor = PlotFunctionMonitor(plot_function)

for i in range(3000):
    diag, out = ice(state, timedelta(minutes=100))
    print('Surface Temperature: ', out['surface_temperature'].values.item())
    print('Sea Ice thickness: ', state['sea_ice_thickness'].values.item())
    print('--------------------------')
    state.update(diag)
    state.update(out)
    monitor.store(state)

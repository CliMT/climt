import climt
from datetime import timedelta
from sympl import PlotFunctionMonitor


ice = climt.IceSheet(maximum_snow_ice_height=30., number_vertical_levels=30)

state = climt.get_default_state([ice])
state['area_type'].values[:] = 'sea_ice'
state['sea_ice_thickness'].values[:] = 5.
state['surface_snow_thickness'].values[:] = 4.
state['surface_temperature'].values[:] = 260.
state['surface_upward_sensible_heat_flux'].values[:] = -0.5

num_levels = ice._output_levels


def plot_function(fig, model_state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(model_state['snow_and_ice_temperature'].values.squeeze(),
            model_state['ice_vertical_levels'].values, 'o-')
    ax.set_title('Evolution of Temperature on Ice-Snow layers')
    ice_thickness = state['sea_ice_thickness'].values.item()
    snow_thickness = state['surface_snow_thickness'].values.item()
    ice_level = ice_thickness/(ice_thickness + snow_thickness)*num_levels
    ax.axhline(ice_level, linestyle='--', label='ice-snow interface')
    ax.set_ylabel('model levels')
    ax.set_xlabel('Kelvin')
    ax.legend()


monitor = PlotFunctionMonitor(plot_function)

for i in range(3000):
    out, diag = ice(state, timedelta(minutes=100))
    print('Surface Temperature: ', out['surface_temperature'].values.item())
    print('Sea Ice thickness: ', state['sea_ice_thickness'].values.item())
    print('--------------------------')
    state.update(diag)
    state.update(out)
    monitor.store(state)

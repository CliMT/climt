import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from climt import RRTMGShortwave, RRTMGLongwave, get_default_state
from sympl import get_constant

prop_cycle = rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

rad_sw = RRTMGShortwave(mcica=True)
state_sw = get_default_state([rad_sw])

rad_lw = RRTMGLongwave(mcica=True)
state_lw = get_default_state([rad_lw])

p = state_lw['air_pressure'][:]
p_interface = state_lw['air_pressure_on_interface_levels'][:]
T = state_lw['air_temperature'][:]
R = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
g = get_constant('gravitational_acceleration', 'm s^-2')
density = p / (R * T)
dz = - np.diff(p_interface, axis=0) / (density * g)  # [m]
z = np.cumsum(dz) * 10**-3  # [km]
ice_density = 0.5 * 10**-3  # [kg m^-3]
cloud_base = 10  # [km]
cloud_top = 15  # [km]
cloud_loc = np.where((z > cloud_base) & (z < cloud_top))
i = 0
for area_fraction in np.arange(0, 1.1, 0.25):
    mass_ice_array = area_fraction * ice_density * dz
    for state in state_sw, state_lw:
        state['mass_content_of_cloud_ice_in_atmosphere_layer'][:] = mass_ice_array
        state['cloud_area_fraction_in_atmosphere_layer'][cloud_loc] = area_fraction
    sw_heating = rad_sw(state_sw)[1]['air_temperature_tendency_from_shortwave']
    lw_output = rad_lw(state_lw)
    lw_heating = lw_output[1]['air_temperature_tendency_from_longwave']
    if i == 0:
        label_sw = f'shortwave heating, area fraction {int(area_fraction*100)} %'
        label_lw = 'longwave cooling'
    else:
        label_sw = f'area fraction: {int(area_fraction*100)} %'
        label_lw = ''
    plt.plot(sw_heating.squeeze(), z.squeeze(), label=label_sw, c=colors[i])
    plt.plot(lw_heating.squeeze(), z.squeeze(), ls='--', label=label_lw, c=colors[i])
    i += 1
plt.axhspan(cloud_base, cloud_top, color='gray', alpha=0.5, label='cloud location')
plt.ylabel('Altitude [km]')
plt.xlabel('Heating rate [K/day]')
plt.legend()
plt.show()

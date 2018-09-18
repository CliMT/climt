import numpy as np
import matplotlib.pyplot as plt
from climt import RRTMGShortwave, get_default_state
from sympl import get_constant

rad_sw = RRTMGShortwave(mcica=True)
state = get_default_state([rad_sw])
p = state['air_pressure'][:]
p_interface = state['air_pressure_on_interface_levels'][:]
T = state['air_temperature'][:]
R = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
g = get_constant('gravitational_acceleration', 'm s^-2')
density = p / (R * T)
dz = - np.diff(p_interface) / (density * g)  # [m]
z = np.cumsum(dz) * 10**-3  # [km]
ice_density = 0.5 * 10**-3  # [kg m^-3]
cloud_base = 10  # [km]
cloud_top = 15  # [km]
cloud_loc = np.where((z > cloud_base) & (z < cloud_top))
for area_fraction in np.arange(0, 1.1, 0.25):
    mass_ice_array = area_fraction * ice_density * dz
    state['mass_content_of_cloud_ice_in_atmosphere_layer'][:] = mass_ice_array
    state['cloud_area_fraction_in_atmosphere_layer'][cloud_loc] = area_fraction
    sw_heating = rad_sw(state)[1]['air_temperature_tendency_from_shortwave']
    plt.plot(sw_heating, z, label=f'area fraction: {int(area_fraction*100)} %')
    plt.fill_betweenx(z[cloud_loc], 0, sw_heating[cloud_loc], alpha=0.4)
plt.ylabel('Altitude [km]')
plt.xlabel('Shortwave heating rate [K/day]')
plt.legend()
plt.show()

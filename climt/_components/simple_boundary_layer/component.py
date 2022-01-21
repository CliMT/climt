from sympl import initialize_numpy_arrays_with_properties, get_constant
from sympl import Stepper
import numpy as np
import numba
from numba import jit


@jit(nopython=True, parallel=True)
def Parallel_columns(air_temperature, surface_temperature, air_pressure,
                     air_pressure_int, surface_pressure, specific_humidity,
                     surface_humidity, northward_wind, eastward_wind,
                     new_air_temperature, new_specific_humidity,
                     new_northward_wind, new_eastward_wind,
                     north_wind_stress, east_wind_stress, boundary_height,
                     Rd_val, Cp_val, g_val, k_val, z0_val, fb_val, P0_val,
                     Ric_val, num_cols, timestep):

    Rd = Rd_val
    Cp_dry = Cp_val
    g = g_val
    k = k_val
    z0 = z0_val
    fb = fb_val
    P0 = P0_val
    Ri_c = Ric_val

    def K_b(Ri, Ri_a, u_fric, C, z):

        if Ri_a <= 0:
            Kb = k*u_fric*np.sqrt(C)*z
        else:
            Kb = k*u_fric*np.sqrt(C)*z/(1+Ri/Ri_c*np.log(z/z0)/(1-Ri/Ri_c))

        return(Kb)

    def TDMAsolver(a, b, c, d):

        n = len(d)
        w = np.zeros(n-1)
        g = np.zeros(n)
        p = np.zeros(n)

        w[0] = c[0]/b[0]
        g[0] = d[0]/b[0]

        for i in range(1, n-1):
            w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
        for i in range(1, n):
            g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
        p[n-1] = g[n-1]
        for i in range(n-1, 0, -1):
            p[i-1] = g[i-1] - w[i-1]*p[i]

        return p

    def diffuse_profile(profile, p, p_int, rho, Diff, dt):

        num_layers = len(profile)

        diag_m = np.zeros(num_layers)
        diag_p = np.zeros(num_layers)
        diag = np.zeros(num_layers)

        for i in range(num_layers):

            if i != 0:
                diag_m[i] = g*g*rho[i-1]*rho[i-1]*Diff[i-1]*dt/(p[i-1]-p[i])\
                    * 1/(p_int[i]-p_int[i+1])

            if i != num_layers-1:
                diag_p[i] = g*g*rho[i]*rho[i]*Diff[i]*dt/(p[i]-p[i+1]) * \
                    1/(p_int[i]-p_int[i+1])

            diag[i] = 1+diag_m[i]+diag_p[i]

        return TDMAsolver(-diag_m[1:], diag, -diag_p[:-1], profile)

    for col in numba.prange(num_cols):

        col_air_temperature = air_temperature[:, col]
        col_surface_temperature = surface_temperature[col]
        col_air_pressure = air_pressure[:, col]
        col_air_pressure_int = air_pressure_int[:, col]
        col_surface_pressure = surface_pressure[col]
        col_specific_humidity = specific_humidity[:, col]
        col_surface_humidity = surface_humidity[col]
        col_north_wind = northward_wind[:, col]
        col_east_wind = eastward_wind[:, col]

        col_north_wind_int = 0.5*(col_north_wind[1:]+col_north_wind[:-1])
        col_east_wind_int = 0.5*(col_east_wind[1:]+col_east_wind[:-1])
        col_air_temperature_int = 0.5*(col_air_temperature[1:] +
                                       col_air_temperature[:-1])
        col_specific_humidity_int = 0.5*(col_specific_humidity[1:] +
                                         col_specific_humidity[:-1])
        col_rho = col_air_pressure_int[1:-1]/(Rd * (1+0.608 *
                                              col_specific_humidity_int) *
                                              col_air_temperature_int)

        n = len(col_air_temperature_int)

        pot_virt_temp = col_air_temperature_int *\
            np.power((P0/col_air_pressure_int[1:-1]), Rd/Cp_dry) * \
            (1+0.608*col_specific_humidity_int)
        pot_virt_temp_surf = col_surface_temperature *\
            np.power((P0/col_surface_pressure), Rd/Cp_dry) *\
            (1+0.608*col_surface_humidity)

        z = np.zeros(n)
        z[0] = (Rd*(1+0.608*col_specific_humidity[0])*col_air_temperature[0] /
                g) * np.log(col_surface_pressure/col_air_pressure_int[1:-1][0])
        for i in range(1, n):
            z[i] = z[i-1]+(Rd*(1+0.608*col_specific_humidity[i]) *
                           col_air_temperature[i]/g) *\
                np.log(col_air_pressure_int[1:-1][i-1] /
                       col_air_pressure_int[1:-1][i])

        wind_int = np.sqrt(np.power(col_north_wind_int, 2) +
                           np.power(col_east_wind_int, 2))
        for i in range(len(wind_int)):
            if wind_int[i] < 1:
                wind_int[i] = 1

        Ri_a = g*z[0]*(pot_virt_temp[0]-pot_virt_temp_surf)/(
            pot_virt_temp_surf*wind_int[0]*wind_int[0])
        if Ri_a < 0:
            C = k*k*np.power(np.log(z[0]/z0), -2)
        elif Ri_a < Ri_c:
            C = k*k*np.power(np.log(z[0]/z0), -2)*np.power((1-Ri_a/Ri_c), 2)
        else:
            C = 0

        count = 0
        Rich = np.zeros(n)
        for i in range(n):
            Rich[i] = g*z[i]*(pot_virt_temp[i]-pot_virt_temp[0])/(
                pot_virt_temp[0]*wind_int[i]*wind_int[i])
            if Rich[i] > Ri_c:
                count = i+1
                break
        h = z[count-1]
        boundary_height[col] = h

        north_wind_stress[col] = col_rho[0]*C*wind_int[0]*col_north_wind_int[0]
        east_wind_stress[col] = col_rho[0]*C*wind_int[0]*col_east_wind_int[0]

        u_fric = wind_int[0]

        diff = np.zeros(n)

        for i in range(count):
            if z[i] < fb*h:
                diff[i] = K_b(Rich[i], Ri_a, u_fric, C, z[i])

            else:
                diff[i] = K_b(Rich[i], Ri_a, u_fric, C, fb*h)*z[i]/(h*fb) *\
                    np.power(1-(z[i]-fb*h)/((1-fb)*h), 2)

        new_temp = diffuse_profile(col_air_temperature, col_air_pressure,
                                   col_air_pressure_int, col_rho, diff,
                                   timestep)
        new_humidity = diffuse_profile(col_specific_humidity, col_air_pressure,
                                       col_air_pressure_int, col_rho, diff,
                                       timestep)
        new_north_wind = diffuse_profile(col_north_wind, col_air_pressure,
                                         col_air_pressure_int, col_rho, diff,
                                         timestep)
        new_east_wind = diffuse_profile(col_east_wind, col_air_pressure,
                                        col_air_pressure_int, col_rho, diff,
                                        timestep)

        new_air_temperature[:, col] = new_temp
        new_specific_humidity[:, col] = new_humidity
        new_northward_wind[:, col] = new_north_wind
        new_eastward_wind[:, col] = new_east_wind


class SimpleBoundaryLayer(Stepper):
    """
    This is a simple boundary layer component that diffuses heat, humidity and
    momemtum upwards from the lowest model level.

    This component assumes that a surface flux component has been already run,
    which has made the changes due to surface fluxes at the lowest model
    level. This component then diffuses heat, humidity and momentum using
    diffusion coefficients calculated using the simplified Monin-Obukhov
    theory.
    """

    input_properties = {
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK ',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'kg/kg',
        },
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa',
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'Pa',
        },
        'northward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'eastward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'surface_air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'surface_specific_humidity': {
            'dims': ['*'],
            'units': 'kg/kg',
        },
    }

    output_properties = {
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK ',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'kg/kg',
        },
        'northward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'eastward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
    }

    diagnostic_properties = {
        'northward_wind_stress': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'eastward_wind_stress': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'boundary_layer_height': {
            'dims': ['*'],
            'units': 'm',
        },
    }

    def __init__(self, von_karman_constant=0.4, roughness_length=0.0000321,
                 specific_fraction=0.1, reference_pressure=100000,
                 critical_richardson_number=1, **kwargs):
        """
        Args:
        roughness_length:
            A measure of the surface roughness.
        specific_fraction:
            A parameter used in the calculation of diffusion coefficients.
        reference_pressure:
            The reference pressure used in the potential temperature
            calculations.
        critical_richardson_number:
            A set threshold value which determines the diffusion coefficients
            and the height of the boundary layer.
        """

        self._k = von_karman_constant
        self._z0 = roughness_length
        self._fb = specific_fraction
        self._P0 = reference_pressure
        self._Ric = critical_richardson_number
        self._update_constants()

        super(SimpleBoundaryLayer, self).__init__(**kwargs)

    def _update_constants(self):

        self._Rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        self._Cp =\
            get_constant('heat_capacity_of_dry_air_at_constant_pressure',
                         'J kg^-1 K^-1')
        self._g = get_constant('gravitational_acceleration', 'm s^-2')

    def array_call(self, state, timestep):
        """
        Takes temperature, humidty and wind profiles for each column and
        returns diffused temperature, humidity and wind profiles.
        """

        num_cols = state['air_temperature'].shape[1]

        new_state = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties
        )

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        Parallel_columns(state['air_temperature'],
                         state['surface_temperature'],
                         state['air_pressure'],
                         state['air_pressure_on_interface_levels'],
                         state['surface_air_pressure'],
                         state['specific_humidity'],
                         state['surface_specific_humidity'],
                         state['northward_wind'], state['eastward_wind'],
                         new_state['air_temperature'],
                         new_state['specific_humidity'],
                         new_state['northward_wind'],
                         new_state['eastward_wind'],
                         diagnostics['northward_wind_stress'],
                         diagnostics['eastward_wind_stress'],
                         diagnostics['boundary_layer_height'],
                         self._Rd, self._Cp, self._g, self._k, self._z0,
                         self._fb, self._P0, self._Ric, num_cols,
                         timestep.total_seconds())

        return diagnostics, new_state

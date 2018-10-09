from sympl import TendencyComponent, get_constant
import numpy as np


class HeldSuarez(TendencyComponent):
    """
    Provide the Held-Suarez forcing.

    Produces the forcings proposed by Held and Suarez for the intercomparison
    of dynamical cores of AGCMs. Relaxes the temperature field to a zonally
    symmetric equilibrium state, and uses Rayleigh damping of low-level winds
    to represent boundary-layer friction. Details can be found in
    `[Held and Suarez (1994)]`_.

    References:
        Held, I. and M. Suarez, 1994:
            A Proposal for the Intercomparison of the
            Dynamical Cores of Atmospheric General Circulation Models.
            Bull. Amer. Meteor. Soc., 75, 1825-1830,
            doi: 10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2.

    .. _[Held and Suarez (1994)]:
        http://journals.ametsoc.org/doi/pdf/10.1175/1520-0477(1994)075%3C1825%3AAPFTIO%3E2.0.CO%3B2
    """

    input_properties = {
        'eastward_wind': {
            'dims': ['*', 'mid_levels'],
            'units': 'm s^-1',
        },
        'northward_wind': {
            'dims': ['*', 'mid_levels'],
            'units': 'm s^-1',
        },
        'air_temperature': {
            'dims': ['*', 'mid_levels'],
            'units': 'degK',
        },
        'air_pressure': {
            'dims': ['*', 'mid_levels'],
            'units': 'Pa',
        },
        'surface_air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'latitude': {
            'dims': ['*'],
            'units': 'degrees_north',
        }
    }

    tendency_properties = {
        'eastward_wind': {'units': 'm s^-2'},
        'northward_wind': {'units': 'm s^-2'},
        'air_temperature': {'units': 'degK s^-1'},
    }

    diagnostic_properties = {}

    def __init__(self,
                 sigma_boundary_layer_top=0.7,
                 k_f=1/86400.,
                 k_a=1/40./86400.,
                 k_s=1/4./86400.,
                 equator_pole_temperature_difference=60,
                 delta_theta_z=10,
                 **kwargs):
        """

        Args:
            sigma_boundary_layer_top (float): The height of the boundary
                layer top in sigma coordinates. Corresponds to $\sigma_b$
                in Held and Suarez, 1994. Default is 0.7.

            k_f (float):
                Velocity damping coefficient at the surface in :math:`s^{-1}`.
                Default is :math:`1\ day^{-1}`.

            k_a (float):
                Parameter used in defining vertical profile of the
                temperature damping in :math:`s^{-1}`, as outlined in
                Held and Suarez, 1994.
                Default is :math:`1/40\ day^{-1}`.

            k_s (float):
                Parameter used in defining vertical profile of the
                temperature damping in :math:`s^{-1}`, as outlined in
                Held and Suarez, 1994.
                Default is :math:`1/4\ day^{-1}`.

            equator_pole_temperature_difference (float):
                Equator to pole
                temperature difference, in K.
                Default is 60K.

            delta_theta_z (float):
                Parameter used in defining the equilibrium
                temperature profile as outlined in Held and Suarez, 1994, in K.
                Default is 10K.

        """

        self._sigma_b = sigma_boundary_layer_top
        self._k_f = k_f
        self._k_a = k_a
        self._k_s = k_s
        self._delta_T_y = equator_pole_temperature_difference
        self._delta_theta_z = delta_theta_z

        self._update_constants()

        super(HeldSuarez, self).__init__(**kwargs)

    def _update_constants(self):
        self._p0 = get_constant('reference_air_pressure', 'Pa')
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
        self._R_d = get_constant('gas_constant_of_dry_air', 'J/kg/degK')
        self._kappa = self._R_d/self._Cpd
        self._Omega = get_constant('planetary_rotation_rate', 's^-1')
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._r_planet = get_constant('planetary_radius', 'm')

    def array_call(self, raw_state):
        """
        Get the Held-Suarez forcing tendencies

        Args:
            raw_state (dict): A model state dictionary of numpy arrays
                satisfying input_properties.

        Returns:
            tendencies (dict), diagnostics (dict):

            * A dictionary whose keys are strings indicating
              state quantities and values are the time derivative of those
              quantities at the time of the input state.
            * A dictionary whose keys are strings indicating
              state quantities and values are the value of those quantities
              at the time of the input state.
        """
        self._update_constants()
        sigma = raw_state['air_pressure'] / raw_state['surface_air_pressure'][:, np.newaxis]

        Teq = self._get_Teq(raw_state['latitude'][:, np.newaxis], raw_state['air_pressure'])
        k_t = self._get_k_t(raw_state['latitude'][:, np.newaxis], sigma)
        k_v = self._get_k_v(sigma)

        tendencies = {
            'eastward_wind': - k_v * raw_state['eastward_wind'],
            'northward_wind': - k_v * raw_state['northward_wind'],
            'air_temperature': - k_t * (raw_state['air_temperature'] - Teq),
        }

        return tendencies, {}

    def _get_Teq(self, latitude, air_pressure):
        return np.maximum(
            200,
            (315 - self._delta_T_y*np.sin(np.radians(latitude))**2 -
             self._delta_theta_z*np.log(air_pressure/self._p0)*np.cos(np.radians(latitude))**2
             ) * (air_pressure/self._p0)**self._kappa
        )

    def _get_k_t(self, latitude, sigma):
        return (
            self._k_a +
            (self._k_s - self._k_a) *
            np.maximum(0, (sigma - self._sigma_b)/(1 - self._sigma_b)) *
            np.cos(np.radians(latitude))**4
        )

    def _get_k_v(self, sigma):
        return self._k_f * np.maximum(0, (sigma - self._sigma_b)/(1 - self._sigma_b))

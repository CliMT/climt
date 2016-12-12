from .._core.base_components import Prognostic
from .._core.util import replace_none_with_default, ensure_3d
from .._core.array import DataArray
import numpy as np


class HeldSuarez(Prognostic):

    def __init__(self, latitude=None, air_pressure=None, sigma=None,
                 sigma_b=0.7,
                 k_f=1/86400., k_a=1/40./86400., k_s=1/4./86400.,
                 delta_T_y=60, delta_theta_z=10,
                 reference_pressure=None, gas_constant_for_dry_air=None,
                 heat_capacity_of_dry_air_at_constant_pressure=None,
                 planetary_rotation_rate=None, gravitational_acceleration=None,
                 planetary_radius=None):
        self._sigma_b = sigma_b
        self._k_f = k_f
        self._k_a = k_a
        self._k_s = k_s
        self._delta_T_y = delta_T_y
        self._delta_theta_z = delta_theta_z
        self._p0 = replace_none_with_default(
            'reference_pressure', reference_pressure)
        self._Cpd = replace_none_with_default(
            'heat_capacity_of_dry_air_at_constant_pressure',
            heat_capacity_of_dry_air_at_constant_pressure)
        self._R_d = replace_none_with_default(
            'gas_constant_for_dry_air', gas_constant_for_dry_air)
        self._kappa = self._R_d/self._Cpd
        self._Omega = replace_none_with_default(
            'planetary_rotation_rate', planetary_rotation_rate)
        self._g = replace_none_with_default(
            'gravitational_acceleration', gravitational_acceleration)
        self._r_planet = replace_none_with_default(
            'planetary_radius', planetary_radius)
        # cache computed profiles if grid coordinates are given
        if air_pressure is not None and latitude is not None:
            self._Teq = self._get_Teq(latitude, air_pressure)
        else:
            self._Teq = None
        if latitude is not None and sigma is not None:
            self._k_t = self._get_k_t(latitude, sigma)
        else:
            self._k_t = None
        if sigma is not None:
            self._k_v = self._get_k_v(sigma)
        else:
            self._k_v = None

    def __call__(self, state):
        """
        Gets tendencies and diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            tendencies (dict): A dictionary whose keys are strings indicating
                state quantities and values are the time derivative of those
                quantities in units/second at the time of the input state.
            diagnostics (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.
        """
        if self._Teq is None:
            Teq = self._get_Teq(state['latitude'], state['air_pressure'])
        else:
            Teq = self._Teq
        if self._k_t is None:
            k_t = self._get_k_t(state['latitude'], state['sigma'])
        else:
            k_t = self._k_t
        if self._k_v is None:
            k_v = self._get_k_v(state['sigma'])
        else:
            k_v = self._k_v
        tendencies = {
            'eastward_wind': DataArray(
                - k_v * state['eastward_wind'].to_units('m s^-1').values,
                dims=state['eastward_wind'].dims, attrs={'units': 'm s^-2'}),
            'northward_wind': DataArray(
                - k_v * state['northward_wind'].to_units('m s^-1').values,
                dims=state['northward_wind'].dims, attrs={'units': 'm s^-2'}),
            'air_temperature': DataArray(
                - k_t * (state['air_temperature'].to_units('K').values - Teq),
                dims=state['air_temperature'].dims, attrs={'units': 'K s^-1'})
        }
        return tendencies, {}

    def _get_Teq(self, latitude, air_pressure):
        latitude = ensure_3d(
            latitude.to_units('degrees_N').values, data_dim_if_1d=2)
        air_pressure = ensure_3d(
            air_pressure.to_units('Pa').values, data_dim_if_1d=3)
        return np.maximum(
            200,
            (315 - self._delta_T_y*np.sin(latitude)**2 -
             self._delta_theta_z*np.log(air_pressure/self._p0)*np.cos(latitude)**2
             ) * (air_pressure/self._p0)**self._kappa
        )

    def _get_k_t(self, latitude, sigma):
        latitude = ensure_3d(
            latitude.to_units('degrees_north').values, data_dim_if_1d=2)
        sigma = ensure_3d(sigma.to_units('').values, data_dim_if_1d=3)
        return (
            self._k_a +
            (self._k_s - self._k_a) *
            np.maximum(0, (sigma - self._sigma_b)/(1 - self._sigma_b)) *
            np.cos(latitude)**4)

    def _get_k_v(self, sigma):
        sigma = ensure_3d(sigma.to_units('').values, data_dim_if_1d=3)
        return self._k_f * np.maximum(
            0,
            (sigma - self._sigma_b)/(1 - self._sigma_b))

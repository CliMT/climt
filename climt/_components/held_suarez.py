from .._core.base_components import Prognostic
from .._core.util import replace_none_with_default, ensure_3d
from .._core.array import DataArray
import numpy as np


class HeldSuarez(Prognostic):
    """
    Produces the forcings proposed by Held and Suarez for the intercomparison
    of dynamical cores of AGCMs. Relaxes the temperature field to a zonally
    symmetric equilibrium state, and uses Rayleigh damping of low-level winds
    to represent boundary-layer friction. Details can be found in
    Held and Suarez (1994).

    Attributes:
        inputs (tuple of str): The quantities required in the state when the
            object is called.
        tendencies (tuple of str): The quantities for which tendencies are
            returned when the object is called.
        diagnostics (tuple of str): The diagnostic quantities returned when
            the object is called.

    References:
        Held, I. and M. Suarez, 1994: A Proposal for the Intercomparison of the
            Dynamical Cores of Atmospheric General Circulation Models.
            Bull. Amer. Meteor. Soc., 75, 1825–1830,
            doi: 10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2.
    """

    inputs = None  # defined at initialization
    tendencies = ('eastward_wind', 'northward_wind', 'air_temperature')
    diagnostics = ()

    def __init__(self, latitude=None, air_pressure=None, sigma=None,
                 sigma_boundary_layer_top=0.7,
                 k_f=1/86400., k_a=1/40./86400., k_s=1/4./86400.,
                 equator_pole_temperature_difference=60, delta_theta_z=10,
                 reference_pressure=None, gas_constant_for_dry_air=None,
                 heat_capacity_of_dry_air_at_constant_pressure=None,
                 planetary_rotation_rate=None, gravitational_acceleration=None,
                 planetary_radius=None):
        """

        Args:
            latitude (DataArray): The latitude coordinate. Passing this in at
                initialization makes it so that this object uses that value
                instead of any value that might be present in the state.
            air_pressure (DataArray): The air pressure coordinate. Passing this in at
                initialization makes it so that this object uses that value
                instead of any value that might be present in the state. Do not
                pass this in if pressure is not a constant coordinate.
            sigma (DataArray): The sigma (air_pressure/surface_pressure) coordinate.
                Passing this in at
                initialization makes it so that this object uses that value
                instead of any value that might be present in the state.
                Do not pass this in if sigma is not a constant coordinate.
            sigma_boundary_layer_top (float): The height of the boundary
                layer top in sigma coordinates. Corresponds to $\sigma_b$
                in Held and Suarez, 1994. Default is 0.7.
            k_f (float): Velocity damping coefficient at the surface in s^{-1}.
                Default is $1 day^{-1}$.
            k_a (float): Parameter used in defining vertical profile of the
                temperature damping in $s^{-1}$, as outlined in
                Held and Suarez, 1994.
                Default is $1/40 day^{-1}$
            k_s (float): Parameter used in defining vertical profile of the
                temperature damping in $s^{-1}$, as outlined in
                Held and Suarez, 1994.
                Default is $1/40 day^{-1}$
            equator_pole_temperature_difference (float): Equator to pole
                temperature difference, in K.
                Default is 60K.
            delta_theta_z (float): Parameter used in defining the equilibrium
                temperature profile as outlined in Held and Suarez, 1994, in K.
                Default is 10K.
            reference_pressure (float): Parameter used to define the
                equilibrium temperature profile, roughly equal to the surface
                pressure, in Pa. Default value is $10^5$ Pa.
            gas_constant_for_dry_air (float): Value in $J kg^{-1} K^{-1}$.
                Default is taken from climt.default_constants.
            heat_capacity_of_dry_air_at_constant_pressure: Value in
                $J kg^{-1} K^{-1}$.
                Default is taken from climt.default_constants.
            planetary_rotation_rate (float) Value in $s^{-1}$.
                Default is taken from climt.default_constants.
            gravitational_acceleration (float): Value in $m s^{-2}$.
                Default is taken from climt.default_constants.
            planetary_radius (float): Value in m.
                Default is taken from climt.default_constants.
        """
        self._sigma_b = sigma_boundary_layer_top
        self._k_f = k_f
        self._k_a = k_a
        self._k_s = k_s
        self._delta_T_y = equator_pole_temperature_difference
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
        self.inputs = ['eastward_wind', 'northward_wind', 'air_temperature']
        if sigma is None:
            self.inputs.append('sigma')
        if latitude is None:
            self.inputs.append('latitude')
        if air_pressure is None:
            self.inputs.append('air_pressure')
        self.inputs = tuple(self.inputs)

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

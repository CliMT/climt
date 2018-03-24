from sympl import (
    Diagnostic, DataArray, jit,
    combine_dimensions, get_numpy_array)
from .._core import ClimtPrognostic, get_constant
import numpy as np


class GrayLongwaveRadiation(ClimtPrognostic):

    _climt_inputs = {
        'longwave_optical_depth_on_interface_levels': 'dimensionless',
        'air_temperature': 'degK',
        'surface_temperature': 'degK',
        'air_pressure': 'Pa',
        'air_pressure_on_interface_levels': 'Pa'}

    _climt_diagnostics = {
        'downwelling_longwave_flux_in_air': 'W m^-2',
        'upwelling_longwave_flux_in_air': 'W m^-2',
        'longwave_heating_rate': 'degK/day'}

    _climt_tendencies = {'air_temperature': 'degK/s'}

    def __init__(self):
        """

        Initialise component.

        """

        self._stefan_boltzmann = get_constant('stefan_boltzmann_constant', 'W/m^2/K^4')
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K')

    def __call__(self, state):
        """
        Get heating tendencies and longwave fluxes.

        Args:

            state (dict): A model state dictionary.

        Returns:

            tendencies (dict), diagnostics (dict):

                * A dictionary whose keys are strings indicating
                  state quantities and values are the time derivative of those
                  quantities in units/second at the time of the input state.
                * A dictionary whose keys are strings indicating
                  state quantities and values are the value of those quantities
                  at the time of the input state.
        """
        tau = get_numpy_array(state[
            'longwave_optical_depth_on_interface_levels'].to_units(''),
            out_dims=('x', 'y', 'z'))
        T = get_numpy_array(
            state['air_temperature'].to_units('degK'), out_dims=('x', 'y', 'z'))
        p_interface = get_numpy_array(
            state['air_pressure_on_interface_levels'].to_units('Pa'),
            out_dims=('x', 'y', 'z'))
        Ts = get_numpy_array(
            state['surface_temperature'].to_units('degK'), out_dims=('x', 'y'))

        (downward_flux, upward_flux, net_lw_flux,
         lw_temperature_tendency, tau) = get_longwave_fluxes(
            T, p_interface, Ts, tau, self._stefan_boltzmann,
            self._g, self._Cpd)

        tendencies = self.create_state_dict_for('_climt_tendencies', state)
        tendencies['air_temperature'].values = lw_temperature_tendency

        diagnostics = self.create_state_dict_for('_climt_diagnostics', state)

        diagnostics['downwelling_longwave_flux_in_air'].values = downward_flux
        diagnostics['upwelling_longwave_flux_in_air'].values = upward_flux
        diagnostics['longwave_heating_rate'].values = \
            tendencies['air_temperature'].to_units('degK/day')

        return tendencies, diagnostics


class Frierson06LongwaveOpticalDepth(Diagnostic):

    inputs = ('sigma_on_interface_levels',)
    diagnostics = ('longwave_optical_depth_on_interface_levels',)

    def __init__(
            self,
            linear_optical_depth_parameter=0.1,
            longwave_optical_depth_at_equator=6,
            longwave_optical_depth_at_poles=1.5):
        """
        Args:

            linear_optical_depth_parameter (float, optional): The constant :math:`f_l` which
                determines how much of the variation of :math:`\\tau` with pressure
                is linear rather than quartic.
                :math:`\\tau = \\tau_0 [f_l \\frac{p}{p_s} + (1 - f_l) (\\frac{p}{p_s})^4]`
                Default is 0.1 as in `[Frierson et al., 2006]`_.

            longwave_optical_depth_at_equator (float, optional): The value of :math:`\\tau_0`
                at the equator.
                Default is 6 as in `[Frierson et al., 2006]`_.

            longwave_optical_depth_at_poles (float, optional): The value of :math:`\\tau_0`
                at the poles.
                Default is 1.5 as in `[Frierson et al., 2006]`_.

        .. _[Frierson et al., 2006]:
            http://journals.ametsoc.org/doi/abs/10.1175/JAS3753.1

        """
        self._fl = linear_optical_depth_parameter
        self._tau0e = longwave_optical_depth_at_equator
        self._tau0p = longwave_optical_depth_at_poles

    def __call__(self, state):
        """
        Calculate longwave optical depth from input state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            tendencies (dict): A dictionary whose keys are strings indicating
                state quantities and values are the time derivative of those
                quantities in units/second at the time of the input state.
            diagnostics (dict): A dictionary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.
        """
        lat = get_numpy_array(
            state['latitude'].to_units('degrees_north'),
            out_dims=('x', 'y', 'z'))

        sigma_interface = get_numpy_array(
            state['sigma_on_interface_levels'].to_units(''),
            out_dims=('x', 'y', 'z'))

        tau = DataArray(
            get_frierson_06_tau(
                lat, sigma_interface, self._tau0e, self._tau0p, self._fl),
            dims=combine_dimensions(
                [state['latitude'], state['sigma_on_interface_levels']],
                out_dims=('x', 'y', 'z')),
            attrs={'units': ''}).squeeze()

        return {
            'longwave_optical_depth_on_interface_levels': tau,
        }


@jit(nopython=True)
def integrate_upward_longwave(T, T_surface, tau, sigma):
    """
    Args:
        T: 3D x-y-z air temperature array in Kelvin where z starts at the
            bottom, and z is on mid levels.
        T_surface: 2D x-y surface temperature array in Kelvin
        tau: 3D x-y-z optical depth array where z starts at the bottom, and z
            is on interface levels.
        sigma: Stefann-Boltzmann constant

    Returns:
        upward_flux: 3D x-y-z longwave radiative flux array where z starts
            at the bottom, and z is on interface levels. Positive means
            upward.
    """
    upward_flux = np.zeros(
        (T.shape[0], T.shape[1], T.shape[2]+1), dtype=np.float32)
    upward_flux[:, :, 0] = sigma*T_surface**4
    for k in range(1, T.shape[2]+1):
        dtau = tau[:, :, k] - tau[:, :, k-1]
        upward_flux[:, :, k] = (
            upward_flux[:, :, k-1] * np.exp(-dtau) +
            sigma * T[:, :, k-1]**4 * (1. - np.exp(-dtau)))
    return upward_flux


@jit(nopython=True)
def integrate_downward_longwave(T, tau, sigma):
    """
    Args:
        T: 3D x-y-z air temperature array in Kelvin where z starts at the
            bottom, and z is on mid levels.
        tau: 3D x-y-z optical depth array where z starts at the bottom, and z
            is on interface levels.
        sigma: Stefann-Boltzmann constant

    Returns:
        downward_flux: 3D x-y-z longwave radiative flux array where z starts
            at the bottom, and z is on interface levels (interfaces). Positive means
            downward.
    """
    downward_flux = np.zeros(
        (T.shape[0], T.shape[1], T.shape[2]+1), dtype=np.float32)
    for k in range(T.shape[2]-1, -1, -1):
        dtau = tau[:, :, k+1] - tau[:, :, k]
        downward_flux[:, :, k] = (
            downward_flux[:, :, k+1]*np.exp(-dtau) +
            sigma * T[:, :, k]**4 * (1 - np.exp(-dtau)))
    return downward_flux


@jit(nopython=True)
def get_longwave_fluxes(
        T, p_interface, T_surface, tau, sigma, g, Cpd):
    upward_flux = integrate_upward_longwave(T, T_surface, tau, sigma)
    downward_flux = integrate_downward_longwave(T, tau, sigma)
    net_lw_flux = upward_flux - downward_flux
    longwave_temperature_tendency = g/Cpd * (
        net_lw_flux[:, :, 1:] - net_lw_flux[:, :, :-1])/(
        p_interface[:, :, 1:] - p_interface[:, :, :-1])
    return (downward_flux, upward_flux, net_lw_flux,
            longwave_temperature_tendency, tau)


@jit(nopython=True)
def get_frierson_06_tau(latitude, sigma, tau0e, tau0p, fl):
    tau_0 = tau0e + (tau0p - tau0e) * np.sin(latitude/360.)**2
    tau = tau_0 * (1 - (fl*sigma + (1 - fl)*sigma**4))
    return tau

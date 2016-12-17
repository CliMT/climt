from .._core.base_components import Prognostic
from .._core.array import DataArray
from .._core.util import jit, replace_none_with_default, ensure_horizontal_only
from .._core.util import get_3d_numpy_array as to_3d_array
import numpy as np


class GrayLongwaveRadiation(Prognostic):

    def __init__(
            self,
            optical_depth=None,
            stefan_boltzmann=None,
            gravitational_acceleration=None,
            heat_capacity_of_dry_air_at_constant_pressure=None):
        """

        Args:
            optical_depth (DataArray): The optical depth $\tau$ of the
                atmosphere, with the greatest values at the surface.
            stefan_boltzmann: Stefan-Boltzmann constant $\sigma$ in
                $W m^{-2} K^{-4}$.
                Default taken from climt.default_constants.
            gravitational_acceleration: Gravitational acceleration in
                $m s^{-2}$.
                Default taken from climt.default_constants.
            heat_capacity_of_dry_air_at_constant_pressure: Heat capacity of
                dry air at constnat pressure in $J kg^{-1} K^{-1}$.
                Default taken from climt.default_constants.
        """
        if optical_depth is not None:
            self._optical_depth = to_3d_array(optical_depth.to_units(''))
        else:
            self._optical_depth = None
        self._stefan_boltzmann = replace_none_with_default(
            'stefan_boltzmann', stefan_boltzmann)
        self._g = replace_none_with_default(
            'gravitational_acceleration', gravitational_acceleration)
        self._Cpd = replace_none_with_default(
            'heat_capacity_of_dry_air_at_constant_pressure',
            heat_capacity_of_dry_air_at_constant_pressure)

    def __call__(self, state):
        """
        Gets tendencies and diagnostics from the passed model state.

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
        if self._optical_depth is None:
            tau = to_3d_array(state['optical_depth'].to_units(''))
        else:
            tau = self._optical_depth
        T = to_3d_array(state['air_temperature'].to_units('degK'))
        p_interface = to_3d_array(
            state['air_pressure_on_interface_levels'].to_units('Pa'))
        Ts = state['surface_temperature'].to_units('degK')  # x, y
        ensure_horizontal_only(Ts, 'surface_temperature')
        Ts = Ts.values
        (downward_flux, upward_flux, net_lw_flux,
         lw_temperature_tendency, tau) = get_longwave_fluxes(
            T, p_interface, Ts, tau, self._stefan_boltzmann,
            self._g, self._Cpd)
        dims_mid = state['air_temperature'].dims
        dims_interface = (
            list(state['air_temperature'].dims[:2]) + ['interface_levels'])
        diagnostics = {
            'downward_longwave_flux': DataArray(
                downward_flux, dims=dims_interface, attrs={'units': 'W m^-2'}),
            'upward_longwave_flux': DataArray(
                upward_flux, dims=dims_interface, attrs={'units': 'W m^-2'}),
            'net_longwave_flux': DataArray(
                net_lw_flux, dims=dims_interface, attrs={'units': 'W m^-2'}),
            'longwave_heating_rate': DataArray(
                lw_temperature_tendency, dims=dims_mid,
                attrs={'units': 'K s^-1'}),
        }
        tendencies = {
            'air_temperature': DataArray(
                lw_temperature_tendency, dims=dims_mid,
                attrs={'units': 'K s^-1'})
        }
        return tendencies, diagnostics


class Frierson06GrayLongwaveRadiation(Prognostic):

    def __init__(
            self, linear_optical_depth_parameter=0.1,
            longwave_optical_depth_at_equator=6,
            longwave_optical_depth_at_poles=1.5,
            stefan_boltzmann=None,
            gravitational_acceleration=None,
            heat_capacity_of_dry_air_at_constant_pressure=None):
        """

        Args:
            linear_optical_depth_parameter (float): The constant $f_l$ which
                determines how much of the variation of $\tau$ with pressure
                is linear rather than quartic.
                $\tau = \tau_0 [f_l \frac{p}{p_s} + (1 - f_l) (\frac{p}{p_s})^4]$
                Default is 0.1 as in Frierson et al., 2006.
            longwave_optical_depth_at_equator (float): The value of $\tau_0$
                at the equator.
                Default is 6 as in Frierson et al. 2006.
            longwave_optical_depth_at_poles (float): The value of $\tau_0$
                at the poles.
                Default is 1.5 as in Frierson et al., 2006.
            stefan_boltzmann: Stefan-Boltzmann constant $\sigma$ in
                $W m^{-2} K^{-4}$.
                Default taken from climt.default_constants.
            gravitational_acceleration: Gravitational acceleration in
                $m s^{-2}$.
                Default taken from climt.default_constants.
            heat_capacity_of_dry_air_at_constant_pressure: Heat capacity of
                dry air at constnat pressure in $J kg^{-1} K^{-1}$.
                Default taken from climt.default_constants.
        """
        self._fl = linear_optical_depth_parameter
        self._tau0e = longwave_optical_depth_at_equator
        self._tau0p = longwave_optical_depth_at_poles
        self._gray_radiation = GrayLongwaveRadiation(
            stefan_boltzmann=stefan_boltzmann,
            gravitational_acceleration=gravitational_acceleration,
            heat_capacity_of_dry_air_at_constant_pressure=heat_capacity_of_dry_air_at_constant_pressure)

    def __call__(self, state):
        """
        Gets tendencies and diagnostics from the passed model state.

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
        lat = to_3d_array(state['latitude'].to_units('degrees_north'))
        sigma_interface = to_3d_array(
            state['sigma_on_interface_levels'].to_units(''))
        tau = DataArray(
            get_frierson_06_tau(
                lat, sigma_interface, self._tau0e, self._tau0p, self._fl),
            dims=state['sigma_on_interface_levels'].dims,
            attrs={'units': ''})
        state_with_tau = state.copy()
        state_with_tau['optical_depth'] = tau
        tendencies, diagnostics = self._gray_radiation(state_with_tau)
        diagnostics['optical_depth'] = tau
        return tendencies, diagnostics


@jit(nopython=True)
def get_interface_pressures(p, ps):
    """Given 3D pressure on model half levels (cell centers) and the 2D surface
    pressure, return the 3D pressure on model full levels (cell interfaces).
    If the z-dimension of p is length K, the returned p_full will have a
    z-dimension of length K+1."""
    interface_pressures = np.zeros(
        (p.shape[0], p.shape[1], p.shape[2]+1), dtype=np.float32)
    interface_pressures[:, :, 1:-1] = 0.5*(p[:, :, 1:] + p[:, :, :-1])
    interface_pressures[:, :, 0] = ps[:, :, 0]
    return interface_pressures


@jit(nopython=True)
def integrate_upward_longwave(T, T_surface, tau, sigma):
    """
    Args:
        T: 3D x-y-z air temperature array in Kelvin where z starts at the
            bottom, and z is on half levels.
        T_surface: 2D x-y surface temperature array in Kelvin
        tau: 3D x-y-z optical depth array where z starts at the bottom, and z
            is on half levels.
        sigma: Stefann-Boltzmann constant

    Returns:
        upward_flux: 3D x-y-z longwave radiative flux array where z starts
            at the bottom, and z is on full levels (interfaces). Positive means
            upward.
    """
    upward_flux = np.zeros(
        (T.shape[0], T.shape[1], T.shape[2]+1), dtype=np.float32)
    upward_flux[:, :, 0] = sigma*T_surface**4
    for k in range(1, T.shape[2]+1):
        dtau = tau[:, :, k] - tau[:, :, k-1]
        upward_flux[:, :, k] = (
            upward_flux[:, :, k-1] * np.exp(dtau) +
            sigma * T[:, :, k-1]**4 * (1. - np.exp(dtau)))
    return upward_flux


@jit(nopython=True)
def integrate_downward_longwave(T, tau, sigma):
    """
    Args:
        T: 3D x-y-z air temperature array in Kelvin where z starts at the
            bottom, and z is on half levels.
        tau: 3D x-y-z optical depth array where z starts at the bottom, and z
            is on half levels.
        sigma: Stefann-Boltzmann constant

    Returns:
        downward_flux: 3D x-y-z longwave radiative flux array where z starts
            at the bottom, and z is on full levels (interfaces). Positive means
            downward.
    """
    downward_flux = np.zeros(
        (T.shape[0], T.shape[1], T.shape[2]+1), dtype=np.float32)
    for k in range(T.shape[2]-2, -1, -1):
        dtau = tau[:, :, k] - tau[:, :, k+1]
        downward_flux[:, :, k] = (
            downward_flux[:, :, k+1]*np.exp(-dtau) -
            sigma * T[:, :, k]**4 * (np.exp(-dtau) - 1.))
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
    tau = tau_0 * (fl*sigma + (1 - fl)*sigma**4)
    return tau

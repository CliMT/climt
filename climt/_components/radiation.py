from sympl import (
    Prognostic, Diagnostic, DataArray, jit, replace_none_with_default,
    combine_dimensions, get_numpy_array)
import numpy as np


class GrayLongwaveRadiation(Prognostic):

    inputs = (
        'longwave_optical_depth_on_interface_levels', 'air_temperature',
        'air_pressure_on_interface_levels', 'surface_temperature')
    diagnostics = (
        'downward_longwave_flux', 'upward_longwave_flux',
        'net_longwave_flux')
    tendencies = ('air_temperature',)

    def __init__(
            self,
            stefan_boltzmann=None,
            gravitational_acceleration=None,
            heat_capacity_of_dry_air_at_constant_pressure=None):
        """

        Args:
            longwave_optical_depth_on_interface_levels (DataArray): The optical
                depth $\tau$ of the atmosphere, with the lowest values at the
                surface.
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

        self._stefan_boltzmann = replace_none_with_default(
            'stefan_boltzmann_constant', stefan_boltzmann)
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
            T, p_interface, Ts, tau, self._stefan_boltzmann.values,
            self._g.values, self._Cpd.values)
        dims_mid = combine_dimensions(
            [state['surface_temperature'],
             state['air_temperature']],
            out_dims=('x', 'y', 'z'))
        dims_interface = combine_dimensions(
            [state['air_pressure_on_interface_levels'],
             state['surface_temperature']],
            out_dims=('x', 'y', 'z'))
        diagnostics = {
            'downward_longwave_flux': DataArray(
                downward_flux, dims=dims_interface, attrs={'units': 'W m^-2'}
            ).squeeze(),
            'upward_longwave_flux': DataArray(
                upward_flux, dims=dims_interface, attrs={'units': 'W m^-2'}
            ).squeeze(),
            'net_longwave_flux': DataArray(
                net_lw_flux, dims=dims_interface, attrs={'units': 'W m^-2'}
            ).squeeze(),
        }
        tendencies = {
            'air_temperature': DataArray(
                lw_temperature_tendency, dims=dims_mid,
                attrs={'units': 'K s^-1'}).squeeze()
        }
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
            latitude (DataArray, optional): The constant latitude coordinate to
                use for calculations. If given, latitude in state is not used.
            sigma_on_interface_levels (DataArray, optional): The constant sigma
                vertical coordinate to use for calculations. If given, its
                value is ignored in input states, so only provide it if the
                model is using constant sigma coordinates.
            linear_optical_depth_parameter (float, optional): The constant $f_l$ which
                determines how much of the variation of $\tau$ with pressure
                is linear rather than quartic.
                $\tau = \tau_0 [f_l \frac{p}{p_s} + (1 - f_l) (\frac{p}{p_s})^4]$
                Default is 0.1 as in Frierson et al., 2006.
            longwave_optical_depth_at_equator (float, optional): The value of $\tau_0$
                at the equator.
                Default is 6 as in Frierson et al. 2006.
            longwave_optical_depth_at_poles (float, optional): The value of $\tau_0$
                at the poles.
                Default is 1.5 as in Frierson et al., 2006.
        """
        self._fl = linear_optical_depth_parameter
        self._tau0e = longwave_optical_depth_at_equator
        self._tau0p = longwave_optical_depth_at_poles

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

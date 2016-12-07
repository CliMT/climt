from ..core.base_components import Prognostic
from ..core.array import DataArray
from ..core.util import(
    ensure_shared_coordinates, ensure_third_dim_is_vertical,
    ensure_horizontal_only, ensure_number_of_dims)
from ..core.exceptions import InvalidStateException
from ..core.units import to_units
from numba import jit
import numpy as np


class Frierson06GrayLongwaveRadiation(Prognostic):

    def __init__(
            self, linear_optical_depth_parameter=0.1,
            longwave_optical_depth_at_equator=6,
            longwave_optical_depth_at_poles=1.5,
            sigma=5.6734e-8,  # W/m^2/K^4
            g=9.80665,  # m/s^2
            Cpd=1004.64,  #J/kg/K
            ):
        """

        Args:
            linear_optical_depth_parameter: The constant $f_l$ which determines
                how much of the variation of $\tau$ with pressure is linear
                rather than quartic.
                $\tau = \tau_0 [f_l \frac{p}{p_s} + (1 - f_l) (\frac{p}{p_s})^4]$
                Default is 0.1 as in Frierson et al., 2006.
            longwave_optical_depth_at_equator: The value of $\tau_0$ at the
                equator.
                Default is 6 as in Frierson et al. 2006.
            longwave_optical_depth_at_poles: The value of $\tau_0$ at the
                poles.
                Default is 1.5 as in Frierson et al., 2006.
            sigma: Stefan-Boltzmann constant $\sigma$ in $W m^{-2} K^{-4}$.
                Default is 5.6734e-8
            g: Gravitational acceleration in $m s^{-2}$. Default is 9.80665.
            Cpd: Heat capacity of dry air in $J kg^{-1} K^{-1}$.
                Default is 1004.64
        """
        self._fl = linear_optical_depth_parameter
        self._tau0e = longwave_optical_depth_at_equator
        self._tau0p = longwave_optical_depth_at_poles
        self._sigma = sigma
        self._g = g
        self._Cpd = Cpd

    def ensure_state_is_valid_input(self, state):
        """
        Args:
            state (dict): A model state dictionary. Must contain the following
                as climt.DataArray objects with valid units:

                * air_temperature
                * air_pressure
                * surface_pressure
                * surface_temperature

                Dimension orders must be (longitude, latitude, height), with
                height starting at the surface.

        Raises:
            KeyError: if a required quantity is missing from the state
            TypeError: if a required quantity is not a climt.DataArray
            InvalidStateException: if state is not a valid input for the
                Prognostic instance.
        """
        try:
            T = state['air_temperature'].to_units('degK').values
            p = state['air_pressure'].to_units('Pa').values  # x, y, z
            ps = state['surface_pressure'].to_units('Pa').values  # x, y
            Ts = state['surface_temperature'].to_units('degK').values  # x, y
            lat = to_units(T.coords['latitude'], 'degrees_north').values  # y
        except AttributeError:
            raise TypeError(
                'provided quantities must be climt.DataArray objects')
        try:
            ensure_shared_coordinates(T, p)
        except AssertionError:
            raise InvalidStateException(
                'air_temperature and air_pressure must share '
                'the same coordinates')
        try:
            ensure_shared_coordinates(Ts, ps)
        except AssertionError:
            raise InvalidStateException(
                'surface_temperature and surface_pressure must share '
                'the same coordinates')
        for name in 'air_temperature', 'air_pressure':
            ensure_third_dim_is_vertical(state[name], name)
        for name in 'surface_pressure', 'surface_temperature':
            ensure_horizontal_only(state[name], name)
        for name in ('surface_temperature', 'surface_pressure'):
            ensure_number_of_dims(state[name], num_dims=2, name=name)

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
        T = state['air_temperature'].to_units('degK').values
        p = state['air_pressure'].to_units('Pa').values  # x, y, z
        ps = state['surface_pressure'].to_units('Pa').values  # x, y
        Ts = state['surface_temperature'].to_units('degK').values  # x, y
        lat = state['latitude'].to_units('degrees_north').values  # y
        (downward_flux, upward_flux, net_lw_flux,
         lw_temperature_tendency, tau) = get_longwave_fluxes(
            T, p, ps[:, :, None], Ts, lat[None, :, None], self._tau0e,
            self._tau0p, self._fl, self._sigma, self._g, self._Cpd)
        dims_half = state['air_temperature'].dims
        dims_full = (
            list(state['air_temperature'].dims[:2]) + ['full_levels'])
        diagnostics = {
            'downward_longwave_flux': DataArray(
                downward_flux, dims=dims_full, attrs={'units': 'W m^-2'}),
            'upward_longwave_flux': DataArray(
                upward_flux, dims=dims_full, attrs={'units': 'W m^-2'}),
            'net_longwave_flux': DataArray(
                net_lw_flux, dims=dims_full, attrs={'units': 'W m^-2'}),
            'optical_depth': DataArray(
                tau, dims=dims_half, attrs={'units': ''}),
            'longwave_heating_rate': DataArray(
                lw_temperature_tendency, dims=dims_half,
                attrs={'units': 'K s^-1'}),
        }
        tendencies = {
            'air_temperature': DataArray(
                lw_temperature_tendency, dims=dims_half,
                attrs={'units': 'K s^-1'})
        }
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
        T, p, ps, T_surface, latitude, tau0e, tau0p, fl, sigma, g, Cpd):
    pi = get_interface_pressures(p, ps)
    tau_0 = tau0e + (tau0p - tau0e) * np.sin(latitude/360.)**2
    tau = tau_0 * (
        fl*pi/ps + (1 - fl)*(pi/ps)**4)
    upward_flux = integrate_upward_longwave(T, T_surface, tau, sigma)
    downward_flux = integrate_downward_longwave(T, tau, sigma)
    net_lw_flux = upward_flux - downward_flux
    longwave_temperature_tendency = g/Cpd * (
        net_lw_flux[:, :, 1:] - net_lw_flux[:, :, :-1])/(
        pi[:, :, 1:] - pi[:, :, :-1])
    return (downward_flux, upward_flux, net_lw_flux,
            longwave_temperature_tendency, tau)

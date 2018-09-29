from sympl import TendencyComponent, DiagnosticComponent, get_constant
import numpy as np


class GrayLongwaveRadiation(TendencyComponent):

    input_properties = {
        'longwave_optical_depth_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'dimensionless',
            'alias': 'tau',
        },
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK',
            'alias': 'sl',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
            'alias': 'T_surface',
        },
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa',
            'alias': 'p',
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'Pa',
            'alias': 'p_interface',
        },
    }

    diagnostic_properties = {
        'downwelling_longwave_flux_in_air': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
            'alias': 'lw_down',
        },
        'upwelling_longwave_flux_in_air': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
            'alias': 'lw_up',
        },
        'longwave_heating_rate': {
            'dims': ['mid_levels', '*'],
            'units': 'degK day^-1'
        }
    }

    tendency_properties = {
        'air_temperature': {'units': 'degK s^-1'},
    }

    def array_call(self, state):
        (downward_flux, upward_flux, net_lw_flux,
         lw_temperature_tendency, tau) = get_longwave_fluxes(
            state['sl'], state['p_interface'], state['T_surface'], state['tau'],
            get_constant('stefan_boltzmann_constant', 'W/m^2/K^4'),
            get_constant('gravitational_acceleration', 'm/s^2'),
            get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K')
        )
        tendencies = {
            'sl': lw_temperature_tendency,
        }
        diagnostics = {
            'lw_down': downward_flux,
            'lw_up': upward_flux,
            'longwave_heating_rate': lw_temperature_tendency * 86400.
        }
        return tendencies, diagnostics


class Frierson06LongwaveOpticalDepth(DiagnosticComponent):

    input_properties = {
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'Pa',
        },
        'surface_air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'latitude': {
            'dims': ['*'],
            'units': 'degrees_N',
        }
    }

    diagnostic_properties = {
        'longwave_optical_depth_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'dimensionless',
        }
    }

    def __init__(
            self,
            linear_optical_depth_parameter=0.1,
            longwave_optical_depth_at_equator=6,
            longwave_optical_depth_at_poles=1.5,
            **kwargs):
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
        super(Frierson06LongwaveOpticalDepth, self).__init__(**kwargs)

    def array_call(self, state):
        return {
            'longwave_optical_depth_on_interface_levels': get_frierson_06_tau(
                state['latitude'],
                state['air_pressure_on_interface_levels'] /
                state['surface_air_pressure'][None, :],
                self._tau0e, self._tau0p, self._fl
            )
        }


# @jit(nopython=True)
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
        (T.shape[0]+1, T.shape[1]), dtype=np.float32)
    upward_flux[0, :] = sigma*T_surface**4
    for k in range(1, T.shape[0]+1):
        dtau = tau[k, :] - tau[k-1, :]
        upward_flux[k, :] = (
            upward_flux[k-1, :] * np.exp(-dtau) +
            sigma * T[k-1, :]**4 * (1. - np.exp(-dtau)))
    return upward_flux


# @jit(nopython=True)
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
        (T.shape[0]+1, T.shape[1]), dtype=np.float32)
    for k in range(T.shape[0]-1, -1, -1):
        dtau = tau[k+1, :] - tau[k, :]
        downward_flux[k, :] = (
            downward_flux[k+1, :]*np.exp(-dtau) +
            sigma * T[k, :]**4 * (1 - np.exp(-dtau)))
    return downward_flux


# @jit(nopython=True)
def get_longwave_fluxes(
        T, p_interface, T_surface, tau, sigma, g, Cpd):
    upward_flux = integrate_upward_longwave(T, T_surface, tau, sigma)
    downward_flux = integrate_downward_longwave(T, tau, sigma)
    net_lw_flux = upward_flux - downward_flux
    longwave_temperature_tendency = g/Cpd * (
        net_lw_flux[1:, :] - net_lw_flux[:-1, :])/(
        p_interface[1:, :] - p_interface[:-1, :])
    return (downward_flux, upward_flux, net_lw_flux,
            longwave_temperature_tendency, tau)


# @jit(nopython=True)
def get_frierson_06_tau(latitude, sigma, tau0e, tau0p, fl):
    tau_0 = tau0e + (tau0p - tau0e) * np.sin(latitude*np.pi/180.)**2
    tau = tau_0 * (1 - (fl*sigma + (1 - fl)*sigma**4))
    return tau

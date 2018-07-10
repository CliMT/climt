from sympl import jit, DataArray
import numpy as np
import functools


def ensure_contiguous_state(func):
    @functools.wraps(func)
    def wrapper(self, state, *args, **kwargs):
        for name, value in state.items():
            if isinstance(value, np.ndarray):
                state[name] = np.ascontiguousarray(value)
        return func(self, state, *args, **kwargs)
    return wrapper


def numpy_version_of(state):
    """
    Return raw numpy arrays from state.

    Args:
        state (dict):
            state dictionary containing DataArrays

    Returns:
        raw_arrays (dict):
            dictionary containing the numpy arrays
            underlying the DataArrays.
    """

    raw_arrays = {}

    for quantity in state.keys():
        if isinstance(state[quantity], DataArray):
            raw_arrays[quantity] = state[quantity].values
        elif isinstance(state[quantity], np.ndarray):
            raw_arrays[quantity] = state[quantity]

    return raw_arrays


def mass_to_volume_mixing_ratio(
        mass_mixing_ratio,
        molecular_weight=None,
        molecular_weight_air=28.964):
    """
    g/g or g/kg to mole/mole

    Converts from mass mixing ratio (mass per unit mass) to volume
    mixing ratio (volume per unit volume)

    Args:

        mass_mixing_ratio (array):
            The quantity to be transformed in units of :math:`g/g`.

        molecular_weight (float):
            The molecular weight of the gas in :math:`g/mol`.

        molecular_weight_air (float,optional):
            The molecular weight of dry air.
            If it is not provided, the value for dry air on earth (28.964 g/mol)
            is used.

    Returns:

        volume_mixing_ratio (array):
            The volume mixing ratio of the gas.

    Raises:

        ValueError:
            if the molecular weight is not provided.

    """

    if molecular_weight is None:
        raise ValueError('The molecular weight must be provided')

    volume_mixing_ratio = mass_mixing_ratio*molecular_weight_air/molecular_weight

    return volume_mixing_ratio


def get_interface_values(
        mid_level_values,
        surface_value,
        mid_level_pressure,
        interface_level_pressure):
    """
    Calculate interface values given mid-level values.

    Given 3D values of a quantity on model mid levels (cell centers) and the 2D surface
    value, return the 3D values of that quantity on model full levels (cell interfaces).
    If the vertical dimension of :code:`mid_level_values` is length K, the returned array will have a
    vertical dimension of length K+1.

    Routine borrowed from CESM (radiation.F90 in rrtmg folder)

    Args:

        mid_level_values (array):
            The values of the quantity on mid-levels.

        surface_value (array):
            The value of the quantity at the surface. Must be in the
            same units as :code:`mid_level_values`

        mid_level_pressure (array):
            Pressure values on mid-levels. Can be in any units.

        interface_level_pressure (array):
            Pressure values on interface levels. Must be in
            in the same units as :code:`mid_level_pressure`.

    Returns:

        interface_values (array):
            values of the quantity on mid-levels.

    """

    interface_values = np.zeros(
        (mid_level_values.shape[0]+1, mid_level_values.shape[1]), dtype=np.double)

    log_mid_p = np.log(mid_level_pressure)

    interp_weight = (
        np.log(interface_level_pressure[1:-1, :]) - log_mid_p[1:, :]) / (
        log_mid_p[:-1, :] - log_mid_p[1::, :])

    interface_values[1:-1, :] = \
        mid_level_values[1:, :] - interp_weight*(
            mid_level_values[1:, :] - mid_level_values[0:-1, :])

    interface_values[0, :] = surface_value[:]
    interface_values[-1, :] = mid_level_values[-1, :]

    return interface_values


def calculate_q_sat(surf_temp, surf_press, Rd, Rv):

    above_freezing = surf_temp > 273
    below_freezing = np.logical_not(above_freezing)

    sat_vap_press = np.zeros(surf_temp.shape)

    eps = Rd/Rv

    sat_vap_press[above_freezing] =\
        (1.0007 + 3.46e-8*surf_press[above_freezing])*611.21*np.exp(
            17.966*(surf_temp[above_freezing] - 273.)/(247.15 + (surf_temp[above_freezing] - 273.)))

    sat_vap_press[below_freezing] = \
        (1.0003 + 4.18e-8*surf_press[below_freezing])*611.15*np.exp(
            22.452*(surf_temp[below_freezing] - 273.)/(272.5 + (surf_temp[below_freezing] - 273.)))

    return eps*sat_vap_press/(surf_press - (1 - eps)*sat_vap_press)


@jit(nopython=True)
def bolton_q_sat(T, p, Rd, Rh2O):
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    epsilon = Rd/Rh2O
    return epsilon*es/(p - (1 - epsilon)*es)


@jit(nopython=True)
def bolton_dqsat_dT(T, Lv, Rh2O, q_sat):
    """Uses the assumptions of equation 12 in Reed and Jablonowski, 2012. In
    particular, assumes d(qsat)/dT is approximately epsilon/p*d(es)/dT"""
    return Lv*q_sat/(Rh2O*T**2)

from pint import UnitRegistry
from sympl import DataArray, jit
import numpy as np

pint_units = UnitRegistry()


def mass_to_volume_mixing_ratio(
        mass_mixing_ratio,
        molecular_weight=None,
        molecular_weight_air=28.964):
    """
    g/g or g/kg to mole/mole

    Converts from mass mixing ratio (mass per unit mass) to volume
    mixing ratio (volume per unit volume)

    Args:

        mass_mixing_ratio (DataArray):
            The quantity to be transformed. It
            must be a DataArray with the units attribute correctly set (i.e, to
            something like 'g/kg' or 'kg/kg').

        molecular_weight (float):
            The molecular weight of the gas in grams/mole.

        molecular_weight_air (float,optional):
            The molecular weight of dry air.
            If it is not provided, the value for dry air on earth (28.964 g/mol)
            is used.

    Returns:

        volume_mixing_ratio (DataArray):
            The volume mixing ratio of the gas.

    Raises:

        ValueError:
            if the molecular weight is not provided.

    """

    if molecular_weight is None:
        raise ValueError('The molecular weight must be provided')

    mass_mixing_ratio_with_units = mass_mixing_ratio.values*pint_units(mass_mixing_ratio.units)

    dims = mass_mixing_ratio.dims

    volume_mixing_ratio = mass_mixing_ratio_with_units*molecular_weight_air/molecular_weight

    volume_mixing_ratio = volume_mixing_ratio.to_base_units()
    volume_mixing_ratio = DataArray(volume_mixing_ratio,
                                    dims=dims,
                                    attrs={'units': str(volume_mixing_ratio.units)})

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
        (mid_level_values.shape[0], mid_level_values.shape[1], mid_level_values.shape[2]+1),
        dtype=np.double, order='F')

    log_mid_p = np.log(mid_level_pressure)
    log_int_p = np.log(interface_level_pressure)

    interp_weight = (log_int_p[:, :, 1:-1] - log_mid_p[:, :, 1::])/(log_mid_p[:, :, :-1] -
                                                                    log_mid_p[:, :, 1::])

    interface_values[:, :, 1:-1] = \
        mid_level_values[:, :, 1::] - interp_weight*(mid_level_values[:, :, 1::] -
                                                     mid_level_values[:, :, 0:-1])

    interface_values[:, :, 0] = surface_value[:, :]
    interface_values[:, :, -1] = mid_level_values[:, :, -1]

    return interface_values


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

from pint import UnitRegistry
from sympl import DataArray

pint_units = UnitRegistry()

def mass_to_volume_mixing_ratio(
    mass_mixing_ratio,
    molecular_weight=None,
    molecular_weight_air=28.964):
    """
    Converts from mass mixing ratio (mass per unit mass) to volume
    mixing ratio (volume per unit volume)

    Args:

    mass_mixing_ratio (DataArray): The quantity to be transformed. It
        must be a DataArray with the units attribute correctly set (i.e, to
        something like 'g/kg' or 'kg/kg').

    molecular_weight (float): The molecular weight of the gas in grams/mole.

    molecular_weight_air (float,optional): The molecular weight of dry air.
        If it is not provided, the value for dry air on earth (28.964 g/mol)
        is used.

    Returns:

    volume_mixing_ratio (DataArray): The volume mixing ratio of the gas.
    """

    if molecular_weight is None:
        raise ValueError('The molecular weight must be provided')

    mass_mixing_ratio_with_units = mass_mixing_ratio.values*pint_units(mass_mixing_ratio.units)

    dims = mass_mixing_ratio.dims

    volume_mixing_ratio = mass_mixing_ratio_with_units*molecular_weight_air/molecular_weight

    volume_mixing_ratio = DataArray(volume_mixing_ratio.to_base_units(), 
                                    dims=dims, attrs={'units': str(volume_mixing_ratio.units)})

    return volume_mixing_ratio

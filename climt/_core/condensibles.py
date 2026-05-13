from typing import NamedTuple

import numpy as np
from sympl import get_constant, set_constant

from .atmospheric_properties import get_active_condensible

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(x, **kwargs):
        return x

SPECIES_ID = {"h2o": 0, "ch4": 1, "co2": 2}

# Register H2O condensible constants as sympl defaults.
# load_atmospheric_properties() overrides these when a profile is loaded.
_H2O_DEFAULTS = {
    "gas_constant_of_condensible_vapor":       (461.5,   "J/kg/K"),
    "heat_capacity_of_condensible_vapor":      (1870.0,  "J/kg/K"),
    "heat_capacity_of_condensible_liquid":     (2500.0,  "J/kg/K"),
    "latent_heat_of_condensible_vaporization": (2.501e6, "J/kg"),
    "density_of_condensible_liquid":           (1000.0,  "kg/m^3"),
    "freezing_point_of_condensible":           (273.15,  "degK"),
}
for _name, (_val, _units) in _H2O_DEFAULTS.items():
    set_constant(_name, _val, _units)


class CondensibleParams(NamedTuple):
    species_id: int
    RV: float        # gas constant of condensible vapor (J/kg/K)
    CPV: float       # heat capacity of condensible vapor (J/kg/K)
    CL: float        # heat capacity of condensible liquid/solid (J/kg/K)
    LV0: float       # latent heat of vaporization at T_freeze (J/kg)
    ROWL: float      # density of condensible liquid/solid (kg/m^3)
    T_freeze: float  # phase transition temperature (K)


def get_condensible_params() -> CondensibleParams:
    """Read active condensible constants from sympl and return a CondensibleParams.

    Call once per component __init__. Requires either H2O defaults (registered
    at module import) or a profile loaded via load_atmospheric_properties().
    """
    species = get_active_condensible()
    return CondensibleParams(
        species_id=SPECIES_ID[species],
        RV=get_constant("gas_constant_of_condensible_vapor", "J/kg/K"),
        CPV=get_constant("heat_capacity_of_condensible_vapor", "J/kg/K"),
        CL=get_constant("heat_capacity_of_condensible_liquid", "J/kg/K"),
        LV0=get_constant("latent_heat_of_condensible_vaporization", "J/kg"),
        ROWL=get_constant("density_of_condensible_liquid", "kg/m^3"),
        T_freeze=get_constant("freezing_point_of_condensible", "degK"),
    )

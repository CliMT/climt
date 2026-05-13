from typing import NamedTuple

import numpy as np
from sympl import get_constant, set_constant

from .atmospheric_properties import get_active_condensible
from .backend import jit_compile, prange

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


@njit
def _sat_vap_pressure(T, species_id):
    """Saturation vapor pressure in hPa for the given species.

    H2O: Magnus formula (above/below freezing).
    CH4, CO2: Clausius-Clapeyron anchored at 1 atm reference point.
    """
    if species_id == 0:  # H2O
        TC = T - 273.15
        if TC >= 0.0:
            return 6.112 * np.exp(17.67 * TC / (243.5 + TC))
        else:
            return np.exp(23.33086 - 6111.72784 / T + 0.15215 * np.log(T))
    elif species_id == 1:  # CH4 — anchored at normal boiling point 111.65 K, 1 atm
        return 1013.25 * np.exp((5.1e5 / 518.3) * (1.0 / 111.65 - 1.0 / T))
    else:  # CO2 — anchored at sublimation point 194.7 K, 1 atm
        return 1013.25 * np.exp((5.71e5 / 188.9) * (1.0 / 194.7 - 1.0 / T))


@njit
def _lv(T, cond):
    """Latent heat of vaporization/sublimation (J/kg) as a function of temperature."""
    return cond.LV0 - (cond.CL - cond.CPV) * (T - cond.T_freeze)


@njit
def _lcl_pressure(P_nk, RH, T_nk, cond):
    """Lifting Condensation Level pressure (hPa).

    H2O: Bolton/Bohren empirical formula.
    CH4, CO2: Clausius-Clapeyron approximation PLCL = P * RH^(Rv*T/Lv).
    """
    if cond.species_id == 0:  # H2O empirical
        CHI = T_nk / (1669.0 - 122.0 * RH - T_nk)
        return P_nk * (RH ** CHI)
    else:  # Clausius-Clapeyron
        chi = cond.RV * T_nk / cond.LV0
        return P_nk * (RH ** chi)


@jit_compile(backend=np, parallel=True)
def compute_qs(T, P, cond, RD):
    """Saturation specific humidity (kg/kg) for a 2-D (nlev, ncol) array.

    Args:
        T:   air temperature (K), shape (nlev, ncol)
        P:   air pressure (hPa), shape (nlev, ncol)
        cond: CondensibleParams for the active condensible
        RD:  specific gas constant of dry air (J/kg/K)

    Returns:
        qs array with same shape as T.
    """
    nlev, ncol = T.shape
    qs = np.zeros_like(T)
    EPS = RD / cond.RV
    for i in prange(ncol):
        for k in range(nlev):
            es = _sat_vap_pressure(T[k, i], cond.species_id)
            qs[k, i] = EPS * es / (P[k, i] - (1.0 - EPS) * es)
    return qs

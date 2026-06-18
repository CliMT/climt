# climt/_components/cork/common.py
import numpy as np

from ..._core.backend import prange

# Molar masses (g/mol) for unit conversion between mole fraction and mass mixing ratio.
# H2O is supplied as specific humidity (kg/kg) and does not need conversion.
# All other gases come in as mole fractions and are multiplied by M_gas / M_AIR.
MOLAR_MASS_DRY_AIR = 28.970  # g/mol
MOLAR_MASS = {
    "h2o": 18.015,
    "co2": 44.010,
    "o3": 47.998,
    "ch4": 16.043,
    "n2o": 44.013,
    "o2": 31.998,
}

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(*args, **kwargs):
        # No-numba fallback. Support both bare ``@njit`` and parametrized
        # ``@njit(parallel=True, ...)`` usage so the pure-Python path imports.
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def _decorator(func):
            return func

        return _decorator


@njit
def compute_heating_rate(net_flux, p_interface, g, cpd):
    """Compute heating rate (K/s) from net flux divergence.

    Args:
        net_flux: (nlev+1, ncol) upward minus downward flux, W/m^2
        p_interface: (nlev+1, ncol) pressure at interfaces, Pa
        g: gravitational acceleration, m/s^2
        cpd: heat capacity of dry air at constant pressure, J/kg/K

    Returns:
        heating_rate: (nlev, ncol) in K/s
    """
    nlev = net_flux.shape[0] - 1
    ncol = net_flux.shape[1]
    heating_rate = np.zeros((nlev, ncol))
    for i in prange(ncol):
        for k in range(nlev):
            dp = p_interface[k + 1, i] - p_interface[k, i]
            dflux = net_flux[k + 1, i] - net_flux[k, i]
            heating_rate[k, i] = g / cpd * dflux / dp
    return heating_rate


@njit
def compute_column_amount(q, p_interface, g):
    """Compute column amount (kg/m^2) of a gas in each layer.

    Args:
        q: (nlev, ncol) mass mixing ratio, kg/kg
        p_interface: (nlev+1, ncol) pressure at interfaces, Pa
        g: gravitational acceleration, m/s^2

    Returns:
        amount: (nlev, ncol) column amount, kg/m^2
    """
    nlev, ncol = q.shape
    amount = np.zeros((nlev, ncol))
    for i in prange(ncol):
        for k in range(nlev):
            dp = abs(p_interface[k + 1, i] - p_interface[k, i])
            amount[k, i] = q[k, i] * dp / g
    return amount

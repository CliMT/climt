# climt/_components/picket_fence/common.py
import numpy as np

from ..._core.backend import prange

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(x, **kwargs):
        return x


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

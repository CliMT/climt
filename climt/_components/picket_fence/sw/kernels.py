# climt/_components/picket_fence/sw/kernels.py
import numpy as np

from ..common import njit, prange


@njit
def _sw_direct_beam(tau_col, mu0, solar_flux_col):
    """Compute direct beam attenuation through a column.

    Args:
        tau_col: (nlev,) optical depth per layer
        mu0: cosine of zenith angle (scalar)
        solar_flux_col: incoming solar flux at TOA (scalar, W/m^2)

    Returns:
        direct_down: (nlev+1,) direct downward flux at each interface
    """
    nlev = tau_col.shape[0]
    direct_down = np.zeros(nlev + 1)
    direct_down[nlev] = solar_flux_col * mu0  # TOA (index nlev = top)

    cumulative_tau = 0.0
    for k in range(nlev - 1, -1, -1):
        cumulative_tau += tau_col[k]
        direct_down[k] = solar_flux_col * mu0 * np.exp(-cumulative_tau / mu0)

    return direct_down


@njit
def sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights):
    """Multi-band, multi-g-point SW radiative transfer.

    Simplified two-stream: direct beam attenuation only (no diffuse for v1).
    Scattering is accounted for by reducing the effective optical depth.

    Args:
        tau: (nband, ngpt, nlev, ncol) extinction optical depth
        ssa: (nband, ngpt, nlev, ncol) single scattering albedo
        asymmetry: (nband, ngpt, nlev, ncol) asymmetry parameter
        zenith: (ncol,) solar zenith angle, radians
        albedo: (ncol,) surface albedo
        solar_flux: (nband, ngpt) TOA solar flux per g-point, W/m^2
        weights: (nband, ngpt) g-point weights

    Returns:
        up_band: (nband, nlev+1, ncol) per-band upward flux
        down_band: (nband, nlev+1, ncol) per-band downward flux
        up_broad: (nlev+1, ncol) broadband upward flux
        down_broad: (nlev+1, ncol) broadband downward flux
    """
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    for b in range(nband):
        for g in range(ngpt):
            w = weights[b, g]
            for i in prange(ncol):
                mu0 = np.cos(zenith[i])
                if mu0 <= 0.0:
                    continue  # nighttime

                # Absorption optical depth: tau * (1 - ssa)
                tau_abs = np.zeros(nlev)
                for k in range(nlev):
                    tau_abs[k] = tau[b, g, k, i] * (1.0 - ssa[b, g, k, i])

                direct = _sw_direct_beam(tau_abs, mu0, solar_flux[b, g])

                # Downward flux = direct beam (diffuse neglected in v1)
                for k in range(nlev + 1):
                    down_band[b, k, i] += w * direct[k]

                # Upward flux = reflected direct at surface
                surface_down = direct[0]
                reflected = albedo[i] * surface_down
                # Reflected flux attenuates going up (simplified)
                for k in range(nlev + 1):
                    up_band[b, k, i] += (
                        w * reflected
                    )  # simplified: no re-absorption of reflected

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    return up_band, down_band, up_broad, down_broad

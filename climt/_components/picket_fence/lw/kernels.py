# climt/_components/picket_fence/lw/kernels.py
import numpy as np

from ..common import njit, prange


@njit
def _lw_transport_single_gpt(
    tau_gpt,
    planck_gpt,
    surface_source_gpt,
    emissivity_band,
    up_flux,
    down_flux,
    nlev,
    ncol,
):
    """LW transport for a single g-point within a single band.

    Uses a diffusivity factor D=1.66 to convert vertical optical depth to
    effective diffuse optical depth, consistent with the Eddington approximation
    used in the Parmentier & Guillot (2014) analytical model.

    Args:
        tau_gpt: (nlev, ncol) optical depth per layer for this g-point
        planck_gpt: (nlev, ncol) Planck source per layer for this g-point
        surface_source_gpt: (ncol,) surface Planck source for this g-point
        emissivity_band: (ncol,) surface emissivity for this band
        up_flux: (nlev+1, ncol) output upward flux (accumulated in-place)
        down_flux: (nlev+1, ncol) output downward flux (accumulated in-place)
        nlev: number of layers
        ncol: number of columns
    """
    DIFFUSIVITY_FACTOR = 1.66
    for i in prange(ncol):
        # Upward sweep: surface to TOA
        up_flux[0, i] = emissivity_band[i] * surface_source_gpt[i]
        for k in range(nlev):
            trans = np.exp(-DIFFUSIVITY_FACTOR * tau_gpt[k, i])
            up_flux[k + 1, i] = up_flux[k, i] * trans + planck_gpt[k, i] * (1.0 - trans)

        # Downward sweep: TOA to surface
        down_flux[nlev, i] = 0.0
        for k in range(nlev - 1, -1, -1):
            trans = np.exp(-DIFFUSIVITY_FACTOR * tau_gpt[k, i])
            down_flux[k, i] = down_flux[k + 1, i] * trans + planck_gpt[k, i] * (
                1.0 - trans
            )


def lw_transport(
    T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma,
    diagnostics_level=0,
):
    """Multi-band, multi-g-point LW radiative transfer.

    Args:
        T: (nlev, ncol) air temperature, K (unused, kept for interface consistency)
        T_surface: (ncol,) surface temperature, K (unused, kept for interface consistency)
        tau: (nband, ngpt, nlev, ncol) optical depth per layer
        planck_source: (nband, ngpt, nlev, ncol) Planck source per layer per g-point, W/m^2
        surface_source: (nband, ngpt, ncol) surface Planck source per g-point, W/m^2
        emissivity: (nband, ncol) surface emissivity per band
        weights: (nband, ngpt) g-point quadrature weights
        sigma: Stefan-Boltzmann constant (unused, kept for interface consistency)
        diagnostics_level: 0 (fluxes only), 1 (per-layer transmittance + per-gpoint fluxes)

    Returns:
        If diagnostics_level == 0:
            (up_band, down_band, up_broad, down_broad)
        If diagnostics_level > 0:
            (up_band, down_band, up_broad, down_broad, diagnostics_dict)
            where diagnostics_dict contains:
                transmittance: (nband, ngpt, nlev, ncol) per-layer diffuse transmittance
                up_per_gpoint: (nband, ngpt, nlev+1, ncol) weighted upward flux per g-point
                down_per_gpoint: (nband, ngpt, nlev+1, ncol) weighted downward flux per g-point
    """
    DIFFUSIVITY_FACTOR = 1.66
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    if diagnostics_level >= 1:
        diag_trans = np.zeros((nband, ngpt, nlev, ncol))
        diag_up_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))
        diag_dn_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))

    for b in range(nband):
        for g in range(ngpt):
            up_gpt = np.zeros((nlev + 1, ncol))
            down_gpt = np.zeros((nlev + 1, ncol))
            _lw_transport_single_gpt(
                tau[b, g, :, :],
                planck_source[b, g, :, :],
                surface_source[b, g, :],
                emissivity[b, :],
                up_gpt,
                down_gpt,
                nlev,
                ncol,
            )

            if diagnostics_level >= 1:
                for k in range(nlev):
                    for i in range(ncol):
                        diag_trans[b, g, k, i] = np.exp(
                            -DIFFUSIVITY_FACTOR * tau[b, g, k, i]
                        )
                diag_up_gpt[b, g, :, :] = up_gpt * weights[b, g]
                diag_dn_gpt[b, g, :, :] = down_gpt * weights[b, g]

            w = weights[b, g]
            for k in range(nlev + 1):
                for i in range(ncol):
                    up_band[b, k, i] += w * up_gpt[k, i]
                    down_band[b, k, i] += w * down_gpt[k, i]

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    if diagnostics_level >= 1:
        diag = {
            "transmittance": diag_trans,
            "up_per_gpoint": diag_up_gpt,
            "down_per_gpoint": diag_dn_gpt,
        }
        return up_band, down_band, up_broad, down_broad, diag

    return up_band, down_band, up_broad, down_broad

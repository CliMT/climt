# climt/_components/picket_fence/lw/kernels.py
import numpy as np

from ..common import njit, prange

DIFFUSIVITY_FACTOR = 1.66


@njit(parallel=True)
def planck_sources_kernel(
    planck_frac, T_grid, T, T_surf, sigma,
    nband, ngpt, is_esft, ngpt_orig, nband_orig,
    planck_src, surf_src,
):
    """Per-(band, g-point, level, column) Planck layer + surface sources.

    Linear-in-T interpolation of the band/g-point Planck fraction times sigma*T^4.
    Parallel over columns; each column writes disjoint planck_src[:,:,:,icol] and
    surf_src[:,:,icol]. Replicates the original pure-Python loop bit-for-bit
    (same searchsorted bracket, [0,1] fraction clamp, and float32 planck_frac
    promotion).
    """
    nT = T_grid.shape[0]
    nlev = T.shape[0]
    ncol = T.shape[1]
    for icol in prange(ncol):
        # --- Surface source ---
        T_s = T_surf[icol]
        iTs = np.searchsorted(T_grid, T_s) - 1
        if iTs < 0:
            iTs = 0
        elif iTs > nT - 2:
            iTs = nT - 2
        fTs = (T_s - T_grid[iTs]) / (T_grid[iTs + 1] - T_grid[iTs])
        if fTs < 0.0:
            fTs = 0.0
        elif fTs > 1.0:
            fTs = 1.0
        surf_planck = sigma * T_s ** 4
        for ib in range(nband):
            ib_orig = ib if ib < nband_orig else nband_orig - 1
            for igp in range(ngpt):
                g_idx_orig = igp % ngpt_orig if is_esft else igp
                frac_s = (planck_frac[ib_orig, g_idx_orig, iTs] * (1.0 - fTs)
                          + planck_frac[ib_orig, g_idx_orig, iTs + 1] * fTs)
                surf_src[ib, igp, icol] = frac_s * surf_planck

        # --- Layer source ---
        for kk in range(nlev):
            T_l = T[kk, icol]
            iTl = np.searchsorted(T_grid, T_l) - 1
            if iTl < 0:
                iTl = 0
            elif iTl > nT - 2:
                iTl = nT - 2
            fTl = (T_l - T_grid[iTl]) / (T_grid[iTl + 1] - T_grid[iTl])
            if fTl < 0.0:
                fTl = 0.0
            elif fTl > 1.0:
                fTl = 1.0
            layer_planck = sigma * T_l ** 4
            for ib in range(nband):
                ib_orig = ib if ib < nband_orig else nband_orig - 1
                for igp in range(ngpt):
                    g_idx_orig = igp % ngpt_orig if is_esft else igp
                    frac_l = (planck_frac[ib_orig, g_idx_orig, iTl] * (1.0 - fTl)
                              + planck_frac[ib_orig, g_idx_orig, iTl + 1] * fTl)
                    planck_src[ib, igp, kk, icol] = frac_l * layer_planck


@njit(parallel=True)
def _lw_transport_kernel(
    tau, planck_source, surface_source, emissivity, weights,
    up_band, down_band, up_broad, down_broad,
    diag_trans, diag_up_gpt, diag_dn_gpt, want_diag,
):
    """Consolidated multi-band, multi-g-point LW transport.

    Loops over columns in parallel; for each (band, g-point) runs the up/down
    diffusivity sweeps and accumulates weighted fluxes into up_band/down_band
    inside the compiled kernel. Accumulation order (g ascending, then b
    ascending for broadband) matches the original python loops bit-for-bit.
    """
    nband, ngpt, nlev, ncol = tau.shape
    for i in prange(ncol):
        for k in range(nlev + 1):
            up_broad[k, i] = 0.0
            down_broad[k, i] = 0.0
        for b in range(nband):
            for k in range(nlev + 1):
                up_band[b, k, i] = 0.0
                down_band[b, k, i] = 0.0
            for g in range(ngpt):
                w = weights[b, g]
                # Upward sweep: surface -> TOA
                up_prev = emissivity[b, i] * surface_source[b, g, i]
                up_band[b, 0, i] += w * up_prev
                if want_diag != 0:
                    diag_up_gpt[b, g, 0, i] = w * up_prev
                for k in range(nlev):
                    trans = np.exp(-DIFFUSIVITY_FACTOR * tau[b, g, k, i])
                    up_cur = up_prev * trans + planck_source[b, g, k, i] * (1.0 - trans)
                    up_band[b, k + 1, i] += w * up_cur
                    if want_diag != 0:
                        diag_trans[b, g, k, i] = trans
                        diag_up_gpt[b, g, k + 1, i] = w * up_cur
                    up_prev = up_cur
                # Downward sweep: TOA -> surface (dn_prev starts at 0 = TOA BC)
                dn_prev = 0.0
                if want_diag != 0:
                    diag_dn_gpt[b, g, nlev, i] = 0.0
                for k in range(nlev - 1, -1, -1):
                    trans = np.exp(-DIFFUSIVITY_FACTOR * tau[b, g, k, i])
                    dn_cur = dn_prev * trans + planck_source[b, g, k, i] * (1.0 - trans)
                    down_band[b, k, i] += w * dn_cur
                    if want_diag != 0:
                        diag_dn_gpt[b, g, k, i] = w * dn_cur
                    dn_prev = dn_cur
            for k in range(nlev + 1):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]


def lw_transport(
    T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma,
    diagnostics_level=0,
):
    """Multi-band, multi-g-point LW radiative transfer (consolidated kernel).

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
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))
    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))

    want_diag = 1 if diagnostics_level >= 1 else 0
    if want_diag:
        diag_trans = np.zeros((nband, ngpt, nlev, ncol))
        diag_up_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))
        diag_dn_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))
    else:
        diag_trans = np.zeros((1, 1, 1, 1))
        diag_up_gpt = np.zeros((1, 1, 1, 1))
        diag_dn_gpt = np.zeros((1, 1, 1, 1))

    _lw_transport_kernel(
        tau, planck_source, surface_source, emissivity, weights,
        up_band, down_band, up_broad, down_broad,
        diag_trans, diag_up_gpt, diag_dn_gpt, want_diag,
    )

    if want_diag:
        diag = {
            "transmittance": diag_trans,
            "up_per_gpoint": diag_up_gpt,
            "down_per_gpoint": diag_dn_gpt,
        }
        return up_band, down_band, up_broad, down_broad, diag

    return up_band, down_band, up_broad, down_broad

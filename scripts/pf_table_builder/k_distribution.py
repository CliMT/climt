"""Band-bin κ(ν) and quadrature to a small number of g-points."""
from __future__ import annotations

import numpy as np


def kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt):
    """Convert line-by-line κ to per-band k-distribution coefficients.

    Args:
        kappa: array shaped (..., nNu), arbitrary leading axes (typically
            (nT, nP) or (nT, nP, nX)).
        nu_grid: (nNu,) wavenumber grid (cm^-1), monotone increasing.
        band_edges: (nband+1,) band edges in cm^-1.
        ngpt: number of Gauss-Legendre nodes per band.

    Returns:
        k_coeffs: (..., nband, ngpt) k-coefficients evaluated at GL nodes.
        weights: (nband, ngpt) GL weights on [0, 1].
    """
    leading_shape = kappa.shape[:-1]
    nNu = kappa.shape[-1]
    assert len(nu_grid) == nNu, "kappa and nu_grid disagree on length"

    nband = len(band_edges) - 1
    # Gauss-Legendre on [-1, 1] -> [0, 1]
    xi_raw, wi_raw = np.polynomial.legendre.leggauss(ngpt)
    g_nodes = 0.5 * (xi_raw + 1.0)
    g_weights = 0.5 * wi_raw

    k_coeffs = np.zeros(leading_shape + (nband, ngpt))

    flat = kappa.reshape(-1, nNu)  # (nLead, nNu)

    for ib in range(nband):
        lo, hi = band_edges[ib], band_edges[ib + 1]
        mask = (nu_grid >= lo) & (nu_grid < hi)
        if not mask.any():
            raise ValueError(f"No nu_grid points in band [{lo}, {hi})")
        # Per-column: sort kappa in band, build CDF over wavenumber width,
        # then evaluate at GL g-nodes.
        for i in range(flat.shape[0]):
            k_band = np.sort(flat[i, mask])
            # Cumulative g axis: uniform in ν since the linepyline grid is uniform
            n_in = len(k_band)
            g_axis = (np.arange(n_in) + 0.5) / n_in
            # Interpolate k(g) onto the GL nodes
            k_coeffs.reshape(-1, nband, ngpt)[i, ib, :] = np.interp(
                g_nodes, g_axis, k_band,
            )

    weights = np.broadcast_to(g_weights, (nband, ngpt)).copy()
    return k_coeffs, weights

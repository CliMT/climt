"""Planck-source distribution across g-points within each LW band."""
import numpy as np

_H_PLANCK = 6.62607015e-34
_C_LIGHT = 2.99792458e8
_K_BOLTZ = 1.380649e-23
_SIGMA_SB = 5.670374419e-8


def _planck_B_wavenumber(wn_cm, T):
    nu_m = wn_cm * 100.0
    x = _H_PLANCK * _C_LIGHT * nu_m / (_K_BOLTZ * T)
    with np.errstate(over="ignore"):
        denom = np.expm1(np.clip(x, 0.0, 700.0))
    denom = np.where(denom == 0, np.inf, denom)
    B = 2.0 * _H_PLANCK * _C_LIGHT ** 2 * nu_m ** 3 / denom
    return B * 100.0  # per cm^-1


def _band_planck_fraction(wn_lo, wn_hi, T, n_quad=64):
    grid = np.linspace(wn_lo, wn_hi, n_quad + 1)
    mid = 0.5 * (grid[1:] + grid[:-1])
    dwn = np.diff(grid)
    band_flux = np.pi * np.sum(_planck_B_wavenumber(mid, T) * dwn)
    return band_flux / (_SIGMA_SB * T ** 4)


def build_uniform_planck_fraction(band_edges, ngpt, T_grid):
    """Per-band Planck integral, replicated across g-points within each band.

    The LW kernel evaluates ``surf_src = planck_fraction * sigma * T**4`` per
    (band, gpt) and then sums over g-points with ``gpoint_weights`` (which sum
    to 1 within a band). For that to recover the correct band emission, the
    stored ``planck_fraction[band, gpt, T]`` must equal the *band's* fraction
    of total blackbody emission at T — the same value across all g-points in
    the band, since the k-sort dimension has no spectral meaning.
    """
    band_edges = np.asarray(band_edges, dtype=np.float64)
    T_grid = np.asarray(T_grid, dtype=np.float64)
    nband = len(band_edges) - 1
    nT = len(T_grid)
    out = np.zeros((nband, ngpt, nT))
    for ib in range(nband):
        for iT, T in enumerate(T_grid):
            frac = _band_planck_fraction(band_edges[ib], band_edges[ib + 1], T)
            out[ib, :, iT] = frac
    return out

"""SW band-integrated solar source and Rayleigh cross-section."""
from __future__ import annotations

import numpy as np


def build_solar_source_per_gpoint(spectrum, band_edges, ngpt):
    """Distribute band-integrated stellar flux uniformly across g-points.

    Args:
        spectrum: dict with "wavenumber" (cm^-1) and "irradiance" (W/m^2/cm^-1).
        band_edges: (nband+1,) cm^-1.
        ngpt: per-band g-point count.

    Returns:
        (nband, ngpt) array, W/m^2 per (band, g-point).
    """
    wn = np.asarray(spectrum["wavenumber"])
    irr = np.asarray(spectrum["irradiance"])
    nband = len(band_edges) - 1
    out = np.zeros((nband, ngpt))
    try:
        _trapz = np.trapezoid
    except AttributeError:
        _trapz = np.trapz
    for ib in range(nband):
        lo, hi = band_edges[ib], band_edges[ib + 1]
        mask = (wn > lo) & (wn < hi)
        wn_b = np.concatenate(([lo], wn[mask], [hi]))
        irr_b = np.concatenate(([np.interp(lo, wn, irr)],
                                irr[mask],
                                [np.interp(hi, wn, irr)]))
        band_flux = _trapz(irr_b, wn_b)
        out[ib, :] = band_flux / ngpt
    return out


def build_rayleigh_per_band(band_edges, mean_molar_mass_g,
                            refractivity_298K=2.97e-4):
    """Band-mean Rayleigh mass scattering coefficient (m²/kg).

    Uses Bodhaine et al. (1999) form. `refractivity_298K` is (n-1) at STP
    for the bulk gas; defaults to dry-air. For N₂ use ~2.97e-4, for CO₂ ~4.5e-4,
    for H₂ ~1.4e-4.

    The wavenumber dependence is ν⁴ (i.e. λ⁻⁴). We evaluate at the band
    geometric mean for the coefficient.
    """
    # Bodhaine: σ(λ) = (24 π³ / N²λ⁴) × ((n²-1)/(n²+2))²
    # mass cross-section = σ × N_A / M
    band_mid_wn = np.sqrt(band_edges[:-1] * band_edges[1:])  # cm^-1
    lam_m = 1.0 / (band_mid_wn * 100.0)  # cm^-1 -> m
    N_loschmidt = 2.547e25  # molecules/m^3 at STP
    f = (refractivity_298K + 1.0) ** 2
    sigma = (24.0 * np.pi**3) / (N_loschmidt**2 * lam_m**4) * \
            ((f - 1) / (f + 2)) ** 2
    # Convert per-molecule cross-section (m²) to per-mass (m²/kg)
    M = mean_molar_mass_g * 1e-3  # kg/mol
    N_A = 6.022e23
    return sigma * N_A / M

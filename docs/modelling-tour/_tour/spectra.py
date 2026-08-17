"""Spectral diagnostics built on cork's per-band output.

The EC2213 notes stop at "the OLR is this integrated over all lambda and over
all theta". Cork performs that integral; these helpers un-do it, turning
``*_per_band`` diagnostics into the spectra and weighting functions the notes
describe analytically.

Wavenumbers are in cm^-1 throughout, so ``c`` is carried in cm/s and spectral
flux densities come out in W m^-2 / cm^-1 — directly comparable to published
IRIS/Nimbus spectra.
"""
import numpy as np

H_PLANCK = 6.62607015e-34    # J s
C_LIGHT = 2.99792458e10      # cm/s  (note: cm, to pair with nu in cm^-1)
K_BOLTZ = 1.380649e-23       # J/K

# 2 h c^2 nu^3 comes out in W cm^-2 sr^-1 / cm^-1; 1e4 converts to per m^-2,
# and pi integrates isotropic intensity into a flux.
_PLANCK_PREFACTOR = np.pi * 1e4 * 2.0 * H_PLANCK * C_LIGHT ** 2


def planck_flux(nu, T):
    """Blackbody spectral flux density pi*B, in W m^-2 / cm^-1.

    Args:
        nu: wavenumber(s) in cm^-1.
        T: temperature(s) in K.
    """
    nu = np.asarray(nu, dtype=float)
    return _PLANCK_PREFACTOR * nu ** 3 / np.expm1(
        H_PLANCK * C_LIGHT * nu / (K_BOLTZ * np.asarray(T, dtype=float)))


def brightness_temperature(flux_density, nu):
    """Temperature of a blackbody radiating ``flux_density`` at ``nu``.

    The exact inverse of :func:`planck_flux`.

    Args:
        flux_density: spectral flux density in W m^-2 / cm^-1.
        nu: wavenumber(s) in cm^-1.
    """
    nu = np.asarray(nu, dtype=float)
    flux_density = np.asarray(flux_density, dtype=float)
    return (H_PLANCK * C_LIGHT * nu / K_BOLTZ) / np.log1p(
        _PLANCK_PREFACTOR * nu ** 3 / flux_density)


def band_limits_of(component):
    """(nband, 2) band edges in cm^-1 for a correlated-k cork component."""
    return np.asarray(component._table["band_wavenumber_limits"], dtype=float)


def band_centres(band_limits):
    """(nband,) band-centre wavenumbers in cm^-1."""
    return 0.5 * (band_limits[:, 0] + band_limits[:, 1])


def band_widths(band_limits):
    """(nband,) band widths in cm^-1."""
    return band_limits[:, 1] - band_limits[:, 0]


def spectral_olr(diagnostics, band_limits, column=0):
    """Per-band OLR and spectral flux density at the top of the atmosphere.

    Args:
        diagnostics: the dict returned by a ``CorkLongwaveRadiation`` call.
        band_limits: (nband, 2) band edges, from :func:`band_limits_of`.
        column: horizontal index to extract.

    Returns:
        (olr_per_band, flux_density) — W m^-2 per band, and W m^-2 / cm^-1.
    """
    up = diagnostics["upwelling_longwave_flux_in_air_per_band"].values
    olr_band = np.asarray(up[-1, column, 0, :] if up.ndim == 4
                          else up[-1, column, :], dtype=float)
    return olr_band, olr_band / band_widths(band_limits)


def tau_star_cumulative(tau_layer, diffusivity_factor):
    """Diffusivity-scaled optical depth measured DOWN from space.

    Args:
        tau_layer: (nz,) per-layer optical depth, surface first.
        diffusivity_factor: D, so the result is D * sum of layers above.

    Returns:
        (nz,) tau* at each mid-level — largest at the surface, ~0 at the top.
    """
    tau_layer = np.asarray(tau_layer, dtype=float)
    # Cumulative sum of everything strictly above each level.
    above = np.cumsum(tau_layer[::-1])[::-1] - tau_layer
    return diffusivity_factor * (above + 0.5 * tau_layer)


def emission_weight(tau_layer, diffusivity_factor):
    """How much of the OLR each layer actually emits.

    The notes derive W ∝ tau* exp(-tau*), peaking at tau* = 1, by assuming
    thin layers and an exponential density profile. This is the exact discrete
    equivalent for the model's own grid: a layer emits ``1 - exp(-D tau_k)`` of
    a blackbody and that emission is attenuated by ``exp(-tau*_above)`` on the
    way out. It reduces to the notes' form as the layers get thin.

    Args:
        tau_layer: (nz,) per-layer optical depth, surface first.
        diffusivity_factor: D.

    Returns:
        (nz,) dimensionless weight per layer, surface first.
    """
    tau_layer = np.asarray(tau_layer, dtype=float)
    above = np.cumsum(tau_layer[::-1])[::-1] - tau_layer
    emitted = 1.0 - np.exp(-diffusivity_factor * tau_layer)
    transmitted = np.exp(-diffusivity_factor * above)
    return emitted * transmitted

import os

import importlib_resources
import numpy as np


def load_stellar_spectrum(name_or_path):
    """Load a stellar spectrum from a built-in name or file path.

    Args:
        name_or_path: "sun", "trappist1", or path to .npz file

    Returns:
        dict with keys "wavenumber" (cm^-1) and "irradiance" (W/m^2/cm^-1)
    """
    if os.path.isfile(name_or_path):
        data = np.load(name_or_path)
    else:
        data_path = importlib_resources.files(
            "climt._data.picket_fence.stellar_spectra"
        ).joinpath(f"{name_or_path}.npz")
        with importlib_resources.as_file(data_path) as f:
            data = np.load(f)

    return {
        "wavenumber": np.array(data["wavenumber"]),
        "irradiance": np.array(data["irradiance"]),
    }


def integrate_spectrum_over_bands(spectrum, band_wavenumber_limits):
    """Integrate a stellar spectrum over spectral bands.

    Args:
        spectrum: dict with "wavenumber" (cm^-1, ascending) and
                  "irradiance" (W/m^2/cm^-1)
        band_wavenumber_limits: (nband, 2) lower and upper wavenumber limits, cm^-1

    Returns:
        flux_per_band: (nband,) integrated flux per band, W/m^2
    """
    wn = spectrum["wavenumber"]
    irr = spectrum["irradiance"]
    nband = band_wavenumber_limits.shape[0]
    flux = np.zeros(nband)

    for b in range(nband):
        wn_lo, wn_hi = band_wavenumber_limits[b]

        mask = (wn > wn_lo) & (wn < wn_hi)
        wn_inner = wn[mask]
        irr_inner = irr[mask]

        irr_lo = np.interp(wn_lo, wn, irr)
        irr_hi = np.interp(wn_hi, wn, irr)

        wn_band = np.concatenate(([wn_lo], wn_inner, [wn_hi]))
        irr_band = np.concatenate(([irr_lo], irr_inner, [irr_hi]))

        flux[b] = np.trapezoid(irr_band, wn_band)

    return flux

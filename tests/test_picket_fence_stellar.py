import numpy as np
import pytest


def test_load_solar_spectrum():
    """Loading the solar spectrum returns arrays with matching shapes."""
    from climt._components.picket_fence.optics.stellar import load_stellar_spectrum

    spec = load_stellar_spectrum("sun")
    assert "wavenumber" in spec
    assert "irradiance" in spec
    assert len(spec["wavenumber"]) == len(spec["irradiance"])
    assert len(spec["wavenumber"]) > 5


def test_solar_spectrum_integrates_to_tsi():
    """The solar spectrum should integrate to approximately 1361 W/m^2."""
    from climt._components.picket_fence.optics.stellar import load_stellar_spectrum

    spec = load_stellar_spectrum("sun")
    total = np.trapezoid(spec["irradiance"], spec["wavenumber"])
    np.testing.assert_allclose(total, 1361.0, rtol=0.02)


def test_integrate_spectrum_over_bands():
    """Band integration distributes total flux across bands."""
    from climt._components.picket_fence.optics.stellar import (
        integrate_spectrum_over_bands,
        load_stellar_spectrum,
    )

    spec = load_stellar_spectrum("sun")
    # 3 equal-width bands spanning the spectrum
    band_limits = np.array([
        [100.0, 10000.0],
        [10000.0, 25000.0],
        [25000.0, 50000.0],
    ])
    flux_per_band = integrate_spectrum_over_bands(spec, band_limits)
    assert flux_per_band.shape == (3,)
    assert np.all(flux_per_band >= 0)
    # Sum should be close to TSI (bands cover entire range)
    np.testing.assert_allclose(flux_per_band.sum(), 1361.0, rtol=0.05)


def test_integrate_spectrum_narrow_band():
    """A very narrow band should get a small fraction of total flux."""
    from climt._components.picket_fence.optics.stellar import (
        integrate_spectrum_over_bands,
        load_stellar_spectrum,
    )

    spec = load_stellar_spectrum("sun")
    band_limits = np.array([[14000.0, 14100.0]])  # 100 cm^-1 narrow band
    flux = integrate_spectrum_over_bands(spec, band_limits)
    assert flux[0] > 0
    assert flux[0] < 100.0  # much less than TSI


def test_trappist1_spectrum_loads():
    from climt._components.picket_fence.optics.stellar import load_stellar_spectrum
    spec = load_stellar_spectrum("trappist1")
    wn = spec["wavenumber"]
    irr = spec["irradiance"]
    assert wn.ndim == 1 and irr.ndim == 1
    assert wn.size == irr.size
    assert wn.size > 100
    assert (irr >= 0).all()
    # TRAPPIST-1 total irradiance at 1 AU is about 4% of solar
    # (0.000553 L_sun, times (1 AU / 1 AU)^2 -- the SED is tabulated at 1 AU)
    total = np.trapezoid(irr, wn)
    assert 20.0 < total < 80.0, f"expected ~55 W/m^2 at 1 AU, got {total:.1f}"

import numpy as np


def test_solar_source_partitions_total_irradiance():
    """Sum of solar_source across (band, gpt) equals toa_irradiance exactly."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.solar_source import (
        build_solar_source_per_gpoint, build_rayleigh_per_band,
    )

    wn = np.linspace(1000.0, 30000.0, 5000)
    irr = np.ones_like(wn)
    spectrum = {"wavenumber": wn, "irradiance": irr}
    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])

    src = build_solar_source_per_gpoint(spectrum, band_edges, ngpt=2,
                                        toa_irradiance=878.0)
    assert src.shape == (3, 2)
    np.testing.assert_allclose(src.sum(), 878.0, rtol=1e-6)


def test_solar_source_spectral_shape_preserved():
    """Relative fractions per band should depend on spectrum shape, not normalization."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.solar_source import build_solar_source_per_gpoint

    # Spectrum with twice the power in the second half
    wn = np.linspace(3250.0, 30000.0, 2000)
    mid = 16625.0
    irr = np.where(wn < mid, 1.0, 2.0)
    spectrum = {"wavenumber": wn, "irradiance": irr}
    band_edges = np.array([3250.0, mid, 30000.0])

    src1 = build_solar_source_per_gpoint(spectrum, band_edges, ngpt=1,
                                         toa_irradiance=100.0)
    src2 = build_solar_source_per_gpoint(spectrum, band_edges, ngpt=1,
                                         toa_irradiance=200.0)  # double toa
    # Fractions per band must be the same regardless of toa_irradiance
    np.testing.assert_allclose(src1 / src1.sum(), src2 / src2.sum(), rtol=1e-9)
    # And total scales with toa_irradiance
    np.testing.assert_allclose(src2.sum() / src1.sum(), 2.0, rtol=1e-9)


def test_solar_source_spectrum_normalization_irrelevant():
    """toa_irradiance controls total flux regardless of spectrum's absolute scale."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.solar_source import build_solar_source_per_gpoint

    wn = np.linspace(3250.0, 30000.0, 1000)
    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])

    # Same spectral shape, different absolute scales
    spec_a = {"wavenumber": wn, "irradiance": np.ones_like(wn)}
    spec_b = {"wavenumber": wn, "irradiance": np.ones_like(wn) * 55.0}

    src_a = build_solar_source_per_gpoint(spec_a, band_edges, ngpt=2,
                                          toa_irradiance=878.0)
    src_b = build_solar_source_per_gpoint(spec_b, band_edges, ngpt=2,
                                          toa_irradiance=878.0)
    np.testing.assert_allclose(src_a, src_b, rtol=1e-9)


def test_rayleigh_per_band_decreases_with_wavelength():
    """Rayleigh coefficient should drop as wavenumber drops (longer λ)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.solar_source import build_rayleigh_per_band

    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])
    r = build_rayleigh_per_band(band_edges, mean_molar_mass_g=28.0,
                                refractivity_298K=2.97e-4)
    assert r.shape == (3,)
    assert r[0] < r[1] < r[2]
    assert (r >= 0).all()

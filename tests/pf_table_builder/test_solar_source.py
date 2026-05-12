import numpy as np


def test_solar_source_partitions_total_irradiance():
    """Sum of solar_source across (band, gpt) equals total stellar flux integral."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.solar_source import (
        build_solar_source_per_gpoint, build_rayleigh_per_band,
    )

    # Simple flat spectrum: 1 W/m²/cm⁻¹ from 1000 to 30000 cm^-1
    wn = np.linspace(1000.0, 30000.0, 5000)
    irr = np.ones_like(wn)
    spectrum = {"wavenumber": wn, "irradiance": irr}
    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])

    src = build_solar_source_per_gpoint(spectrum, band_edges, ngpt=2)
    assert src.shape == (3, 2)
    # Total = 30000 - 3250 (flat unit spectrum)
    np.testing.assert_allclose(src.sum(), 30000.0 - 3250.0, rtol=1e-3)


def test_rayleigh_per_band_decreases_with_wavelength():
    """Rayleigh coefficient should drop as wavenumber drops (longer λ)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.solar_source import build_rayleigh_per_band

    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])
    r = build_rayleigh_per_band(band_edges, mean_molar_mass_g=28.0,
                                refractivity_298K=2.97e-4)
    assert r.shape == (3,)
    assert r[0] < r[1] < r[2]
    assert (r >= 0).all()

import numpy as np


def test_uniform_planck_fraction_band_integral():
    """planck_fraction[band, gpt, T] must equal the band's blackbody share,
    constant across g-points within the band, so that summing
    sum_b sum_g gpoint_weights[b,g] * planck_fraction[b,g,T] * sigma*T^4
    recovers the band-resolved Planck emission. gpoint_weights sum to 1
    within a band, so cork is the same value at every g-point of the band."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.planck_fraction import build_uniform_planck_fraction

    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    ngpt = 2
    T_grid = np.array([200.0, 250.0, 300.0, 350.0, 400.0])
    cork = build_uniform_planck_fraction(band_edges, ngpt, T_grid)

    nband = len(band_edges) - 1
    assert cork.shape == (nband, ngpt, len(T_grid))
    # Same value across g-points within each band.
    np.testing.assert_allclose(cork[:, 0, :], cork[:, 1, :], rtol=1e-12)
    # Per-band positivity and ≤ 1.
    assert (cork >= 0).all()
    assert (cork <= 1.0 + 1e-12).all()
    # 10–3250 cm^-1 captures essentially all the LW Planck flux at these
    # temperatures, so the per-T sum across bands should be ~1.
    band_sums = cork[:, 0, :].sum(axis=0)
    np.testing.assert_allclose(band_sums, 1.0, atol=1e-2)

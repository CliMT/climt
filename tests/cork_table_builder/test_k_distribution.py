import numpy as np
import pytest


def test_kappa_to_k_coeffs_uniform_recovers_value():
    """A uniform kappa(nu) in a band should give k_coeffs = kappa at every g-point."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.k_distribution import kappa_to_k_coeffs

    nu_grid = np.linspace(10.0, 3250.0, 1000)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    # Constant kappa = 1e-3 everywhere
    kappa = 1e-3 * np.ones((1, 1, len(nu_grid)))  # (nT=1, nP=1, nNu)

    k_coeffs, weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt=2)
    # Shape: (nT, nP, nband, ngpt)
    assert k_coeffs.shape == (1, 1, 4, 2)
    np.testing.assert_allclose(k_coeffs, 1e-3, rtol=1e-10)
    # Gauss-Legendre weights on [0,1] sum to 1
    np.testing.assert_allclose(weights.sum(axis=-1), 1.0, rtol=1e-10)


def test_kappa_to_k_coeffs_two_peaks_orders_correctly():
    """Within a band, k_coeffs at successive g-points are monotone increasing."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.k_distribution import kappa_to_k_coeffs

    nu_grid = np.linspace(10.0, 3250.0, 5000)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    # Pretend lines: kappa large at every 50 cm^-1, small elsewhere
    kappa = 1e-5 * np.ones_like(nu_grid)
    kappa[::50] = 1.0
    kappa = kappa.reshape(1, 1, -1)

    k_coeffs, _ = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt=4)
    # In every band, k must be non-decreasing in g
    for ib in range(4):
        for ig in range(3):
            assert k_coeffs[0, 0, ib, ig] <= k_coeffs[0, 0, ib, ig + 1] + 1e-12


def test_kappa_to_k_coeffs_with_extra_axes():
    """Builder is shape-agnostic in the leading kappa axes."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.k_distribution import kappa_to_k_coeffs

    nu_grid = np.linspace(10.0, 3250.0, 200)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    kappa = np.random.RandomState(0).uniform(1e-6, 1e-2,
                                              size=(3, 4, 5, len(nu_grid)))
    k_coeffs, weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt=2)
    assert k_coeffs.shape == (3, 4, 5, 4, 2)
    assert weights.shape == (4, 2)

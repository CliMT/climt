# tests/test_picket_fence_optics.py
import numpy as np
import pytest


def test_parmentier_coefficients_grey_limit():
    """When gamma_P=1 and R=1, kappa_1 should equal kappa_2 (grey limit)."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_thermal_opacities,
    )

    kappa_R = 1e-3  # m^2/kg
    gamma_P = 1.0
    beta = 0.5
    R = 1.0

    kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)
    np.testing.assert_allclose(kappa_1, kappa_2, rtol=1e-10)
    np.testing.assert_allclose(kappa_1, kappa_R, rtol=1e-10)


def test_parmentier_coefficients_nongrey():
    """When R > 1, kappa_1 > kappa_R > kappa_2."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_thermal_opacities,
    )

    kappa_R = 1e-3
    gamma_P = 10.0
    beta = 0.8
    R = 100.0

    kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)
    assert kappa_1 > kappa_R
    assert kappa_2 < kappa_R
    assert kappa_1 / kappa_2 == pytest.approx(R, rel=1e-6)


def test_parmentier_rosseland_mean():
    """kappa_1 and kappa_2 should reconstruct kappa_R via the Rosseland mean formula."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_thermal_opacities,
    )

    kappa_R = 1e-2
    gamma_P = 5.0
    beta = 0.7
    R = 50.0

    kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)
    # Rosseland mean: 1/kappa_R = beta/kappa_1 + (1-beta)/kappa_2
    # => kappa_R = kappa_1 * kappa_2 / (beta*kappa_2 + (1-beta)*kappa_1)
    kappa_R_reconstructed = kappa_1 * kappa_2 / (beta * kappa_2 + (1 - beta) * kappa_1)
    np.testing.assert_allclose(kappa_R_reconstructed, kappa_R, rtol=1e-10)


def test_parmentier_lookup_coefficients():
    """Lookup of ratio coefficients from T_eff should return plausible values."""
    from climt._components.picket_fence.optics.parmentier import (
        load_parmentier_coefficients,
        lookup_ratio_coefficients,
    )

    coeffs = load_parmentier_coefficients("solar_composition")
    # T_eff = 1500 K is in the 600-1400 K range for solar composition
    gamma_v1, gamma_v2, gamma_v3, beta, gamma_P, R = lookup_ratio_coefficients(
        coeffs, 1500.0
    )

    assert gamma_v1 > 0
    assert gamma_v2 > 0
    assert gamma_v3 > 0
    assert 0 < beta <= 1
    assert gamma_P >= 1
    assert R >= 1

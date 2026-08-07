"""Surface-flux behaviour of climt.SimpleBoundaryLayer.

Covers the three ``surface_fluxes`` modes ('bulk', 'external', None), the
tridiagonal surface boundary condition that implements them, and the column
budgets they must close.
"""
from datetime import timedelta

import numpy as np
import pytest
from sympl import get_constant

import climt
from climt._components.simple_boundary_layer.component import _diffuse_profile


# ---------------------------------------------------------------- fixtures

def _test_column():
    """A small, well-conditioned 4-layer column for _diffuse_profile tests."""
    p_int = np.array([1.0e5, 9.0e4, 8.0e4, 7.0e4, 6.0e4])
    p = 0.5 * (p_int[:-1] + p_int[1:])
    rho = np.array([1.0, 0.9, 0.8])      # length num_layers - 1
    diff = np.array([10.0, 8.0, 6.0])    # length num_layers - 1
    return p, p_int, rho, diff


# -------------------------------------------- _diffuse_profile boundary terms

def test_diffuse_profile_zero_boundary_terms_reproduce_no_flux():
    """Both scalars zero must reproduce the pre-existing no-flux solve."""
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              0.0, 0.0)
    dp = p_int[:-1] - p_int[1:]
    assert np.isclose(np.sum(result * dp), np.sum(profile * dp), rtol=1e-12)


def test_diffuse_profile_does_not_mutate_profile():
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    before = profile.copy()
    _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665, 0.5, 140.0)
    assert np.array_equal(profile, before)


def test_diffuse_profile_bulk_term_matches_analytic_two_layer():
    """With diff == 0 the system is diagonal, so layer 0 has a closed form."""
    p, p_int, rho, _ = _test_column()
    diff = np.zeros(3)
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    beta = 0.25
    surface_value = 300.0
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              beta, beta * surface_value)
    expected0 = (profile[0] + beta * surface_value) / (1.0 + beta)
    assert np.isclose(result[0], expected0, rtol=1e-12)
    assert np.allclose(result[1:], profile[1:], rtol=1e-12)


def test_diffuse_profile_neumann_source_adds_exact_column_amount():
    """A pure Neumann source changes the mass-weighted integral by dp0 * S."""
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    source = 3.0
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              0.0, source)
    dp = p_int[:-1] - p_int[1:]
    change = np.sum(result * dp) - np.sum(profile * dp)
    assert np.isclose(change, dp[0] * source, rtol=1e-10)

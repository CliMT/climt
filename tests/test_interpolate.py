"""Tests for the pure-numpy not-a-knot cubic spline used to replace scipy.

scipy.interpolate.CubicSpline is used ONLY as a dev-environment oracle here
to verify bit-closeness; it is never imported by shipped climt code.
"""

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from climt._core.interpolate import cubic_spline_interpolate


def _random_strictly_increasing(rng, n):
    steps = rng.uniform(0.1, 2.0, size=n - 1)
    x = np.concatenate([[rng.uniform(-5, 5)], steps]).cumsum()
    return x


@pytest.mark.parametrize("n", [4, 5, 12, 30])
def test_matches_scipy_cubic_spline_including_extrapolation(n):
    rng = np.random.default_rng(1234 + n)
    for trial in range(5):
        x = _random_strictly_increasing(rng, n)
        y = rng.normal(size=n)

        span = x[-1] - x[0]
        x_new = np.concatenate(
            [
                [x[0] - 0.5 * span, x[0] - 0.05 * span],  # below range
                np.linspace(x[0], x[-1], 25),  # in range (includes knots)
                [x[-1] + 0.05 * span, x[-1] + 0.5 * span],  # above range
            ]
        )

        expected = CubicSpline(x, y)(x_new)
        actual = cubic_spline_interpolate(x_new, x, y)

        np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-12)


def test_matches_scipy_on_real_ozone_profile_with_extrapolation():
    import importlib.resources as importlib_resources

    data_path = importlib_resources.files("climt._data").joinpath(
        "ozone_profile.npy"
    )
    with importlib_resources.as_file(data_path) as f:
        ozone_ref = np.load(f)
    p_ref = 1e5 * np.linspace(0.998, 0.001, 30)

    # np.interp/spline both require increasing xp; p_ref is descending.
    x = p_ref[::-1]
    y = ozone_ref[::-1]

    assert x[0] == pytest.approx(100.0)
    assert x[-1] == pytest.approx(99800.0)

    # Deliberately include points above x[-1] (99800 Pa) and below x[0]
    # (100 Pa) to exercise scipy's default extrapolate=True behavior.
    p_new = np.linspace(1.01e5, 50.0, 28)

    expected = CubicSpline(x, y)(p_new)
    actual = cubic_spline_interpolate(p_new, x, y)

    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-12)

    max_abs_diff = np.max(np.abs(actual - expected))
    max_rel_diff = np.max(
        np.abs(actual - expected) / np.maximum(np.abs(expected), 1e-30)
    )
    assert max_abs_diff < 1e-12 or max_rel_diff < 1e-9

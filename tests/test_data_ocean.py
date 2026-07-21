import numpy as np
from datetime import datetime
from climt._components.data_ocean._sst_interpolation import (
    mid_month_values, interp_time,
)


def test_mid_month_preserves_monthly_means():
    rng = np.random.default_rng(1)
    means = 285.0 + rng.normal(0, 3, size=12)
    mm = mid_month_values(means)
    # Independent geometric reconstruction: month m spans from the midpoint of
    # (mm[m-1],mm[m]) to the midpoint of (mm[m],mm[m+1]); average the two-segment
    # piecewise-linear curve (boundary -> peak mm[m] -> boundary) over the month.
    recon = np.empty(12)
    for m in range(12):
        a = 0.5 * (mm[(m - 1) % 12] + mm[m])   # start-of-month boundary value
        b = 0.5 * (mm[m] + mm[(m + 1) % 12])   # end-of-month boundary value
        recon[m] = 0.5 * (0.5 * (a + mm[m]) + 0.5 * (mm[m] + b))
    assert np.allclose(recon, means, atol=1e-6)


def test_mid_month_preserves_means_numerically():
    rng = np.random.default_rng(2)
    means = 285.0 + rng.normal(0, 3, size=12)
    mm = mid_month_values(means)

    # Independent numerical-integration check: sample the piecewise-linear
    # interpolant through mid-month values (mid-months at fractions
    # (m+0.5)/12 of the year) on a fine uniform grid, and numerically
    # average it over each equal month-interval [m/12, (m+1)/12].
    n = 12
    mid_fracs = (np.arange(n) + 0.5) / n
    # Extend cyclically so interpolation covers the full [0, 1] range.
    ext_fracs = np.concatenate(
        [[mid_fracs[-1] - 1.0], mid_fracs, [mid_fracs[0] + 1.0]]
    )
    ext_vals = np.concatenate([[mm[-1]], mm, [mm[0]]])

    n_samples = 20000
    x = np.linspace(0.0, 1.0, n_samples, endpoint=False) + 0.5 / n_samples
    interp_vals = np.interp(x, ext_fracs, ext_vals)

    recon = np.empty(n)
    for m in range(n):
        lo, hi = m / n, (m + 1) / n
        mask = (x >= lo) & (x < hi)
        recon[m] = interp_vals[mask].mean()

    assert np.allclose(recon, means, atol=1e-2)


def test_interp_time_within_bounds():
    means = np.linspace(280, 300, 12)
    mm = mid_month_values(means)
    for month in range(1, 13):
        val = interp_time(mm, datetime(2000, month, 15, 12))
        assert mm.min() - 1e-6 <= val <= mm.max() + 1e-6

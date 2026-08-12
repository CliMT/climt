"""Mid-month (Taylor, Williamson & Zwiers 2000) SST time interpolation.

Monthly-mean SSTs are treated as the average, over each month, of a
piecewise-linear profile through mid-month boundary values. With mid-months
placed at the center of each (equal-length) month, month ``m`` runs from the
midpoint of ``(mm[m-1], mm[m])`` to the midpoint of ``(mm[m], mm[m+1])``.
Averaging that two-segment linear curve over the month gives the cyclic
tridiagonal relation

    mean_m = (1/8)*mm[(m-1) % n] + (3/4)*mm[m] + (1/8)*mm[(m+1) % n]

We solve this cyclic tridiagonal system for the mid-month values ``mm`` whose
linear reconstruction reproduces the input monthly means exactly. The
resulting matrix is strictly diagonally dominant (3/4 > 1/8 + 1/8) and hence
non-singular, unlike a naive bidiagonal averaging discretization
(``mean_m = 0.5*(mm[m] + mm[m+1])``), which is a singular circulant matrix
for n=12 (it has a zero eigenvalue at the alternating mode) and cannot be
solved directly.

This assumes equal-length months (idealized calendar). `interp_time` below
uses the true calendar day-of-year for interpolation, which is a minor,
accepted inconsistency (AMIP-grade simplification).
"""
import numpy as np

_DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], float)
_MID = np.cumsum(_DAYS) - _DAYS + _DAYS / 2.0   # day-of-year of each mid-month


def mid_month_values(monthly_means):
    """Solve for mid-month values from monthly means.

    Parameters
    ----------
    monthly_means : array_like
        Array of shape ``(n, ...)`` giving the monthly mean values (n is
        typically 12, but the solve is general for any n).

    Returns
    -------
    np.ndarray
        Array of shape ``(n, ...)`` of mid-month values whose piecewise-linear
        reconstruction reproduces ``monthly_means`` (see module docstring for
        the cyclic tridiagonal system that is solved).
    """
    means = np.asarray(monthly_means, dtype=float)
    n = means.shape[0]
    lead = means.shape[1:]
    flat = means.reshape(n, -1)

    # Cyclic tridiagonal: (1/8) on (m-1)%n, (3/4) on m, (1/8) on (m+1)%n.
    A = np.zeros((n, n))
    for m in range(n):
        A[m, (m - 1) % n] += 0.125
        A[m, m] += 0.75
        A[m, (m + 1) % n] += 0.125

    x = np.linalg.solve(A, flat)
    return x.reshape((n,) + lead)


def _day_of_year(dt):
    return (dt - type(dt)(dt.year, 1, 1)).days + \
        (dt.hour * 3600 + dt.minute * 60 + dt.second) / 86400.0


def interp_time(mid_month, model_time):
    doy = _day_of_year(model_time)
    mids = _MID
    # find bracketing mid-month indices on the cyclic year (365-day)
    if doy <= mids[0]:
        m0, m1 = 11, 0
        span = mids[0] + (365.0 - mids[11])
        frac = (doy + (365.0 - mids[11])) / span
    elif doy >= mids[-1]:
        m0, m1 = 11, 0
        span = mids[0] + (365.0 - mids[11])
        frac = (doy - mids[11]) / span
    else:
        m1 = int(np.searchsorted(mids, doy))
        m0 = m1 - 1
        frac = (doy - mids[m0]) / (mids[m1] - mids[m0])
    return (1 - frac) * mid_month[m0] + frac * mid_month[m1]

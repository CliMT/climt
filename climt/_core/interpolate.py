"""Pure-numpy not-a-knot cubic spline interpolation.

Reproduces ``scipy.interpolate.CubicSpline(x, y)(x_new)`` under scipy's
default settings (``bc_type="not-a-knot"``, ``extrapolate=True``) without
depending on scipy, so it can run in scipy-free environments such as
Pyodide/JupyterLite.

Only numpy and :func:`climt._core.tridiagonal.solve_tridiagonal` are used.

Construction
------------
The spline is built in terms of the second derivatives ``M_i = S''(x_i)``
at each breakpoint. Continuity of the first derivative at each interior
breakpoint gives the classic tridiagonal system (Numerical-Recipes form)

    h[i-1]*M[i-1] + 2*(h[i-1]+h[i])*M[i] + h[i]*M[i+1]
        = 6*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])

for interior knots ``i = 1 .. n-2``, where ``h[i] = x[i+1] - x[i]``. That
is ``n-2`` equations for ``n`` unknowns ``M[0..n-1]``; two more equations
come from the not-a-knot boundary condition, which requires the third
derivative to be continuous across the second and second-to-last knots
(equivalently: the first two polynomial pieces are literally the same
cubic, and likewise the last two):

    h[1]*M[0]     - (h[0]+h[1])*M[1]         + h[0]*M[2]     = 0
    -h[n-2]*M[n-3] + (h[n-3]+h[n-2])*M[n-2]  - h[n-3]*M[n-1] = 0

Each boundary equation is used to express the boundary unknown (``M[0]``
or ``M[n-1]``) as a linear combination of its two neighbors, which is then
substituted into the adjacent interior equation. This reduces the system
to a genuine tridiagonal system of size ``n-2`` in the interior unknowns
``M[1..n-2]``, solved with :func:`solve_tridiagonal`; ``M[0]`` and
``M[n-1]`` are then recovered by back-substitution into the (already
solved) boundary relations.
"""

import numpy as np

from .tridiagonal import solve_tridiagonal


def cubic_spline_interpolate(x_new, x, y):
    """Evaluate a not-a-knot cubic spline through ``(x, y)`` at ``x_new``.

    Parameters
    ----------
    x_new : array_like
        Points at which to evaluate the spline. May contain values outside
        ``[x[0], x[-1]]``; those are extrapolated using the cubic
        polynomial of the nearest end segment, matching scipy's default
        ``extrapolate=True`` behavior.
    x : array_like
        Strictly increasing breakpoints (caller guarantees this).
    y : array_like
        Values at ``x``.

    Returns
    -------
    numpy.ndarray
        Spline values at ``x_new``, shaped like ``x_new``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_new = np.asarray(x_new, dtype=float)

    n = x.size
    if n < 2:
        raise ValueError("cubic_spline_interpolate requires at least 2 points")

    h = np.diff(x)  # h[i] = x[i+1] - x[i], length n-1

    if n == 2:
        slope = (y[1] - y[0]) / h[0]
        return y[0] + slope * (x_new - x[0])

    m = _second_derivatives(x, y, h)
    return _evaluate_piecewise_cubic(x, y, h, m, x_new)


def _second_derivatives(x, y, h):
    """Solve for M[i] = S''(x[i]) under not-a-knot boundary conditions."""
    n = x.size

    if n == 3:
        # Only one interior knot exists, so the "second" and "second-to-
        # last" not-a-knot knots coincide and the two boundary equations
        # degenerate into the same equation. Fall back to a natural-like
        # closure (M[0] = M[2] = 0); climt never calls this with n == 3
        # (ozone init uses n == 30), so bit-exact agreement with scipy's
        # internal n == 3 special case is not required here.
        h0, h1 = h[0], h[1]
        a_diag = np.array([2.0 * (h0 + h1)])
        rhs = np.array(
            [6.0 * ((y[2] - y[1]) / h1 - (y[1] - y[0]) / h0)]
        )
        m1 = rhs[0] / a_diag[0]
        return np.array([0.0, m1, 0.0])

    # Interior rows for knots i = 1 .. n-2, reindexed locally as
    # k = i - 1 = 0 .. n-3 (size n-2 reduced system for M[1..n-2]).
    size = n - 2
    a_lower = np.zeros(size)
    a_diag = np.zeros(size)
    a_upper = np.zeros(size)
    rhs = np.zeros(size)

    i = np.arange(1, n - 1)
    a_lower[:] = h[i - 1]
    a_diag[:] = 2.0 * (h[i - 1] + h[i])
    a_upper[:] = h[i]
    rhs[:] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

    # Boundary row at x[1]: h1*M0 - (h0+h1)*M1 + h0*M2 = 0
    h0, h1 = h[0], h[1]
    c00, c01, c02 = h1, -(h0 + h1), h0

    # Eliminate M0 from the first interior row (k=0, i.e. i=1), whose
    # M0 coefficient is a_lower[0] = h0.
    factor0 = a_lower[0] / c00
    a_diag[0] -= factor0 * c01
    if size > 1:
        a_upper[0] -= factor0 * c02
    else:
        # size == 1 only happens for n == 3, already handled above.
        raise AssertionError("unreachable: n == 3 handled separately")
    a_lower[0] = 0.0
    rhs[0] -= factor0 * 0.0  # boundary rhs is 0

    # Boundary row at x[n-2]:
    #   -h[n-2]*M[n-3] + (h[n-3]+h[n-2])*M[n-2] - h[n-3]*M[n-1] = 0
    hnm3, hnm2 = h[n - 3], h[n - 2]
    cN0, cN1, cN2 = -hnm2, (hnm3 + hnm2), -hnm3

    # Eliminate M[n-1] from the last interior row (k=size-1, i.e.
    # i=n-2), whose M[n-1] coefficient is a_upper[size-1] = h[n-2].
    factor1 = a_upper[size - 1] / cN2
    a_diag[size - 1] -= factor1 * cN1
    if size > 1:
        a_lower[size - 1] -= factor1 * cN0
    a_upper[size - 1] = 0.0
    rhs[size - 1] -= factor1 * 0.0

    m_interior = solve_tridiagonal(a_lower, a_diag, a_upper, rhs)

    m = np.zeros(n)
    m[1 : n - 1] = m_interior
    m[0] = -(c01 * m[1] + c02 * m[2]) / c00
    m[n - 1] = -(cN0 * m[n - 3] + cN1 * m[n - 2]) / cN2
    return m


def _evaluate_piecewise_cubic(x, y, h, m, x_new):
    n = x.size
    # Locate the containing segment for each x_new; clamp so that points
    # outside [x[0], x[-1]] are evaluated with the first/last segment's
    # cubic polynomial (matches scipy's default extrapolate=True).
    idx = np.searchsorted(x, x_new, side="right") - 1
    idx = np.clip(idx, 0, n - 2)

    xi = x[idx]
    xi1 = x[idx + 1]
    hi = h[idx]
    mi = m[idx]
    mi1 = m[idx + 1]
    yi = y[idx]
    yi1 = y[idx + 1]

    a = xi1 - x_new
    b = x_new - xi

    return (
        (mi * a**3 + mi1 * b**3) / (6.0 * hi)
        + (yi / hi - mi * hi / 6.0) * a
        + (yi1 / hi - mi1 * hi / 6.0) * b
    )

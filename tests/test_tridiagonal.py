import numpy as np
from climt._core.tridiagonal import solve_tridiagonal


def _dense(a_lower, a_diag, a_upper):
    n = a_diag.size
    A = np.diag(a_diag).astype(float)
    for i in range(1, n):
        A[i, i - 1] = a_lower[i]
    for i in range(n - 1):
        A[i, i + 1] = a_upper[i]
    return A


def test_matches_numpy_solve_on_random_systems():
    rng = np.random.default_rng(0)
    for n in (2, 5, 20):
        a_lower = np.zeros(n); a_lower[1:] = rng.normal(size=n - 1)
        a_upper = np.zeros(n); a_upper[:-1] = rng.normal(size=n - 1)
        a_diag = rng.normal(size=n) + 5.0  # diagonally dominant
        rhs = rng.normal(size=n)
        x = solve_tridiagonal(a_lower, a_diag, a_upper, rhs)
        np.testing.assert_allclose(x, np.linalg.solve(_dense(a_lower, a_diag, a_upper), rhs), rtol=1e-10)

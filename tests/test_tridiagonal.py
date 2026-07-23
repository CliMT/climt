import numpy as np
import pytest

from climt._core.tridiagonal import solve_tridiagonal, tridiagonal_matvec


def _dense(lower, diag, upper):
    n = diag.shape[0]
    A = np.diag(diag).astype(float)
    for i in range(n - 1):
        A[i + 1, i] = lower[i]
        A[i, i + 1] = upper[i]
    return A


@pytest.mark.parametrize("n", [1, 2, 3, 5, 30])
def test_solve_matches_numpy(n):
    rng = np.random.RandomState(0)
    diag = rng.rand(n) + 3.0            # diagonally dominant => well-conditioned
    lower = rng.rand(max(n - 1, 0))
    upper = rng.rand(max(n - 1, 0))
    rhs = rng.rand(n)
    x = solve_tridiagonal(lower, diag, upper, rhs)
    if n == 1:
        expected = rhs / diag
    else:
        expected = np.linalg.solve(_dense(lower, diag, upper), rhs)
    assert np.allclose(x, expected, atol=1e-10)


def test_solve_does_not_mutate_inputs():
    lower = np.array([1.0, 1.0])
    diag = np.array([4.0, 4.0, 4.0])
    upper = np.array([1.0, 1.0])
    rhs = np.array([1.0, 2.0, 3.0])
    lower0, diag0, upper0, rhs0 = (a.copy() for a in (lower, diag, upper, rhs))
    solve_tridiagonal(lower, diag, upper, rhs)
    for a, a0 in ((lower, lower0), (diag, diag0), (upper, upper0), (rhs, rhs0)):
        assert np.array_equal(a, a0)


def test_matvec_matches_dense():
    rng = np.random.RandomState(1)
    n = 6
    diag = rng.rand(n)
    lower = rng.rand(n - 1)
    upper = rng.rand(n - 1)
    x = rng.rand(n)
    y = tridiagonal_matvec(lower, diag, upper, x)
    assert np.allclose(y, _dense(lower, diag, upper) @ x)


def test_solve_matvec_roundtrip():
    rng = np.random.RandomState(2)
    n = 10
    diag = rng.rand(n) + 3.0
    lower = rng.rand(n - 1)
    upper = rng.rand(n - 1)
    x = rng.rand(n)
    rhs = tridiagonal_matvec(lower, diag, upper, x)
    assert np.allclose(solve_tridiagonal(lower, diag, upper, rhs), x, atol=1e-10)

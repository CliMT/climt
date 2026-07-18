import numpy as np

from .backend import jit_compile


@jit_compile(backend=np)
def solve_tridiagonal(a_lower, a_diag, a_upper, rhs):
    """Solve the tridiagonal system A x = rhs via the Thomas algorithm.

    All band arrays have length n. a_lower[i] is A[i, i-1] (a_lower[0] unused);
    a_upper[i] is A[i, i+1] (a_upper[n-1] unused); a_diag[i] is A[i, i].
    """
    n = rhs.shape[0]
    cp = np.zeros(n)
    dp = np.zeros(n)
    x = np.zeros(n)
    cp[0] = a_upper[0] / a_diag[0]
    dp[0] = rhs[0] / a_diag[0]
    for i in range(1, n):
        m = a_diag[i] - a_lower[i] * cp[i - 1]
        cp[i] = a_upper[i] / m
        dp[i] = (rhs[i] - a_lower[i] * dp[i - 1]) / m
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x

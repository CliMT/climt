# -*- coding: utf-8 -*-
"""Canonical dense tridiagonal solver for climt.

This is THE dense tridiagonal solver for the library. Diffusion-type
components (boundary layer, snow/ice/soil columns, subsurface transport)
should call ``solve_tridiagonal`` rather than hand-rolling a Thomas sweep or
pulling in ``scipy.sparse``/``spsolve`` (which cannot be numba-jitted or, in
future, JAX-traced).

Tridiagonal families in climt, for reference:
  * dense, non-cyclic  -> ``solve_tridiagonal`` (this module) -- default.
  * cyclic/periodic    -> ``data_ocean._sst_interpolation`` (12x12
    ``np.linalg.solve``, cold path); needs a Sherman-Morrison variant if it
    ever moves onto a jit path.

JAX note: the forward/backward sweeps are loop-carried. A JAX backend would
express them via ``jax.lax.scan`` behind the ``jit_compile(backend=...)``
hook in ``backend.py``. Both functions here are pure (they allocate and
return; inputs are never mutated), which keeps that swap mechanical.
"""
import numpy as np

from .backend import jit_compile


@jit_compile
def solve_tridiagonal(lower, diag, upper, rhs):
    """Solve ``A x = rhs`` for a tridiagonal ``A`` via the Thomas algorithm.

    ``lower`` (sub-diagonal, ``A[i, i-1]``) and ``upper`` (super-diagonal,
    ``A[i, i+1]``) have length ``n-1``; ``diag`` and ``rhs`` have length ``n``.
    Returns the length-``n`` solution. Inputs are not mutated.
    """
    n = rhs.shape[0]
    x = np.zeros(n)
    if n == 1:
        x[0] = rhs[0] / diag[0]
        return x

    c_prime = np.zeros(n - 1)
    d_prime = np.zeros(n)
    c_prime[0] = upper[0] / diag[0]
    d_prime[0] = rhs[0] / diag[0]
    for i in range(1, n - 1):
        m = diag[i] - lower[i - 1] * c_prime[i - 1]
        c_prime[i] = upper[i] / m
        d_prime[i] = (rhs[i] - lower[i - 1] * d_prime[i - 1]) / m
    m = diag[n - 1] - lower[n - 2] * c_prime[n - 2]
    d_prime[n - 1] = (rhs[n - 1] - lower[n - 2] * d_prime[n - 2]) / m

    x[n - 1] = d_prime[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x


@jit_compile
def tridiagonal_matvec(lower, diag, upper, x):
    """Return ``A @ x`` for the banded layout used by ``solve_tridiagonal``.

    Same argument layout as ``solve_tridiagonal``. Pure.
    """
    n = x.shape[0]
    y = diag * x
    if n > 1:
        y[1:] = y[1:] + lower * x[:-1]
        y[:-1] = y[:-1] + upper * x[1:]
    return y

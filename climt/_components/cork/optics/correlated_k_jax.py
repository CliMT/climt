"""Differentiable JAX correlated-k optics (7-D additive-CO2 path).

Mirrors the Numba oracle ``_ck_tau_additive_co2_kernel`` in
``correlated_k.py``. Integer bracket indices are gathered from the static
k-table; gradients flow through the fractional interpolation weights, which
are smooth functions of the atmospheric state.
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

FLOOR = 1e-40


def _bracket(grid, v):
    """Index + fraction for multilinear interpolation. Matches _ck_bracket."""
    n = grid.shape[0]
    idx = jnp.clip(jnp.searchsorted(grid, v) - 1, 0, n - 2)
    lo = grid[idx]
    frac = jnp.clip((v - lo) / (grid[idx + 1] - lo), 0.0, 1.0)
    return idx, frac


def _txx7(k, iT, fT, iP, fP, iX, fX, iC):
    """8-corner trilinear over (T,P,X) at fixed C index. Returns
    (ngas,nband,ngpt,nlev,ncol)."""
    def g(a, b, c):
        return k[:, :, :, iT + a, iP + b, iX + c, iC]
    x0 = (g(0, 0, 0) * (1 - fT) * (1 - fP) + g(1, 0, 0) * fT * (1 - fP)
          + g(0, 1, 0) * (1 - fT) * fP + g(1, 1, 0) * fT * fP)
    x1 = (g(0, 0, 1) * (1 - fT) * (1 - fP) + g(1, 0, 1) * fT * (1 - fP)
          + g(0, 1, 1) * (1 - fT) * fP + g(1, 1, 1) * fT * fP)
    return x0 * (1 - fX) + x1 * fX


def _txx_cont(log_cont, iT, fT, iP, fP, iX, fX):
    """8-corner trilinear over (T,P,X) of log continuum. Returns (nband,nlev,ncol)."""
    def g(a, b, c):
        return log_cont[:, iT + a, iP + b, iX + c]
    x0 = (g(0, 0, 0) * (1 - fT) * (1 - fP) + g(1, 0, 0) * fT * (1 - fP)
          + g(0, 1, 0) * (1 - fT) * fP + g(1, 1, 0) * fT * fP)
    x1 = (g(0, 0, 1) * (1 - fT) * (1 - fP) + g(1, 0, 1) * fT * (1 - fP)
          + g(0, 1, 1) * (1 - fT) * fP + g(1, 1, 1) * fT * fP)
    return x0 * (1 - fX) + x1 * fX


def compute_tau_jax(T, log_p, log_x, log_c, gas_amounts,
                    k, T_grid, p_grid_log, log_x_grid, log_c_grid,
                    has_cont, log_cont, co2_logk):
    iT, fT = _bracket(T_grid, T)
    iP, fP = _bracket(p_grid_log, log_p)
    iX, fX = _bracket(log_x_grid, log_x)
    iC, fC = _bracket(log_c_grid, log_c)

    c0 = _txx7(k, iT, fT, iP, fP, iX, fX, iC)       # (ngas,nband,ngpt,nlev,ncol)
    c1 = _txx7(k, iT, fT, iP, fP, iX, fX, iC + 1)
    if co2_logk:
        l0 = jnp.log(jnp.maximum(c0, FLOOR))
        l1 = jnp.log(jnp.maximum(c1, FLOOR))
        kv = jnp.exp(l0 * (1 - fC) + l1 * fC)
    else:
        kv = c0 * (1 - fC) + c1 * fC

    # acc[ib,igp,nlev,ncol] = sum over gas of kv * gas_amount
    tau = jnp.sum(kv * gas_amounts[:, None, None, :, :], axis=0)

    if has_cont:
        cont_val = jnp.exp(_txx_cont(log_cont, iT, fT, iP, fP, iX, fX))  # (nband,nlev,ncol)
        tau = tau + cont_val[:, None, :, :] * gas_amounts[0][None, None, :, :]
    return tau

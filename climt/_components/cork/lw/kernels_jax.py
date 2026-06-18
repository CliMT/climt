"""Differentiable JAX longwave kernels (additive correlated-k path).

Mirrors the Numba oracles ``planck_sources_kernel`` and
``_lw_transport_kernel``. Valid for the non-ESFT additive case where the
g-point grid is not expanded (nband==nband_orig, ngpt==ngpt_orig).
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax

from ..optics.correlated_k_jax import _bracket

DIFFUSIVITY_FACTOR = 1.66


def planck_sources_jax(planck_frac, T_grid, T, T_surf, sigma):
    iT, fT = _bracket(T_grid, T)            # (nlev,ncol)
    pf0 = planck_frac[:, :, iT]             # (nband,ngpt,nlev,ncol)
    pf1 = planck_frac[:, :, iT + 1]
    frac_l = pf0 * (1.0 - fT) + pf1 * fT
    planck_src = frac_l * (sigma * T ** 4)

    iTs, fTs = _bracket(T_grid, T_surf)     # (ncol,)
    pfs0 = planck_frac[:, :, iTs]           # (nband,ngpt,ncol)
    pfs1 = planck_frac[:, :, iTs + 1]
    frac_s = pfs0 * (1.0 - fTs) + pfs1 * fTs
    surf_src = frac_s * (sigma * T_surf ** 4)
    return planck_src, surf_src

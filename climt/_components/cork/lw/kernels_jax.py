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


def lw_transport_jax(tau, planck_source, surface_source, emissivity, weights):
    """Differentiable LW radiative-transfer solver using ``jax.lax.scan``.

    Mirrors the Numba oracle ``lw_transport`` / ``_lw_transport_kernel`` for
    the forward (flux) computation. Supports jit/vmap/grad.

    Args:
        tau: (nband, ngpt, nlev, ncol) optical depth per layer
        planck_source: (nband, ngpt, nlev, ncol) Planck layer source, W/m^2
        surface_source: (nband, ngpt, ncol) surface Planck source, W/m^2
        emissivity: (nband, ncol) surface emissivity per band
        weights: (nband, ngpt) g-point quadrature weights

    Returns:
        up_band: (nband, nlev+1, ncol) weighted upward flux per band
        down_band: (nband, nlev+1, ncol) weighted downward flux per band
        up_broad: (nlev+1, ncol) broadband upward flux
        down_broad: (nlev+1, ncol) broadband downward flux
    """
    nband, ngpt, nlev, ncol = tau.shape
    D = DIFFUSIVITY_FACTOR
    trans = jnp.exp(-D * tau)                       # (nband,ngpt,nlev,ncol)
    trans_lev = jnp.moveaxis(trans, 2, 0)           # (nlev,nband,ngpt,ncol)
    planck_lev = jnp.moveaxis(planck_source, 2, 0)  # (nlev,nband,ngpt,ncol)

    def step(prev, layer):
        tr, pl = layer
        cur = prev * tr + pl * (1.0 - tr)
        return cur, cur

    # Upward sweep: surface -> TOA. BC = emissivity * surface_source.
    up0 = emissivity[:, None, :] * surface_source   # (nband,ngpt,ncol)
    _, up_stack = lax.scan(step, up0, (trans_lev, planck_lev))   # (nlev,nband,ngpt,ncol)
    up_iface = jnp.concatenate([up0[None], up_stack], axis=0)    # (nlev+1,nband,ngpt,ncol)

    # Downward sweep: TOA -> surface. BC = 0 at TOA.
    dn0 = jnp.zeros((nband, ngpt, ncol))
    _, dn_stack = lax.scan(step, dn0, (trans_lev, planck_lev), reverse=True)
    down_iface = jnp.concatenate([dn_stack, dn0[None]], axis=0)  # (nlev+1,nband,ngpt,ncol)

    up_g = jnp.moveaxis(up_iface, 0, 2)             # (nband,ngpt,nlev+1,ncol)
    down_g = jnp.moveaxis(down_iface, 0, 2)
    w = weights[:, :, None, None]
    up_band = jnp.sum(w * up_g, axis=1)             # (nband,nlev+1,ncol)
    down_band = jnp.sum(w * down_g, axis=1)
    up_broad = jnp.sum(up_band, axis=0)             # (nlev+1,ncol)
    down_broad = jnp.sum(down_band, axis=0)
    return up_band, down_band, up_broad, down_broad

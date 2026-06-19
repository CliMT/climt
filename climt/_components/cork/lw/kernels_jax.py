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
from ..optics.correlated_k import _CO2_INTERP_LOGK

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


# ---------------------------------------------------------------------------
# Phase 3 — static table container + differentiable orchestrator
# ---------------------------------------------------------------------------
from typing import NamedTuple
import numpy as np

from ..optics.correlated_k_jax import compute_tau_jax
from ..common import MOLAR_MASS, MOLAR_MASS_DRY_AIR

# Mass ratio M_H2O / M_dry_air, for specific humidity -> H2O VMR conversion.
_M_H2O_RATIO = MOLAR_MASS["h2o"] / MOLAR_MASS_DRY_AIR


class CorkLWTable(NamedTuple):
    k: jnp.ndarray
    T_grid: jnp.ndarray
    p_grid_log: jnp.ndarray
    log_x_grid: jnp.ndarray
    log_c_grid: jnp.ndarray
    planck_frac: jnp.ndarray
    gpoint_weights: jnp.ndarray
    has_cont: bool
    log_cont: jnp.ndarray
    co2_logk: bool


def build_cork_lw_table(table):
    """Host-side: pack a loaded numpy k-table dict into static jnp arrays.

    Only valid for the 7-D additive-CO2 premixed-background case.
    """
    k = np.asarray(table["k_coefficients"])
    assert k.ndim == 7, "JAX LW path requires a 7-D additive-CO2 k-table"
    x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
    c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
    has_cont = ("continuum_kappa" in table
                and np.asarray(table["continuum_kappa"]).ndim == 4)
    if has_cont:
        log_cont = np.log(np.maximum(
            np.asarray(table["continuum_kappa"], dtype=np.float64), 1e-40))
    else:
        log_cont = np.zeros((k.shape[1], 1, 1, 1), dtype=np.float64)
    return CorkLWTable(
        k=jnp.asarray(k, dtype=jnp.float64),
        T_grid=jnp.asarray(table["temperature_grid"], dtype=jnp.float64),
        p_grid_log=jnp.asarray(table["pressure_grid_log"], dtype=jnp.float64),
        log_x_grid=jnp.log(jnp.maximum(jnp.asarray(x_grid), 1e-30)),
        log_c_grid=jnp.log(jnp.maximum(jnp.asarray(c_grid), 1e-30)),
        planck_frac=jnp.asarray(table["planck_fraction"], dtype=jnp.float64),
        gpoint_weights=jnp.asarray(table["gpoint_weights"], dtype=jnp.float64),
        has_cont=bool(has_cont),
        log_cont=jnp.asarray(log_cont),
        co2_logk=bool(_CO2_INTERP_LOGK),
    )


def _column_amount_jax(q, p_int, g):
    dp = jnp.abs(p_int[1:] - p_int[:-1])
    return q * dp / g


def _heating_rate_jax(net_flux, p_int, g, cpd):
    dp = p_int[1:] - p_int[:-1]
    dflux = net_flux[1:] - net_flux[:-1]
    return g / cpd * dflux / dp


def cork_lw_jax(T, p, p_int, T_surf, h2o, co2, emissivity, tau_cloud,
                jtable, sigma, g, cpd):
    """Differentiable LW orchestrator: optics -> Planck -> transport -> heating.

    Args:
        T: (nlev, ncol) layer temperature [K]
        p: (nlev, ncol) layer pressure [Pa]
        p_int: (nlev+1, ncol) interface pressure [Pa]
        T_surf: (ncol,) surface temperature [K]
        h2o: (nlev, ncol) specific humidity [kg/kg]
        co2: (nlev, ncol) CO2 volume mixing ratio [mol/mol]
        emissivity: (nband, ncol) surface emissivity
        tau_cloud: (nlev, ncol, nband) cloud optical depth
        jtable: CorkLWTable (static, not traced)
        sigma: Stefan-Boltzmann constant [W/m^2/K^4]
        g: gravitational acceleration [m/s^2]
        cpd: dry-air heat capacity [J/kg/K]

    Returns:
        tendency: (nlev, ncol) temperature tendency [K/s]
        up_broad: (nlev+1, ncol) broadband upward flux [W/m^2]
        down_broad: (nlev+1, ncol) broadband downward flux [W/m^2]
        up_band: (nband, nlev+1, ncol) per-band upward flux [W/m^2]
        down_band: (nband, nlev+1, ncol) per-band downward flux [W/m^2]
        tau: (nband, ngpt, nlev, ncol) total optical depth
    """
    nlev, ncol = T.shape
    # gas amount is column air mass (premixed-bg); H2O enters via the X axis.
    air = _column_amount_jax(jnp.ones((nlev, ncol)), p_int, g)
    gas_amounts = air[None, :, :]                       # (1,nlev,ncol)

    h2o_vmr = h2o / jnp.maximum(h2o + (1.0 - h2o) * _M_H2O_RATIO, 1e-30)
    x_lo = jnp.exp(jtable.log_x_grid[0]); x_hi = jnp.exp(jtable.log_x_grid[-1])
    c_lo = jnp.exp(jtable.log_c_grid[0]); c_hi = jnp.exp(jtable.log_c_grid[-1])
    log_p = jnp.log(jnp.maximum(p, 1.0))
    log_x = jnp.log(jnp.maximum(jnp.clip(h2o_vmr, x_lo, x_hi), 1e-30))
    log_c = jnp.log(jnp.maximum(jnp.clip(co2, c_lo, c_hi), 1e-30))

    tau = compute_tau_jax(
        T, log_p, log_x, log_c, gas_amounts,
        jtable.k, jtable.T_grid, jtable.p_grid_log,
        jtable.log_x_grid, jtable.log_c_grid,
        jtable.has_cont, jtable.log_cont, jtable.co2_logk)

    # cloud optical depth: (nlev,ncol,nband) -> (nband,1,nlev,ncol)
    tau = tau + jnp.transpose(tau_cloud, (2, 0, 1))[:, None, :, :]

    planck_src, surf_src = planck_sources_jax(
        jtable.planck_frac, jtable.T_grid, T, T_surf, sigma)

    up_band, down_band, up_broad, down_broad = lw_transport_jax(
        tau, planck_src, surf_src, emissivity, jtable.gpoint_weights)

    net_flux = up_broad - down_broad
    tendency = _heating_rate_jax(net_flux, p_int, g, cpd)
    return tendency, up_broad, down_broad, up_band, down_band, tau

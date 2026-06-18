# tests/test_cork_lw_kernels_jax_parity.py
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest


def test_planck_sources_jax_matches_oracle():
    from climt._components.cork.lw.kernels import planck_sources_kernel
    from climt._components.cork.lw.kernels_jax import planck_sources_jax
    rng = np.random.RandomState(0)
    nband, ngpt, nT, nlev, ncol = 4, 3, 6, 7, 5
    planck_frac = rng.uniform(0.0, 1.0, size=(nband, ngpt, nT))
    T_grid = np.linspace(180.0, 330.0, nT)
    T = rng.uniform(200.0, 315.0, size=(nlev, ncol))
    T_surf = rng.uniform(260.0, 305.0, size=(ncol,))
    sigma = 5.670374419e-8

    planck_src = np.zeros((nband, ngpt, nlev, ncol))
    surf_src = np.zeros((nband, ngpt, ncol))
    planck_sources_kernel(
        np.ascontiguousarray(planck_frac),
        np.ascontiguousarray(T_grid, dtype=np.float64),
        np.ascontiguousarray(T, dtype=np.float64),
        np.ascontiguousarray(T_surf, dtype=np.float64),
        float(sigma), nband, ngpt, False, ngpt, nband,
        planck_src, surf_src,
    )
    ps_jax, ss_jax = planck_sources_jax(
        jnp.asarray(planck_frac), jnp.asarray(T_grid),
        jnp.asarray(T), jnp.asarray(T_surf), sigma)
    np.testing.assert_allclose(np.asarray(ps_jax), planck_src, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(np.asarray(ss_jax), surf_src, rtol=1e-6, atol=1e-9)

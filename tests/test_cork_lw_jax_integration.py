import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest
from climt._components.cork.optics.correlated_k import load_k_table


def test_cork_lw_table_roundtrip_shapes():
    from climt._components.cork.lw.kernels_jax import build_cork_lw_table
    table = load_k_table("earth_low_res_lw")
    jt = build_cork_lw_table(table)
    assert jt.k.shape == table["k_coefficients"].shape
    assert jt.planck_frac.shape == np.asarray(table["planck_fraction"]).shape
    assert jt.gpoint_weights.shape == np.asarray(table["gpoint_weights"]).shape
    assert jt.has_cont is True and jt.co2_logk is True

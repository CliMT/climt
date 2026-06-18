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


def _make_state(nlev, ncol, as_jax):
    from climt import get_grid, get_default_state
    from climt._components.cork import CorkLongwaveRadiation
    from climt._core.util import numpy_version_of
    comp = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    state = get_default_state([comp], grid_state=grid)
    raw = numpy_version_of(state)
    nband = comp.num_longwave_bands
    xp = jnp if as_jax else np

    def arr(name, shape):
        return xp.asarray(raw[name].reshape(shape))

    s = {
        "T": arr("air_temperature", (nlev, ncol)),
        "p": arr("air_pressure", (nlev, ncol)),
        "p_int": arr("air_pressure_on_interface_levels", (nlev + 1, ncol)),
        "T_surf": arr("surface_temperature", (ncol,)),
        "emissivity": arr("surface_longwave_emissivity", (nband, ncol)),
        "h2o": arr("specific_humidity", (nlev, ncol)),
        "co2": arr("mole_fraction_of_carbon_dioxide_in_air", (nlev, ncol)),
        "tau_cloud_lw": arr("longwave_optical_thickness_due_to_cloud", (nlev, ncol, nband)),
    }
    return comp, s


def test_array_call_jax_matches_numpy():
    nlev, ncol = 20, 4
    comp, s_np = _make_state(nlev, ncol, as_jax=False)
    _, s_jx = _make_state(nlev, ncol, as_jax=True)
    tend_np, diag_np = comp.array_call(s_np)
    tend_jx, diag_jx = comp.array_call(s_jx)
    np.testing.assert_allclose(np.asarray(tend_jx["T"]), tend_np["T"], rtol=1e-6, atol=1e-9)
    for key in ("upwelling_longwave_flux_in_air", "downwelling_longwave_flux_in_air"):
        np.testing.assert_allclose(
            np.asarray(diag_jx[key]), diag_np[key], rtol=1e-6, atol=1e-6)

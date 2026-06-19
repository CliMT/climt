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


def _core_args(nlev, ncol):
    from climt._components.cork.lw.kernels_jax import build_cork_lw_table
    comp, s = _make_state(nlev, ncol, as_jax=True)
    jt = comp._jax_table
    args = dict(
        T=s["T"], p=s["p"], p_int=s["p_int"], T_surf=s["T_surf"],
        h2o=s["h2o"], co2=s["co2"], emissivity=s["emissivity"],
        tau_cloud=s["tau_cloud_lw"], sigma=5.670374419e-8, g=9.80665, cpd=1004.64)
    return jt, args


def test_cork_lw_jax_jit_matches_eager():
    import jax
    from climt._components.cork.lw.kernels_jax import cork_lw_jax
    jt, args = _core_args(20, 4)
    eager = cork_lw_jax(jtable=jt, **args)[0]
    jitted = jax.jit(lambda **kw: cork_lw_jax(jtable=jt, **kw))(**args)[0]
    np.testing.assert_allclose(np.asarray(jitted), np.asarray(eager), rtol=1e-9, atol=1e-12)


def test_cork_lw_jax_vmap_over_columns_matches_loop():
    import jax
    from climt._components.cork.lw.kernels_jax import cork_lw_jax
    jt, args = _core_args(20, 3)

    def one_col(T, p, p_int, T_surf, h2o, co2, emis, tau_cloud):
        # promote a single column (…,1) and take tendency
        t, *_ = cork_lw_jax(T[:, None], p[:, None], p_int[:, None], T_surf[None],
                            h2o[:, None], co2[:, None], emis[:, None],
                            tau_cloud[:, None, :], jt, args["sigma"], args["g"], args["cpd"])
        return t[:, 0]

    vmapped = jax.vmap(one_col, in_axes=(1, 1, 1, 0, 1, 1, 1, 1))(
        args["T"], args["p"], args["p_int"], args["T_surf"],
        args["h2o"], args["co2"], args["emissivity"], args["tau_cloud"])
    full = cork_lw_jax(jtable=jt, **args)[0]
    np.testing.assert_allclose(np.asarray(vmapped.T), np.asarray(full), rtol=1e-7, atol=1e-9)


def test_cork_lw_jax_grad_wrt_state_is_finite():
    import jax
    from climt._components.cork.lw.kernels_jax import cork_lw_jax
    jt, args = _core_args(20, 2)

    def olr(T, T_surf):
        kw = dict(args); kw.pop("T"); kw.pop("T_surf")
        _, up_broad, *_ = cork_lw_jax(T=T, T_surf=T_surf, jtable=jt, **kw)
        return up_broad[-1].sum()

    gT, gTs = jax.grad(olr, argnums=(0, 1))(args["T"], args["T_surf"])
    assert np.all(np.isfinite(np.asarray(gT)))
    assert np.all(np.isfinite(np.asarray(gTs)))
    # OLR increases with surface temperature
    assert float(np.asarray(gTs).sum()) > 0.0

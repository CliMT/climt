# tests/test_cork_optics_jax_parity.py
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest


def _build(tmp_path, k, cont=None):
    from scripts.cork_table_builder.netcdf_writer import write_lw_table
    from climt._components.cork.optics.correlated_k import load_k_table
    ngas, nband, ngpt, nT, nP, nX, nC = k.shape
    write_lw_table(
        str(tmp_path / "t.nc"),
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(200.0, 320.0, nT),
        log_p_grid=np.linspace(np.log(1e2), np.log(1e5), nP),
        band_edges=np.linspace(10.0, 3250.0, nband + 1),
        planck_fraction=np.full((nband, ngpt, nT), 1.0 / ngpt),
        h2o_vmr_grid=np.array([1e-6, 1e-3, 1e-1])[:nX],
        co2_vmr_grid=np.logspace(-5, -2, nC),
        continuum_kappa=cont,
        gas_names=("effective",),
    )
    return load_k_table(str(tmp_path / "t.nc"))


def _inputs(nlev=6, ncol=4, seed=0):
    rng = np.random.RandomState(seed)
    T = rng.uniform(210.0, 310.0, size=(nlev, ncol))
    p = np.exp(rng.uniform(np.log(2e2), np.log(8e4), size=(nlev, ncol)))
    h2o = np.exp(rng.uniform(np.log(1e-6), np.log(1e-1), size=(nlev, ncol)))
    co2 = np.exp(rng.uniform(np.log(1e-5), np.log(1e-2), size=(nlev, ncol)))
    return T, p, h2o, co2


def _prep(table, T, p, h2o, co2):
    """Replicate _additive_co2_fast host prep so JAX gets identical inputs."""
    x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
    c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
    log_p = np.log(np.maximum(p, 1.0))
    x_clamped = np.clip(h2o, float(x_grid[0]), float(x_grid[-1]))
    log_x = np.log(np.maximum(x_clamped, 1e-30))
    c_clamped = np.clip(co2, float(c_grid[0]), float(c_grid[-1]))
    log_c = np.log(np.maximum(c_clamped, 1e-30))
    return log_p, log_x, log_c


def _jax_args(table, T, log_p, log_x, log_c, gas):
    has_cont = ("continuum_kappa" in table
                and np.asarray(table["continuum_kappa"]).ndim == 4)
    log_cont = (np.log(np.maximum(np.asarray(table["continuum_kappa"], np.float64), 1e-40))
                if has_cont else np.zeros((table["k_coefficients"].shape[1], 1, 1, 1)))
    return dict(
        T=jnp.asarray(T), log_p=jnp.asarray(log_p), log_x=jnp.asarray(log_x),
        log_c=jnp.asarray(log_c), gas_amounts=jnp.asarray(gas),
        k=jnp.asarray(table["k_coefficients"], dtype=jnp.float64),
        T_grid=jnp.asarray(table["temperature_grid"], dtype=jnp.float64),
        p_grid_log=jnp.asarray(table["pressure_grid_log"], dtype=jnp.float64),
        log_x_grid=jnp.log(jnp.maximum(jnp.asarray(table["h2o_vmr_grid"], jnp.float64), 1e-30)),
        log_c_grid=jnp.log(jnp.maximum(jnp.asarray(table["co2_vmr_grid"], jnp.float64), 1e-30)),
        has_cont=has_cont, log_cont=jnp.asarray(log_cont), co2_logk=True,
    )


@pytest.mark.parametrize("with_cont", [False, True])
def test_tau_jax_matches_oracle(tmp_path, with_cont):
    from climt._components.cork.optics.correlated_k import compute_ck_optical_depth
    from climt._components.cork.optics.correlated_k_jax import compute_tau_jax
    rng = np.random.RandomState(1)
    k = rng.uniform(1e-6, 1e-2, size=(1, 5, 3, 4, 5, 3, 6))
    cont = rng.uniform(1e-7, 1e-3, size=(5, 4, 5, 3)) if with_cont else None
    table = _build(tmp_path, k, cont=cont)
    T, p, h2o, co2 = _inputs()
    gas = np.full((1,) + T.shape, 90.0)
    tau_oracle = compute_ck_optical_depth(table, T, p, gas, h2o_vmr=h2o, co2_vmr=co2)
    log_p, log_x, log_c = _prep(table, T, p, h2o, co2)
    tau_jax = compute_tau_jax(**_jax_args(table, T, log_p, log_x, log_c, gas))
    np.testing.assert_allclose(np.asarray(tau_jax), tau_oracle, rtol=1e-6, atol=1e-12)


def test_tau_jax_grad_matches_finite_difference(tmp_path):
    from climt._components.cork.optics.correlated_k_jax import compute_tau_jax
    rng = np.random.RandomState(7)
    k = rng.uniform(1e-6, 1e-2, size=(1, 3, 2, 4, 5, 3, 6))
    table = _build(tmp_path, k)
    T, p, h2o, co2 = _inputs(nlev=5, ncol=2, seed=4)
    gas = np.full((1,) + T.shape, 100.0)
    log_p, log_x, log_c = _prep(table, T, p, h2o, co2)
    args = _jax_args(table, T, log_p, log_x, log_c, gas)
    static = {key: args[key] for key in args if key != "T"}

    def scalar(Tj):
        return jnp.sum(compute_tau_jax(T=Tj, **static))

    g = np.asarray(jax.grad(scalar)(args["T"]))
    assert np.all(np.isfinite(g))
    # central finite difference at one interior point well inside grid brackets
    eps = 1e-3
    i, j = 2, 1
    Tp = args["T"].at[i, j].add(eps)
    Tm = args["T"].at[i, j].add(-eps)
    fd = (float(scalar(Tp)) - float(scalar(Tm))) / (2 * eps)
    np.testing.assert_allclose(g[i, j], fd, rtol=1e-3, atol=1e-6)

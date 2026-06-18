import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


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


def _reference_tau(table, T, p, gas, h2o, co2):
    """Oracle: per-level interpolate_k (+ interpolate_continuum) assembly."""
    from climt._components.cork.optics.correlated_k import (
        interpolate_k, interpolate_continuum)
    ngas, nband, ngpt = table["k_coefficients"].shape[:3]
    nlev, ncol = T.shape
    tau = np.zeros((nband, ngpt, nlev, ncol))
    has_cont = ("continuum_kappa" in table
                and np.asarray(table["continuum_kappa"]).ndim == 4)
    for kk in range(nlev):
        ki = interpolate_k(table, T[kk], p[kk], h2o_vmr=h2o[kk], co2_vmr=co2[kk])
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    for ic in range(ncol):
                        tau[ib, igp, kk, ic] += ki[ig, ib, igp, ic] * gas[ig, kk, ic]
        if has_cont:
            cont = interpolate_continuum(table, T[kk], p[kk], h2o[kk])
            for ib in range(nband):
                for igp in range(ngpt):
                    for ic in range(ncol):
                        tau[ib, igp, kk, ic] += cont[ib, ic] * gas[0, kk, ic]
    return tau


def _inputs(nlev=6, ncol=4, seed=0):
    rng = np.random.RandomState(seed)
    T = rng.uniform(210.0, 310.0, size=(nlev, ncol))
    p = np.exp(rng.uniform(np.log(2e2), np.log(8e4), size=(nlev, ncol)))
    h2o = np.exp(rng.uniform(np.log(1e-6), np.log(1e-1), size=(nlev, ncol)))
    co2 = np.exp(rng.uniform(np.log(1e-5), np.log(1e-2), size=(nlev, ncol)))
    return T, p, h2o, co2


def test_njit_co2_path_matches_oracle_no_continuum(tmp_path):
    from climt._components.cork.optics.correlated_k import (
        compute_ck_optical_depth)
    rng = np.random.RandomState(1)
    k = rng.uniform(1e-6, 1e-2, size=(1, 5, 3, 4, 5, 3, 6))
    table = _build(tmp_path, k)
    T, p, h2o, co2 = _inputs()
    gas = np.full((1,) + T.shape, 80.0)
    tau_fast = compute_ck_optical_depth(table, T, p, gas, h2o_vmr=h2o, co2_vmr=co2)
    tau_ref = _reference_tau(table, T, p, gas, h2o, co2)
    np.testing.assert_allclose(tau_fast, tau_ref, rtol=1e-6, atol=1e-12)


def test_njit_co2_path_matches_oracle_with_continuum(tmp_path):
    from climt._components.cork.optics.correlated_k import (
        compute_ck_optical_depth)
    rng = np.random.RandomState(2)
    nband, nT, nP, nX = 5, 4, 5, 3
    k = rng.uniform(1e-6, 1e-2, size=(1, nband, 3, nT, nP, nX, 6))
    cont = rng.uniform(1e-7, 1e-3, size=(nband, nT, nP, nX))
    table = _build(tmp_path, k, cont=cont)
    T, p, h2o, co2 = _inputs(seed=3)
    gas = np.full((1,) + T.shape, 120.0)
    tau_fast = compute_ck_optical_depth(table, T, p, gas, h2o_vmr=h2o, co2_vmr=co2)
    tau_ref = _reference_tau(table, T, p, gas, h2o, co2)
    np.testing.assert_allclose(tau_fast, tau_ref, rtol=1e-6, atol=1e-12)

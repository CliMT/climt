import numpy as np
import pytest
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def _build_co2_table(tmp_path, k):
    from scripts.pf_table_builder.netcdf_writer import write_lw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP, nX, nC = k.shape
    # band_edges must have nband+1 entries
    band_edges = np.linspace(10.0, 3250.0, nband + 1)
    write_lw_table(
        str(tmp_path / "t.nc"),
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(200.0, 300.0, nT),
        log_p_grid=np.linspace(np.log(1e2), np.log(1e5), nP),
        band_edges=band_edges,
        planck_fraction=np.full((nband, ngpt, nT), 1.0 / ngpt),
        h2o_vmr_grid=np.array([1e-6, 1e-3, 1e-1]),
        co2_vmr_grid=np.array([1e-4, 1e-3, 1e-2]),
        gas_names=("effective",),
    )
    return load_k_table(str(tmp_path / "t.nc"))


def test_interpolate_k_recovers_node_value(tmp_path):
    """At an exact grid node, quadrilinear interp returns the stored k."""
    from climt._components.picket_fence.optics.correlated_k import interpolate_k

    rng = np.random.RandomState(0)
    k = rng.uniform(1e-6, 1e-2, size=(1, 4, 2, 3, 4, 3, 3))
    table = _build_co2_table(tmp_path, k)

    # Pick node iT=1 (T=250), iP=2 (mid), iX=1 (1e-3), iC=2 (1e-2).
    T = np.array([250.0])
    p = np.array([np.exp(table["pressure_grid_log"][2])])
    h2o = np.array([1e-3])
    co2 = np.array([1e-2])
    out = interpolate_k(table, T, p, h2o_vmr=h2o, co2_vmr=co2)
    # out shape (ngas, nband, ngpt, ncol)
    np.testing.assert_allclose(
        out[:, :, :, 0], k[:, :, :, 1, 2, 1, 2].astype(np.float32), rtol=2e-5
    )


def test_interpolate_k_co2_geometric_midpoint(tmp_path):
    """Halfway (in log-X_CO2) between two CO2 nodes, log-k interp returns the
    geometric mean of the two node values (T,P,X_H2O all on-node)."""
    from climt._components.picket_fence.optics.correlated_k import interpolate_k

    k = np.ones((1, 1, 1, 3, 4, 3, 3))
    # CO2 nodes 1e-4,1e-3,1e-2; set k=2 at iC=1 and k=8 at iC=2 for one node.
    k[0, 0, 0, 1, 2, 1, 1] = 2.0
    k[0, 0, 0, 1, 2, 1, 2] = 8.0
    table = _build_co2_table(tmp_path, k)

    T = np.array([250.0])
    p = np.array([np.exp(table["pressure_grid_log"][2])])
    h2o = np.array([1e-3])
    co2 = np.array([np.sqrt(1e-3 * 1e-2)])  # geometric midpoint: fC = 0.5
    out = interpolate_k(table, T, p, h2o_vmr=h2o, co2_vmr=co2)
    # geometric mean of 2 and 8 = 4
    np.testing.assert_allclose(out[0, 0, 0, 0], 4.0, rtol=1e-4)

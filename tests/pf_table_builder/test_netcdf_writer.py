import numpy as np
import pytest


def test_write_lw_table_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.netcdf_writer import write_lw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP = 1, 4, 2, 3, 5
    k = np.random.RandomState(0).uniform(1e-7, 1e-2,
                                          size=(ngas, nband, ngpt, nT, nP))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    planck_fraction = np.full((nband, ngpt, nT), 1.0 / ngpt)

    path = tmp_path / "test_lw.nc"
    write_lw_table(
        str(path),
        k_coefficients=k,
        gpoint_weights=weights,
        T_grid=T_grid,
        log_p_grid=log_p_grid,
        band_edges=band_edges,
        planck_fraction=planck_fraction,
        gas_names=("effective",),
        source="linepyline:test",
    )
    loaded = load_k_table(str(path))
    np.testing.assert_allclose(loaded["k_coefficients"], k, rtol=1e-5)
    np.testing.assert_allclose(loaded["gpoint_weights"], weights, rtol=1e-5)
    np.testing.assert_allclose(loaded["temperature_grid"], T_grid, rtol=1e-5)
    np.testing.assert_allclose(loaded["pressure_grid_log"], log_p_grid, rtol=1e-5)


def test_write_lw_table_with_h2o_axis_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.netcdf_writer import write_lw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP, nX = 1, 4, 2, 3, 5, 4
    k = np.random.RandomState(1).uniform(1e-7, 1e-2,
                                          size=(ngas, nband, ngpt, nT, nP, nX))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    h2o_vmr = np.array([1e-6, 1e-4, 1e-2, 1e-1])
    planck_fraction = np.full((nband, ngpt, nT), 1.0 / ngpt)

    path = tmp_path / "test_lw_hmoist.nc"
    write_lw_table(
        str(path),
        k_coefficients=k,
        gpoint_weights=weights,
        T_grid=T_grid,
        log_p_grid=log_p_grid,
        band_edges=band_edges,
        planck_fraction=planck_fraction,
        h2o_vmr_grid=h2o_vmr,
        gas_names=("effective",),
        source="linepyline:test_hmoist",
    )
    loaded = load_k_table(str(path))
    assert loaded["k_coefficients"].shape == (ngas, nband, ngpt, nT, nP, nX)
    np.testing.assert_allclose(loaded["h2o_vmr_grid"], h2o_vmr, rtol=1e-5)


def test_write_sw_table_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.netcdf_writer import write_sw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP = 1, 3, 2, 3, 5
    k = np.random.RandomState(2).uniform(0, 1e-3,
                                          size=(ngas, nband, ngpt, nT, nP))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])
    solar = np.array([[100.0, 100.0], [200.0, 200.0], [300.0, 300.0]])
    rayleigh = np.array([1e-6, 5e-6, 2e-5])

    path = tmp_path / "test_sw.nc"
    write_sw_table(
        str(path),
        k_coefficients=k, gpoint_weights=weights,
        T_grid=T_grid, log_p_grid=log_p_grid, band_edges=band_edges,
        solar_source_per_gpoint=solar, rayleigh_coefficient=rayleigh,
        gas_names=("effective",), source="linepyline:test_sw",
    )
    loaded = load_k_table(str(path))
    np.testing.assert_allclose(loaded["solar_source_per_gpoint"], solar, rtol=1e-5)
    np.testing.assert_allclose(loaded["rayleigh_coefficient"], rayleigh, rtol=1e-5)

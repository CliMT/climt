"""Tests for the linepyline kappa sampler.

These tests are skipped unless linepyline is importable in the active env.
"""
import numpy as np
import pytest

linepyline = pytest.importorskip("linepyline")


def test_sample_kappa_grid_shapes(monkeypatch):
    """sample_kappa_grid returns (nT, nP, nNu) for a single absorber."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.kappa_sampling import sample_kappa_grid

    T_grid = np.array([200.0, 250.0, 300.0])
    p_grid = np.array([1e3, 1e4, 1e5])
    nu_grid = np.arange(500.0, 600.0, 1.0)  # 100 points, narrow LW band

    # Pure CO2 atmosphere, no background gas
    kappa = sample_kappa_grid(
        background_gas=None,
        absorbers={"CO2": 1.0},
        T_grid=T_grid,
        p_grid=p_grid,
        nu_grid=nu_grid,
        line_shape="lorentz",
        binning=True,
    )
    assert kappa.shape == (3, 3, 100)
    assert (kappa >= 0).all()
    # CO2 ν₂ bending mode is in this band — kappa should be non-trivial
    assert kappa.max() > 1e-5


def test_sample_kappa_grid_with_h2o_vmr_axis():
    """sample_kappa_grid with X axis returns (nT, nP, nX, nNu)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.kappa_sampling import sample_kappa_grid

    T_grid = np.array([250.0, 300.0])
    p_grid = np.array([1e3, 1e5])
    X_h2o = np.array([1e-6, 1e-3, 1e-1])  # 3 humidity levels
    nu_grid = np.arange(100.0, 200.0, 1.0)

    kappa = sample_kappa_grid(
        background_gas="N2",
        absorbers={"CO2": 4e-4},
        h2o_vmr_grid=X_h2o,
        T_grid=T_grid,
        p_grid=p_grid,
        nu_grid=nu_grid,
        line_shape="lorentz",
        binning=True,
    )
    assert kappa.shape == (2, 2, 3, 100)
    # Wet column must have more opacity than dry one
    assert kappa[:, :, -1, :].mean() > kappa[:, :, 0, :].mean()

import sys
import numpy as np
from climt._components.cork.optics.correlated_k import load_k_table


def test_npz_table_loads_without_scipy(monkeypatch):
    monkeypatch.setitem(sys.modules, "scipy", None)  # force ImportError on scipy use
    for name in ("single_band_unit_lw", "earth_low_res_lw", "earth_low_res_sw"):
        table = load_k_table(name)
        assert "k_coefficients" in table
        assert np.asarray(table["k_coefficients"]).ndim >= 5

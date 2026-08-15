import sys
import numpy as np
from climt._components.cork.optics.correlated_k import load_k_table


def test_npz_table_loads_without_scipy(monkeypatch):
    monkeypatch.setitem(sys.modules, "scipy", None)  # force ImportError on scipy use
    for name in ("single_band_unit_lw", "earth_low_res_lw", "earth_low_res_sw"):
        table = load_k_table(name)
        assert "k_coefficients" in table
        assert np.asarray(table["k_coefficients"]).ndim >= 5


def test_npz_table_is_materialised_not_lazy(tmp_path):
    """load_k_table must return arrays in memory, never a lazy NpzFile.

    An NpzFile re-inflates the requested array out of the zip on *every*
    __getitem__, and the radiation components index the table repeatedly per
    call, so a lazy return value costs tens of ms per radiation step. Guard
    both the shipped-by-name and the explicit-path branches.
    """
    for table in (load_k_table("earth_low_res_lw"),
                  load_k_table("earth_low_res_sw")):
        assert isinstance(table, dict)
        assert isinstance(table["k_coefficients"], np.ndarray)
        # Repeated lookups hand back the same object rather than re-reading.
        assert table["k_coefficients"] is table["k_coefficients"]

    path = str(tmp_path / "t.npz")
    np.savez(path, k_coefficients=np.ones((1, 2, 2, 3, 3)))
    from_path = load_k_table(path)
    assert isinstance(from_path, dict)
    assert from_path["k_coefficients"] is from_path["k_coefficients"]

import numpy as np
import pytest

SHIPPED = [
    "earth_low_res_lw", "earth_low_res_sw",
    "mars_lw",          "mars_sw",
    "venus_lw",         "venus_sw",
    "trappist1e_hab1_lw", "trappist1e_hab1_sw",
    "trappist1e_hab2_lw", "trappist1e_hab2_sw",
    "titan_lw",         "titan_sw",
]


@pytest.mark.parametrize("name", SHIPPED)
def test_shipped_table_loads_and_has_sane_values(name):
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    table = load_k_table(name)
    k = table["k_coefficients"]
    assert k.shape[0] >= 1
    assert k.shape[1] >= 1
    assert k.shape[2] >= 1
    assert (k >= 0).all()

    w = table["gpoint_weights"]
    np.testing.assert_allclose(w.sum(axis=-1), 1.0, rtol=1e-5)

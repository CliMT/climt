import pytest


def test_condensible_params_h2o_defaults():
    """get_condensible_params() returns H2O values without loading any profile."""
    from climt._core.condensibles import get_condensible_params, SPECIES_ID

    cond = get_condensible_params()
    assert cond.species_id == SPECIES_ID["h2o"]  # 0
    assert abs(cond.RV - 461.5) < 0.1
    assert abs(cond.CPV - 1870.0) < 1.0
    assert abs(cond.CL - 2500.0) < 1.0
    assert abs(cond.LV0 - 2.501e6) < 1e3
    assert abs(cond.ROWL - 1000.0) < 0.1
    assert abs(cond.T_freeze - 273.15) < 0.01


def test_species_id_mapping():
    from climt._core.condensibles import SPECIES_ID

    assert SPECIES_ID["h2o"] == 0
    assert SPECIES_ID["ch4"] == 1
    assert SPECIES_ID["co2"] == 2

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


def test_sat_vap_pressure_h2o_above_freezing():
    """H2O saturation pressure at 273.15 K (0°C) is ~6.112 hPa."""
    from climt._core.condensibles import _sat_vap_pressure

    es = _sat_vap_pressure(273.15, 0)
    assert abs(es - 6.112) < 0.01


def test_sat_vap_pressure_h2o_below_freezing():
    """H2O ice saturation formula is active below 273.15 K."""
    from climt._core.condensibles import _sat_vap_pressure

    # At -20°C (253.15 K), ice saturation should be ~1.03 hPa
    es = _sat_vap_pressure(253.15, 0)
    assert 0.9 < es < 1.2


def test_sat_vap_pressure_ch4_at_boiling():
    """CH4 saturation pressure at 111.65 K equals 1 atm (1013.25 hPa) by construction."""
    from climt._core.condensibles import _sat_vap_pressure

    es = _sat_vap_pressure(111.65, 1)
    assert abs(es - 1013.25) < 1.0


def test_sat_vap_pressure_co2_at_sublimation():
    """CO2 saturation pressure at 194.7 K equals 1 atm by construction."""
    from climt._core.condensibles import _sat_vap_pressure

    es = _sat_vap_pressure(194.7, 2)
    assert abs(es - 1013.25) < 1.0


def test_lv_at_freeze():
    """_lv at T_freeze returns LV0."""
    from climt._core.condensibles import _lv, get_condensible_params

    cond = get_condensible_params()
    lv = _lv(cond.T_freeze, cond)
    assert abs(lv - cond.LV0) < 1.0


def test_lcl_pressure_saturated():
    """_lcl_pressure with RH=1.0 returns P unchanged (both H2O and CC formula)."""
    from climt._core.condensibles import _lcl_pressure, get_condensible_params

    cond = get_condensible_params()
    P, T = 1000.0, 300.0
    plcl = _lcl_pressure(P, 1.0, T, cond)
    assert abs(plcl - P) < 0.1


def test_lcl_pressure_subsaturated():
    """_lcl_pressure with RH < 1.0 returns PLCL < P."""
    from climt._core.condensibles import _lcl_pressure, get_condensible_params

    cond = get_condensible_params()
    plcl = _lcl_pressure(1000.0, 0.7, 300.0, cond)
    assert plcl < 1000.0
    assert plcl > 500.0  # sanity bound


def test_emanuel_titan_ch4_smoke():
    """EmanuelConvectionPython runs without error on a Titan-like CH4 atmosphere."""
    import numpy as np
    from datetime import timedelta
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPython

    load_atmospheric_properties("titan")
    try:
        nlev, ncol = 30, 1
        # Titan surface ~94 K, tropopause ~70 K, lapse rate ~1 K/km
        T = np.linspace(94.0, 70.0, nlev).reshape(nlev, ncol)
        # CH4 specific humidity: ~5% near surface, decreasing with altitude
        q = np.linspace(0.05, 0.001, nlev).reshape(nlev, ncol)
        u = np.zeros((nlev, ncol))
        v = np.zeros((nlev, ncol))
        # Titan surface pressure ~1467 hPa
        p = np.linspace(1467.0, 10.0, nlev).reshape(nlev, ncol)
        ph = np.linspace(1500.0, 5.0, nlev + 1).reshape(nlev + 1, ncol)
        cbmf = np.zeros(ncol)

        conv = EmanuelConvectionPython()
        state = {
            "air_temperature": T,
            "specific_humidity": q,
            "eastward_wind": u,
            "northward_wind": v,
            "air_pressure": p,
            "air_pressure_on_interface_levels": ph,
            "cloud_base_mass_flux": cbmf,
        }
        # Should complete without exception; output values not checked here
        tendencies, diagnostics = conv.array_call(state, timedelta(minutes=15))
        assert "air_temperature" in tendencies
        assert tendencies["air_temperature"].shape == (nlev, ncol)
    finally:
        reset_atmospheric_properties()
    """compute_qs must match bolton_q_sat for H2O above freezing to within rtol=1e-5."""
    import numpy as np
    from climt._core.condensibles import compute_qs, get_condensible_params
    from climt._core.util import bolton_q_sat

    cond = get_condensible_params()
    # All temperatures above 273.15 K so both formulas use the same (liquid) branch
    T = np.array([[280.0], [290.0], [300.0], [310.0]])   # (4, 1)
    P_hpa = np.array([[950.0], [900.0], [850.0], [800.0]])

    qs_new = compute_qs(T, P_hpa, cond, 287.04)
    qs_ref = bolton_q_sat(T, P_hpa * 100.0, 287.04, 461.5)
    np.testing.assert_allclose(qs_new, qs_ref, rtol=1e-5)

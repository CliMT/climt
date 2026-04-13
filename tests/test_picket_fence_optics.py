# tests/test_picket_fence_optics.py
import numpy as np
import pytest


def test_parmentier_coefficients_grey_limit():
    """When gamma_P=1 and R=1, kappa_1 should equal kappa_2 (grey limit)."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_thermal_opacities,
    )

    kappa_R = 1e-3  # m^2/kg
    gamma_P = 1.0
    beta = 0.5
    R = 1.0

    kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)
    np.testing.assert_allclose(kappa_1, kappa_2, rtol=1e-10)
    np.testing.assert_allclose(kappa_1, kappa_R, rtol=1e-10)


def test_parmentier_coefficients_nongrey():
    """When R > 1, kappa_1 > kappa_R > kappa_2."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_thermal_opacities,
    )

    kappa_R = 1e-3
    gamma_P = 10.0
    beta = 0.8
    R = 100.0

    kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)
    assert kappa_1 > kappa_R
    assert kappa_2 < kappa_R
    assert kappa_1 / kappa_2 == pytest.approx(R, rel=1e-6)


def test_parmentier_rosseland_mean():
    """kappa_1 and kappa_2 should reconstruct kappa_R via the Rosseland mean formula."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_thermal_opacities,
    )

    kappa_R = 1e-2
    gamma_P = 5.0
    beta = 0.7
    R = 50.0

    kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)
    # Rosseland mean: 1/kappa_R = beta/kappa_1 + (1-beta)/kappa_2
    # => kappa_R = kappa_1 * kappa_2 / (beta*kappa_2 + (1-beta)*kappa_1)
    kappa_R_reconstructed = kappa_1 * kappa_2 / (beta * kappa_2 + (1 - beta) * kappa_1)
    np.testing.assert_allclose(kappa_R_reconstructed, kappa_R, rtol=1e-10)


def test_freedman2014_opacity_physical_range():
    """Freedman 2014 kappa_R should be in a physically plausible range for gas opacity."""
    from climt._components.picket_fence.optics.parmentier import (
        compute_rosseland_mean_opacity,
        load_freedman2014_coefficients,
    )

    coeffs = load_freedman2014_coefficients()

    # At T=1500 K, P=1e5 Pa (1 bar): typical hot Jupiter photosphere
    # Freedman et al. 2014 gas opacity is roughly 0.001-0.1 m^2/kg in this regime
    kappa = compute_rosseland_mean_opacity(1500.0, 1e5, coeffs)
    assert 1e-4 < kappa < 1.0, f"kappa_R={kappa:.2e} m^2/kg outside plausible range"

    # kappa should increase with temperature (more thermal opacity at higher T)
    kappa_hot = compute_rosseland_mean_opacity(2000.0, 1e5, coeffs)
    kappa_cool = compute_rosseland_mean_opacity(1000.0, 1e5, coeffs)
    assert kappa_hot > kappa_cool, (
        f"opacity should increase with T: kappa(2000K)={kappa_hot:.2e}, "
        f"kappa(1000K)={kappa_cool:.2e}"
    )

    # kappa should increase with pressure (collision-induced absorption)
    kappa_highp = compute_rosseland_mean_opacity(1500.0, 1e6, coeffs)
    kappa_lowp = compute_rosseland_mean_opacity(1500.0, 1e4, coeffs)
    assert kappa_highp > kappa_lowp, (
        f"opacity should increase with P: kappa(1e6)={kappa_highp:.2e}, "
        f"kappa(1e4)={kappa_lowp:.2e}"
    )


def test_parmentier_lookup_coefficients():
    """Lookup of ratio coefficients from T_eff should return plausible values."""
    from climt._components.picket_fence.optics.parmentier import (
        load_parmentier_coefficients,
        lookup_ratio_coefficients,
    )

    coeffs = load_parmentier_coefficients("solar_composition")
    # T_eff = 1500 K is in the 600-1400 K range for solar composition
    gamma_v1, gamma_v2, gamma_v3, beta, gamma_P, R = lookup_ratio_coefficients(
        coeffs, 1500.0
    )

    assert gamma_v1 > 0
    assert gamma_v2 > 0
    assert gamma_v3 > 0
    assert 0 < beta <= 1
    assert gamma_P >= 1
    assert R >= 1


def test_correlated_k_load_table():
    """Loading a test table returns expected dimensions."""
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    table = load_k_table("test_2band_lw")
    assert table["k_coefficients"].shape[0] == 1  # 1 gas
    assert table["k_coefficients"].shape[1] == 2  # 2 bands
    assert table["k_coefficients"].shape[2] == 2  # 2 g-points
    assert table["gpoint_weights"].shape == (2, 2)
    # Weights sum to 1 per band
    for b in range(2):
        np.testing.assert_allclose(table["gpoint_weights"][b].sum(), 1.0)


def test_correlated_k_interpolation():
    """k-coefficient interpolation returns plausible values."""
    from climt._components.picket_fence.optics.correlated_k import (
        interpolate_k,
        load_k_table,
    )

    table = load_k_table("test_2band_lw")
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]

    # Interpolate at a grid point — should return exact tabulated value
    T_val = T_grid[1]
    p_log_val = p_grid_log[1]
    k = table["k_coefficients"]  # (ngas, nband, ngpt, nT, np)

    result = interpolate_k(table, np.array([T_val]), np.array([np.exp(p_log_val)]))
    # result shape: (ngas, nband, ngpt, 1)
    expected = k[:, :, :, 1, 1]
    np.testing.assert_allclose(result[:, :, :, 0], expected, rtol=1e-6)


def test_correlated_k_lw_co2_units_are_mole_fraction(tmp_path):
    """CO2 and other non-H2O gases must be declared as mole/mole, not kg/kg."""
    from climt._components.picket_fence import PicketFenceLongwave

    # Build a minimal 2-gas table (h2o + co2) for this test
    table_file = str(tmp_path / "two_gas_lw.npz")
    ngas, nband, ngpt, nT, nP = 2, 2, 2, 5, 5
    np.savez(
        table_file,
        k_coefficients=np.ones((ngas, nband, ngpt, nT, nP)) * 0.01,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        planck_fraction=np.full((nband, ngpt, nT), 0.5),
        band_wavenumber_limits=np.array([[200.0, 800.0], [800.0, 1400.0]]),
        temperature_grid=np.linspace(200.0, 400.0, nT),
        pressure_grid_log=np.linspace(np.log(100.0), np.log(1e5), nP),
        gas_names=np.array(["h2o", "co2"]),
        overlap_method=np.array("additive"),
        resolution=np.array("low"),
    )

    lw = PicketFenceLongwave(optics="correlated_k", table=table_file)
    props = lw.input_properties

    # H2O uses mass mixing ratio
    assert props["specific_humidity"]["units"] == "kg/kg", (
        f"H2O units should be kg/kg, got {props['specific_humidity']['units']}"
    )
    # CO2 is a mole fraction — must NOT be kg/kg
    co2_key = "mole_fraction_of_carbon_dioxide_in_air"
    assert co2_key in props, "CO2 input property missing"
    assert props[co2_key]["units"] == "mole/mole", (
        f"CO2 units should be mole/mole, got {props[co2_key]['units']}"
    )


def test_correlated_k_optical_depth():
    """Correlated-k optics produces optical depths with correct shape."""
    from climt._components.picket_fence.optics.correlated_k import (
        compute_ck_optical_depth,
        load_k_table,
    )

    table = load_k_table("test_2band_lw")
    nlev = 10
    ncol = 2
    ngas = 1
    nband = 2
    ngpt = 2

    T = 250.0 * np.ones((nlev, ncol))
    p = np.linspace(100.0, 1e5, nlev).reshape(nlev, 1) * np.ones((1, ncol))
    gas_amounts = 1e-3 * np.ones((ngas, nlev, ncol))  # kg/m^2

    tau = compute_ck_optical_depth(table, T, p, gas_amounts)
    assert tau.shape == (nband, ngpt, nlev, ncol)
    assert np.all(tau >= 0)

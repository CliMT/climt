import numpy as np
import pytest
import sympl


def test_picket_fence_lw_parmentier_runs():
    """PicketFenceLongwave in Parmentier mode produces non-zero fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([lw], grid_state=grid)

    tendencies, diagnostics = lw(state)

    assert np.any(tendencies["air_temperature"].values != 0.0)
    assert not np.any(np.isnan(tendencies["air_temperature"].values))
    assert np.any(diagnostics["upwelling_longwave_flux_in_air"].values != 0.0)
    assert np.any(diagnostics["downwelling_longwave_flux_in_air"].values != 0.0)


def test_picket_fence_lw_per_band_sum():
    """Per-band fluxes should sum to broadband."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=4, ny=1, nz=20)
    state = get_default_state([lw], grid_state=grid)

    tendencies, diagnostics = lw(state)

    up_band = diagnostics["upwelling_longwave_flux_in_air_per_band"].values
    up_broad = diagnostics["upwelling_longwave_flux_in_air"].values

    # Sum over band dimension (last axis)
    np.testing.assert_allclose(up_band.sum(axis=-1), up_broad, rtol=1e-10)


def test_picket_fence_lw_isothermal_equilibrium():
    """An isothermal atmosphere should have near-zero heating rate."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([lw], grid_state=grid)

    # Set everything to the same temperature
    state["air_temperature"].values[:] = 280.0
    state["surface_temperature"].values[:] = 280.0

    tendencies, diagnostics = lw(state)

    # Heating rate should be small (not exactly zero due to optical depth structure,
    # but should be much smaller than a non-isothermal case)
    max_hr = np.abs(diagnostics["air_temperature_tendency_from_longwave"].values).max()
    assert max_hr < 5.0  # K/day — loose bound for isothermal


def test_surface_emissivity_reduces_olr():
    """Setting surface emissivity < 1 should reduce OLR compared to emissivity = 1."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=20)

    state_full = get_default_state([lw], grid_state=grid)
    # Default emissivity should come from state (ones after our fix)
    _, diag_full = lw(state_full)
    olr_full = diag_full["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    state_low = get_default_state([lw], grid_state=grid)
    state_low["surface_longwave_emissivity"].values[:] = 0.5
    _, diag_low = lw(state_low)
    olr_low = diag_low["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    assert olr_low < olr_full, (
        f"Lower emissivity should reduce OLR: olr_full={olr_full:.2f}, olr_low={olr_low:.2f}"
    )


def test_picket_fence_lw_correlated_k_runs():
    """PicketFenceLongwave in correlated-k mode produces non-zero fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="correlated_k", table="test_2band_lw")
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([lw], grid_state=grid)

    state["specific_humidity"].values[:] = 1e-3

    tendencies, diagnostics = lw(state)

    assert np.any(tendencies["air_temperature"].values != 0.0)
    assert not np.any(np.isnan(tendencies["air_temperature"].values))


def test_picket_fence_lw_correlated_k_per_band_sum():
    """Per-band fluxes sum to broadband in correlated-k mode."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="correlated_k", table="test_2band_lw")
    grid = get_grid(nx=2, ny=1, nz=20)
    state = get_default_state([lw], grid_state=grid)

    state["specific_humidity"].values[:] = 1e-3

    tendencies, diagnostics = lw(state)

    up_band = diagnostics["upwelling_longwave_flux_in_air_per_band"].values
    up_broad = diagnostics["upwelling_longwave_flux_in_air"].values
    np.testing.assert_allclose(up_band.sum(axis=-1), up_broad, rtol=1e-10)


def test_lw_per_band_optical_depth_diagnostics():
    """LW component returns per-band optical depth and transmittance."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([lw], grid_state=grid)

    tend, diag = lw(state)
    assert "longwave_optical_depth_per_band" in diag
    assert "longwave_transmittance_per_band" in diag
    assert "air_temperature_tendency_from_longwave_per_band" in diag

    tau = diag["longwave_optical_depth_per_band"].values
    trans = diag["longwave_transmittance_per_band"].values

    # Optical depth must be non-negative
    assert np.all(tau >= 0)
    # Transmittance must be in [0, 1]
    assert np.all(trans >= 0)
    assert np.all(trans <= 1)
    # Transmittance should be consistent with optical depth (exp(-D*tau))
    D = 1.66
    expected_trans = np.exp(-D * tau)
    np.testing.assert_allclose(trans, expected_trans, rtol=1e-10)


def test_lw_per_band_heating_rate_sums_to_broadband():
    """Per-band heating rates should sum to the broadband heating rate."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([lw], grid_state=grid)

    tend, diag = lw(state)
    hr_broad = diag["air_temperature_tendency_from_longwave"].values
    hr_band = diag["air_temperature_tendency_from_longwave_per_band"].values
    # Sum over band dimension (last axis)
    hr_band_sum = hr_band.sum(axis=-1)
    np.testing.assert_allclose(hr_band_sum, hr_broad, rtol=1e-10)


def test_picket_fence_lw_co2_axis_changes_olr(tmp_path):
    """A premixed-bg table with a CO2 axis reads CO2 from state: more CO2 ->
    less OLR (more trapping)."""
    import sys
    sys.path.insert(0, str(tmp_path.parent.parent))
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from scripts.pf_table_builder.netcdf_writer import write_lw_table
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    # Tiny 2-band CO2-axis table; band 0 k rises steeply with CO2 index.
    nband, ngpt, nT, nP, nX, nC = 2, 2, 4, 4, 3, 3
    k = np.full((1, nband, ngpt, nT, nP, nX, nC), 1e-5)
    # band 0 k increases with CO2 index
    k[0, 0, :, :, :, :, :] = np.array([1e-5, 1e-3, 1e-1])[None, None, None, None, :]
    write_lw_table(
        str(tmp_path / "co2tab.nc"),
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(180.0, 320.0, nT),
        log_p_grid=np.linspace(np.log(1e2), np.log(1e5), nP),
        band_edges=np.array([10.0, 1000.0, 3250.0]),
        planck_fraction=np.full((nband, ngpt, nT), 1.0 / nband),
        h2o_vmr_grid=np.array([1e-5, 1e-3, 1e-1]),
        co2_vmr_grid=np.array([1e-4, 1e-3, 1e-2]),
        gas_names=("effective",),
    )

    lw = PicketFenceLongwave(optics="correlated_k", table=str(tmp_path / "co2tab.nc"))
    assert lw._has_co2_axis is True
    grid = get_grid(nx=1, ny=1, nz=20)

    state_lo = get_default_state([lw], grid_state=grid)
    state_lo["specific_humidity"].values[:] = 1e-3
    state_lo["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1e-4
    _, diag_lo = lw(state_lo)
    olr_lo = diag_lo["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    state_hi = get_default_state([lw], grid_state=grid)
    state_hi["specific_humidity"].values[:] = 1e-3
    state_hi["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1e-2
    _, diag_hi = lw(state_hi)
    olr_hi = diag_hi["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    assert olr_hi < olr_lo, f"more CO2 should trap more: lo={olr_lo}, hi={olr_hi}"


def test_shipped_earth_hifi_lw_is_co2_adjustable():
    """The default earth_low_res_lw table loads by name, is 14-band,
    exposes a CO2 axis, and OLR responds to CO2 (quadrupling traps more LW)."""
    table_name = "earth_low_res_lw"
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())
    lw = PicketFenceLongwave(optics="correlated_k", table=table_name)
    assert lw._has_co2_axis is True
    assert lw.num_longwave_bands == 14

    grid = get_grid(nx=1, ny=1, nz=30)

    s_lo = get_default_state([lw], grid_state=grid)
    s_lo["specific_humidity"].values[:] = 5e-3
    s_lo["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 280e-6
    _, d_lo = lw(s_lo)
    olr_lo = d_lo["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    s_hi = get_default_state([lw], grid_state=grid)
    s_hi["specific_humidity"].values[:] = 5e-3
    s_hi["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1120e-6
    _, d_hi = lw(s_hi)
    olr_hi = d_hi["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    assert olr_hi < olr_lo, f"CO2 quadrupling should trap more: lo={olr_lo}, hi={olr_hi}"
    assert not np.isnan(olr_lo) and not np.isnan(olr_hi)

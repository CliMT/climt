# tests/test_picket_fence_integration.py
import numpy as np
import pytest
import sympl


def test_lw_sw_combined_energy_balance():
    """Combined LW+SW fluxes produce physically sensible results."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave, PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="parmentier")
    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([lw, sw], grid_state=grid)

    state["zenith_angle"].values[:] = np.pi / 3  # daytime
    state["surface_albedo_for_direct_shortwave"].values[:] = 0.3
    state["irradiation_temperature"].values[:] = (
        300.0  # SW flux to ensure net heating at TOA
    )

    tend_lw, diag_lw = lw(state)
    tend_sw, diag_sw = sw(state)

    # OLR should be positive (upward LW at TOA)
    olr = diag_lw["upwelling_longwave_flux_in_air"].values[-1, 0, 0]
    assert olr > 0

    # Downward SW at surface should be positive
    sw_sfc = diag_sw["downwelling_shortwave_flux_in_air"].values[0, 0, 0]
    assert sw_sfc > 0

    # Net tendency should have cooling at some levels and heating at others
    total_tend = tend_lw["air_temperature"].values + tend_sw["air_temperature"].values
    assert np.any(total_tend > 0)  # some heating
    assert np.any(total_tend < 0)  # some cooling


def test_cloud_optical_depth_modifies_fluxes():
    """Adding cloud optical depth changes the LW fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave

    sympl.set_backend(sympl.DataArrayBackend())

    lw = PicketFenceLongwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=20)
    state_clear = get_default_state([lw], grid_state=grid)

    _, diag_clear = lw(state_clear)
    olr_clear = diag_clear["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    # Now add cloud in the middle of the atmosphere
    state_cloudy = get_default_state([lw], grid_state=grid)
    assert "longwave_optical_thickness_due_to_cloud" in state_cloudy, (
        "Cloud tau not in default state — input_properties missing it"
    )
    cloud_tau = state_cloudy["longwave_optical_thickness_due_to_cloud"].values
    cloud_tau[10, ...] = 5.0  # thick cloud at level 10, all columns, all bands

    _, diag_cloudy = lw(state_cloudy)
    olr_cloudy = diag_cloudy["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    # Cloud should reduce OLR (more LW trapped)
    assert olr_cloudy < olr_clear, (
        f"Cloud should reduce OLR: clear={olr_clear:.2f}, cloudy={olr_cloudy:.2f}"
    )

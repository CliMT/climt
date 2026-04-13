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
    max_hr = np.abs(diagnostics["longwave_heating_rate"].values).max()
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

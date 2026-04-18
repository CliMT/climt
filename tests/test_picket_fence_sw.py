# tests/test_picket_fence_sw.py
import numpy as np
import pytest
import sympl


def test_picket_fence_sw_parmentier_runs():
    """PicketFenceShortwave in Parmentier mode produces non-zero fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([sw], grid_state=grid)

    # Set non-zero zenith angle (daytime)
    state["zenith_angle"].values[:] = np.pi / 4

    tendencies, diagnostics = sw(state)

    assert np.any(diagnostics["downwelling_shortwave_flux_in_air"].values != 0.0)
    assert not np.any(np.isnan(tendencies["air_temperature"].values))


def test_picket_fence_sw_nighttime_zero():
    """With zenith angle >= pi/2, all SW fluxes should be zero."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)

    state["zenith_angle"].values[:] = np.pi / 2 + 0.1  # nighttime

    tendencies, diagnostics = sw(state)

    np.testing.assert_allclose(
        diagnostics["downwelling_shortwave_flux_in_air"].values, 0.0, atol=1e-10
    )
    np.testing.assert_allclose(tendencies["air_temperature"].values, 0.0, atol=1e-10)


def test_flux_adjustment_scales_toa_solar():
    """flux_adjustment_for_earth_sun_distance should scale TOA solar flux proportionally."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=20)

    state_ref = get_default_state([sw], grid_state=grid)
    state_ref["zenith_angle"].values[:] = np.pi / 4
    state_ref["flux_adjustment_for_earth_sun_distance"].values[()] = 1.0

    state_near = get_default_state([sw], grid_state=grid)
    state_near["zenith_angle"].values[:] = np.pi / 4
    state_near["flux_adjustment_for_earth_sun_distance"].values[()] = (
        1.1  # 10% more flux
    )

    _, diag_ref = sw(state_ref)
    _, diag_near = sw(state_near)

    toa_ref = diag_near["downwelling_shortwave_flux_in_air"].values[-1, 0, 0]
    toa_near = diag_ref["downwelling_shortwave_flux_in_air"].values[-1, 0, 0]
    np.testing.assert_allclose(toa_ref / toa_near, 1.1, rtol=1e-6)


def test_picket_fence_sw_per_band_sum():
    """Per-band fluxes sum to broadband in SW."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=4, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 3

    tendencies, diagnostics = sw(state)

    down_band = diagnostics["downwelling_shortwave_flux_in_air_per_band"].values
    down_broad = diagnostics["downwelling_shortwave_flux_in_air"].values
    np.testing.assert_allclose(down_band.sum(axis=-1), down_broad, rtol=1e-10)


def test_picket_fence_sw_correlated_k_runs():
    """PicketFenceShortwave in correlated-k mode produces non-zero fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="correlated_k", table="test_2band_sw")
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)

    state["zenith_angle"].values[:] = np.pi / 4

    tendencies, diagnostics = sw(state)

    assert np.any(diagnostics["downwelling_shortwave_flux_in_air"].values != 0.0)
    assert not np.any(np.isnan(tendencies["air_temperature"].values))


def test_picket_fence_sw_correlated_k_band_sum():
    """Per-band fluxes sum to broadband in correlated-k SW mode."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="correlated_k", table="test_2band_sw")
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 3

    _, diagnostics = sw(state)

    down_band = diagnostics["downwelling_shortwave_flux_in_air_per_band"].values
    down_broad = diagnostics["downwelling_shortwave_flux_in_air"].values
    np.testing.assert_allclose(down_band.sum(axis=-1), down_broad, rtol=1e-10)


def test_picket_fence_sw_correlated_k_nighttime():
    """Correlated-k SW produces zero flux at nighttime."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave

    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="correlated_k", table="test_2band_sw")
    grid = get_grid(nx=1, ny=1, nz=15)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 2 + 0.1

    tendencies, diagnostics = sw(state)

    np.testing.assert_allclose(
        diagnostics["downwelling_shortwave_flux_in_air"].values, 0.0, atol=1e-10
    )


def test_picket_fence_sw_stellar_spectrum_loads():
    """PicketFenceShortwave with stellar_spectrum='sun' uses non-default solar flux."""
    import numpy as np

    from climt._components.picket_fence import PicketFenceShortwave

    sw = PicketFenceShortwave(optics="parmentier", stellar_spectrum="sun")
    # The component should have loaded the spectrum and have non-trivial flux per band
    # (not exactly 1361/3 per band, since the spectrum is non-uniform)
    flux = sw._solar_flux_per_band
    assert flux.shape == (3,)
    assert np.all(flux > 0)
    # Total should be approximately 1361 W/m^2
    np.testing.assert_allclose(flux.sum(), 1361.0, rtol=0.05)
    # Bands should NOT be exactly equal (spectrum is non-uniform)
    assert not np.allclose(flux, flux[0], rtol=0.01)


def test_sw_per_band_optical_depth_diagnostics():
    """SW component returns per-band optical depth."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave
    import sympl
    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 4

    tend, diag = sw(state)
    assert "shortwave_optical_depth_per_band" in diag
    assert "shortwave_heating_rate_per_band" in diag

    tau = diag["shortwave_optical_depth_per_band"].values
    assert np.all(tau >= 0)


def test_sw_per_band_heating_rate_sums_to_broadband():
    """Per-band SW heating rates should sum to the broadband heating rate."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceShortwave
    import sympl
    sympl.set_backend(sympl.DataArrayBackend())

    sw = PicketFenceShortwave(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 4

    tend, diag = sw(state)
    hr_broad = diag["shortwave_heating_rate"].values
    hr_band = diag["shortwave_heating_rate_per_band"].values
    hr_band_sum = hr_band.sum(axis=-1)
    np.testing.assert_allclose(hr_band_sum, hr_broad, rtol=1e-10)


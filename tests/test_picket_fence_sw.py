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

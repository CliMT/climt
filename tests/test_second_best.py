import numpy as np
from datetime import timedelta
from climt import SecondBEST, get_default_state, get_grid


def test_soil_grid_present():
    grid = get_grid(nx=4, ny=2, nz=10, n_soil_interface_levels=4)
    assert "height_on_soil_interface_levels" in grid
    assert grid["height_on_soil_interface_levels"].shape[0] == 4


class _StubFluxes:
    called = False
    def __call__(self, drag, atmos, soil, soil_props, timestep):
        _StubFluxes.called = True
        return {"sensible_heat_flux": 12.0, "latent_heat_flux": 34.0,
                "momentum_flux": np.array([0.0, 0.0]), "evaporation": 0.0}


def test_second_best_uses_injected_process():
    comp = SecondBEST(fluxes=_StubFluxes())
    state = get_default_state([comp], grid_state=get_grid(
        nx=1, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    diag, new = comp(state, timedelta(seconds=100))
    assert _StubFluxes.called
    assert np.isclose(diag["surface_upward_sensible_heat_flux"].values[0], 12.0)


def test_second_best_passes_through_non_land_columns_unchanged():
    comp = SecondBEST()
    state = get_default_state([comp], grid_state=get_grid(
        nx=2, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    state["area_type"].values[0, 1] = "sea"

    # Give the sea column distinct, non-zero values that must survive a
    # step untouched since SecondBEST only computes land/land_ice physics.
    state["surface_temperature"].values[0, 1] = 290.0
    state["soil_temperature"].values[:, 0, 1] = 295.0

    land_T_before = float(state["surface_temperature"].values[0, 0])

    diag, new = comp(state, timedelta(seconds=100))

    # Sea column: pass-through, unchanged.
    assert np.isclose(new["surface_temperature"].values[0, 1], 290.0)
    assert np.allclose(new["soil_temperature"].values[:, 0, 1], 295.0)

    # Land column: actually stepped (should differ from its initial value
    # given nonzero fluxes/temperature gradients in the default state).
    assert not np.isclose(new["surface_temperature"].values[0, 0], land_T_before)


def test_second_best_reaches_reasonable_state():
    comp = SecondBEST()
    state = get_default_state([comp], grid_state=get_grid(
        nx=2, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    state["downwelling_shortwave_flux_in_air"].values[:] = 300.0
    state["downwelling_longwave_flux_in_air"].values[:] = 300.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 60.0
    state["upwelling_longwave_flux_in_air"].values[:] = 350.0
    for _ in range(50):
        diag, new = comp(state, timedelta(seconds=600))
        state.update(new); state.update(diag)
    T = state["surface_temperature"].values
    assert np.all(np.isfinite(T))
    assert np.all((T > 200.0) & (T < 350.0))


def test_extending_second_best_swap_example():
    """Documents how to swap a process — a constant-drag SurfaceLayer."""
    from climt._components.second_best.processes import SurfaceLayer

    class ConstantDrag(SurfaceLayer):
        def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
            return {"C_Dm": 0.002, "C_Dh": 0.002, "C_DN": 0.002, "Ri": 0.0}

    comp = SecondBEST(surface_layer=ConstantDrag())
    state = get_default_state([comp], grid_state=get_grid(
        nx=1, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    diag, _ = comp(state, timedelta(seconds=100))
    assert np.isclose(diag["surface_drag_coefficient_for_heat"].values[0], 0.002)


def _land_state_for_drag():
    comp = SecondBEST()
    state = get_default_state([comp], grid_state=get_grid(
        nx=1, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    # unstable-neutral-ish, nonzero wind so drag is well defined
    state["eastward_wind"].values[:] = 5.0
    state["surface_temperature"].values[:] = 290.0
    return comp, state


def test_surface_layer_reference_height_uses_lowest_model_level():
    # The drag reference height is the hypsometric height of the lowest model
    # level above the surface, derived from surface_air_pressure vs the lowest
    # mid-level pressure. The old code used a fixed scale height (Rd*T/g,
    # ~8-9 km) independent of the pressure gap, so changing surface_air_pressure
    # would NOT move the drag. Here we make the lowest layer thicker (larger
    # p_surface/p_lowest gap) and require the neutral drag coefficient to drop,
    # which can only happen if the reference height tracks the pressure gap.
    comp, state = _land_state_for_drag()
    p_lowest = float(state["air_pressure"].values[0, 0, 0])

    # thin lowest layer: surface pressure just above lowest mid-level
    state["surface_air_pressure"].values[:] = p_lowest + 200.0
    d_thin, _ = comp(state, timedelta(seconds=1))
    cdn_thin = d_thin["surface_drag_coefficient_for_momentum"].values.item(0)

    # thick lowest layer: much larger gap -> larger reference height
    comp2, state2 = _land_state_for_drag()
    state2["surface_air_pressure"].values[:] = p_lowest + 3000.0
    d_thick, _ = comp2(state2, timedelta(seconds=1))
    cdn_thick = d_thick["surface_drag_coefficient_for_momentum"].values.item(0)

    # taller reference height -> smaller neutral drag coefficient
    assert cdn_thick < cdn_thin
    # and both are physically plausible surface drag coefficients (~1e-3),
    # not the ~9e-4 the km-scale reference height produced
    assert 1e-3 < cdn_thin < 6e-3


def test_screen_level_diagnostics_are_physical():
    comp = SecondBEST()
    state = get_default_state([comp], grid_state=get_grid(
        nx=1, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    state["eastward_wind"].values[:] = 6.0
    state["northward_wind"].values[:] = -2.0
    state["surface_temperature"].values[:] = 298.0
    diag, _ = comp(state, timedelta(seconds=100))

    T_surf = 298.0
    T_air = float(state["air_temperature"].values[0, 0, 0])
    t2m = diag["air_temperature_at_2m"].values.item(0)
    # 2m temperature lies between the surface and the lowest-level air value
    assert min(T_surf, T_air) - 1e-9 <= t2m <= max(T_surf, T_air) + 1e-9

    q2m = diag["specific_humidity_at_2m"].values.item(0)
    assert np.isfinite(q2m) and q2m >= 0.0

    u10 = diag["eastward_wind_at_10m"].values.item(0)
    v10 = diag["northward_wind_at_10m"].values.item(0)
    spd_lo = np.sqrt(6.0 ** 2 + 2.0 ** 2)
    spd10 = np.sqrt(u10 ** 2 + v10 ** 2)
    # 10 m is below the lowest model level, so the 10 m speed is reduced
    assert 0.0 < spd10 < spd_lo
    # direction is preserved (eastward positive, northward negative)
    assert u10 > 0.0 and v10 < 0.0
    assert np.isclose(u10 / v10, 6.0 / -2.0)


def test_screen_level_diagnostics_zero_on_non_land():
    comp = SecondBEST()
    state = get_default_state([comp], grid_state=get_grid(
        nx=2, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    state["area_type"].values[0, 1] = "sea"
    diag, _ = comp(state, timedelta(seconds=100))
    # SecondBEST does not own "sea": its screen diagnostics stay at the default
    assert diag["air_temperature_at_2m"].values[0, 1] == 0.0
    assert diag["eastward_wind_at_10m"].values[0, 1] == 0.0
    # the land column was populated
    assert diag["air_temperature_at_2m"].values[0, 0] > 0.0

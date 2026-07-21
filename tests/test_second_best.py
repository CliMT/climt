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

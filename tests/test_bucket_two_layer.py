from datetime import timedelta
import numpy as np
from climt import BucketHydrology, get_default_state, get_grid


def _one_layer_state(smax):
    comp = BucketHydrology(soil_moisture_max=smax)
    state = get_default_state([comp], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["downwelling_shortwave_flux_in_air"].values[:] = 40.0
    state["downwelling_longwave_flux_in_air"].values[:] = 40.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 40.0
    state["upwelling_longwave_flux_in_air"].values[:] = 40.0
    state["soil_layer_thickness"].values[:] = 1.0
    state["surface_temperature"].values[:] = 300.0
    state["stratiform_precipitation_rate"].values[:] = 0.5
    state["convective_precipitation_rate"].values[:] = 0.5
    state["lwe_thickness_of_soil_moisture_content"].values[:] = 0.4
    return comp, state


def test_soil_moisture_clamped_to_configured_max():
    comp, state = _one_layer_state(smax=0.4)
    _, new = comp(state, timedelta(seconds=1))
    # heavy precip; moisture must be capped at the configured max (0.4), not 0.15
    assert np.all(new["lwe_thickness_of_soil_moisture_content"].values <= 0.4 + 1e-12)
    assert np.all(new["lwe_thickness_of_soil_moisture_content"].values > 0.15)


def _two_layer_state():
    comp = BucketHydrology(num_layers=2, soil_moisture_max=0.15,
                           deep_soil_moisture_max=0.5,
                           moisture_diffusion_timescale=86400.0)
    state = get_default_state([comp], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["downwelling_shortwave_flux_in_air"].values[:] = 200.0
    state["downwelling_longwave_flux_in_air"].values[:] = 300.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 50.0
    state["upwelling_longwave_flux_in_air"].values[:] = 300.0
    state["soil_layer_thickness"].values[:] = 1.0
    state["surface_temperature"].values[:] = 295.0
    state["deep_soil_temperature"].values[:] = 290.0
    state["stratiform_precipitation_rate"].values[:] = 0.001
    state["convective_precipitation_rate"].values[:] = 0.0
    state["lwe_thickness_of_soil_moisture_content"].values[:] = 0.08
    state["deep_soil_moisture_content"].values[:] = 0.20
    return comp, state


def test_two_layer_total_water_conserved():
    comp, state = _two_layer_state()
    dt = timedelta(seconds=100)
    diag, new = comp(state, dt)
    w0 = (state["lwe_thickness_of_soil_moisture_content"].values
          + state["deep_soil_moisture_content"].values)
    w1 = (new["lwe_thickness_of_soil_moisture_content"].values
          + new["deep_soil_moisture_content"].values)
    P = (state["stratiform_precipitation_rate"].values
         + state["convective_precipitation_rate"].values)
    E = diag["evaporation_rate"].values
    R = diag["runoff_rate"].values
    # no drainage set -> total water change = (P - E - R) * dt
    assert np.isclose((w1 - w0), (P - E - R) * dt.total_seconds(), rtol=0, atol=1e-9)


def test_two_layer_deep_store_is_slower_than_shallow():
    comp, state = _two_layer_state()
    dt = timedelta(seconds=3600)
    # Drive with an oscillating precip signal and compare the step-to-step
    # variability of the shallow vs. deep store.
    #
    # Deviation from the brief's literal numbers, and why: the brief's
    # original amplitude (0.001 m/s applied to convective_precipitation_rate,
    # on top of the fixture's 0.001 m/s stratiform baseline) saturates the
    # 0.15 m shallow bucket on the very first step (0.001 m/s * 3600 s = 3.6 m
    # >> 0.15 m capacity) and pins it there for the whole run, so
    # std(diff(shallow)) == 0 for any correct implementation -- this is a
    # units/magnitude issue in the test fixture, not the two-layer physics.
    # Separately, convective_precipitation_rate's default state array carries
    # units "mm day^-1" (climt/_core/initialization.py) while the component
    # consumes "m s^-1"; writing raw m/s-scale floats into `.values` is
    # silently divided by ~8.64e7 by sympl's unit conversion on the next
    # `comp()` call, so the intended oscillation barely registers regardless.
    # Here the oscillation instead drives `stratiform_precipitation_rate`
    # (units already "m s^-1", no silent conversion) at a magnitude that
    # perturbs the shallow store without permanently saturating it, which
    # preserves the property under test: the shallow store responds more
    # strongly to the fast oscillation than the deep store does.
    shallow, deep = [], []
    for step in range(200):
        # subseasonal (fast) precip signal
        fast = 2e-6 * (1 + np.sin(2 * np.pi * step / 24))
        state["stratiform_precipitation_rate"].values[:] = max(fast, 0.0)
        state["convective_precipitation_rate"].values[:] = 0.0
        diag, new = comp(state, dt)
        state.update(new); state.update(diag)
        shallow.append(float(new["lwe_thickness_of_soil_moisture_content"].values[0]))
        deep.append(float(new["deep_soil_moisture_content"].values[0]))
    # shallow store tracks the fast forcing more strongly than the deep store
    assert np.std(np.diff(shallow)) > np.std(np.diff(deep))

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
    # Drive with a sustained oscillating precip signal and compare the
    # step-to-step variability of the shallow vs. deep store, once both are
    # tracking the periodic forcing in quasi-equilibrium rather than
    # responding to initial transients.
    #
    # Two artifacts to avoid, found by instrumenting an earlier version of
    # this test:
    #
    # 1. Saturation clipping. This fixture has zero evaporation (wind speed
    #    defaults to 0, so potential_evaporation == 0), and the shared
    #    `_two_layer_state()` fixture sets no `deep_drainage_timescale`, so
    #    the only way water leaves the shallow store is the shallow<->deep
    #    exchange flux, and nothing leaves the deep store at all. Under any
    #    positive-mean forcing both stores climb monotonically toward their
    #    caps; once the shallow store pins at 0.15 m, std(diff(shallow))
    #    stops measuring "tracks the oscillation" and starts measuring
    #    "ran into a wall". Fix: build a local two-layer state (not the
    #    shared fixture) with `deep_drainage_timescale` set, so the deep
    #    store can lose water and the whole system settles into a periodic
    #    quasi-equilibrium instead of saturating.
    # 2. Transient drainage. The shared fixture starts the shallow and deep
    #    stores at different relative saturations (0.08/0.15 = 53% vs.
    #    0.20/0.5 = 40%), so the first many steps are dominated by the
    #    shallow<->deep exchange flux draining that initial imbalance, not
    #    by the oscillating forcing. Fix: start both layers at the same
    #    relative saturation, and discard an initial spin-up window before
    #    comparing variability.
    #
    # Amplitude was checked for robustness across roughly 8e-8-3e-7 m/s
    # (holding mean == amplitude, forcing always >= 0): the shallow/deep
    # std(diff(.)) ratio stays close to ~2.8-2.9 with the shallow store
    # never pinned at its cap in the post-spin-up window; at 1.6e-5 m/s
    # (an amplitude tried during investigation of the original test) both
    # stores are pinned at saturation almost immediately and the comparison
    # becomes meaningless, confirming this must stay well below the
    # store's fill rate rather than being an arbitrarily large "stronger
    # signal is better" choice.
    comp = BucketHydrology(num_layers=2, soil_moisture_max=0.15,
                           deep_soil_moisture_max=0.5,
                           moisture_diffusion_timescale=86400.0,
                           deep_drainage_timescale=30 * 86400.0)
    state = get_default_state([comp], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["downwelling_shortwave_flux_in_air"].values[:] = 200.0
    state["downwelling_longwave_flux_in_air"].values[:] = 300.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 50.0
    state["upwelling_longwave_flux_in_air"].values[:] = 300.0
    state["soil_layer_thickness"].values[:] = 1.0
    state["surface_temperature"].values[:] = 295.0
    state["deep_soil_temperature"].values[:] = 290.0
    state["convective_precipitation_rate"].values[:] = 0.0
    # both layers start at the same 40% relative saturation
    state["lwe_thickness_of_soil_moisture_content"].values[:] = 0.4 * 0.15
    state["deep_soil_moisture_content"].values[:] = 0.4 * 0.5

    dt = timedelta(seconds=3600)
    amplitude = 1.5e-7  # m/s
    nsteps, spinup = 600, 300
    shallow, deep = [], []
    for step in range(nsteps):
        # sustained diurnal-ish oscillation, always non-negative
        precip = amplitude * (1 + np.sin(2 * np.pi * step / 24))
        state["stratiform_precipitation_rate"].values[:] = precip
        diag, new = comp(state, dt)
        state.update(new); state.update(diag)
        shallow.append(float(new["lwe_thickness_of_soil_moisture_content"].values[0]))
        deep.append(float(new["deep_soil_moisture_content"].values[0]))

    # discard the spin-up window; only compare quasi-equilibrium tracking
    shallow = np.array(shallow[spinup:])
    deep = np.array(deep[spinup:])

    # guard against the saturation-clipping artifact silently creeping back
    # in (e.g. via a future change to amplitude or drainage timescale)
    assert np.mean(shallow >= 0.15 - 1e-9) < 0.5, (
        "shallow store is pinned at its cap for most of the measurement "
        "window; the variability comparison below would reflect runoff "
        "clipping, not timescale separation"
    )

    # shallow store tracks the fast oscillation substantially more than the
    # low-pass-filtered deep store, once both are in quasi-equilibrium
    assert np.std(np.diff(shallow)) > 2 * np.std(np.diff(deep))

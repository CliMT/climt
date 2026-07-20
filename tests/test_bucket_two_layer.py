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

from datetime import timedelta
import numpy as np
from climt import (BucketHydrology, get_default_state)

def get_quantities(state):
    heat = state['surface_material_density'].values*state['soil_layer_thickness'].values*\
           state['surface_temperature'].values*state['heat_capacity_of_soil'].values
    moisture = state['lwe_thickness_of_soil_moisture_content'].values

    return heat, moisture

time_step = timedelta(seconds=1)
state = get_default_state([BucketHydrology()])

state['upwelling_shortwave_flux_in_air'].values[:] = 40.
state['upwelling_longwave_flux_in_air'].values[:] = 40.
state['downwelling_shortwave_flux_in_air'].values[:] = 40.
state['downwelling_longwave_flux_in_air'].values[:] = 40.
state['soil_layer_thickness'].values[:] = 1.
state['surface_specific_humidity'].values[:] = 1
state['northward_wind'].values[:] = 3
state['eastward_wind'].values[:] = 4
state['surface_temperature'].values[:] = 300
state['stratiform_precipitation_rate'].values[:] = 0.005
state['convective_precipitation_rate'].values[:] = 0.5
state['lwe_thickness_of_soil_moisture_content'].values[:] = 0.07


diag, new_state = BucketHydrology()(state, time_step)


surf_forcing = 0

if 'upwelling_shortwave_flux_in_air' in state:
    surf_forcing -= state['upwelling_shortwave_flux_in_air'].to_units('W/m^2')[0].values

if 'upwelling_longwave_flux_in_air' in state:
    surf_forcing -= state['upwelling_longwave_flux_in_air'].to_units('W/m^2')[0].values

if 'downwelling_shortwave_flux_in_air' in state:
    surf_forcing += state['downwelling_shortwave_flux_in_air'].to_units('W/m^2')[0].values

if 'downwelling_longwave_flux_in_air' in state:
    surf_forcing += state['downwelling_longwave_flux_in_air'].to_units('W/m^2')[0].values

if 'surface_upward_sensible_heat_flux' in diag:
    surf_forcing -= diag['surface_upward_sensible_heat_flux'].to_units('W/m^2').values

if 'surface_upward_latent_heat_flux' in diag:
    surf_forcing -= diag['surface_upward_latent_heat_flux'].to_units('W/m^2').values


mois_forcing = 0

if 'convective_precipitation_rate' in state:
    mois_forcing += state['convective_precipitation_rate'].to_units('m s^-1').values

if 'stratiform_precipitation_rate' in state:
    mois_forcing += state['stratiform_precipitation_rate'].to_units('m s^-1').values

if 'evaporation_rate' in diag:
    mois_forcing -= diag['evaporation_rate'].to_units('m s^-1').values

old_heat, old_moisture = get_quantities(state)

state.update(diag)
state.update(new_state)

new_heat, new_moisture = get_quantities(state)


#Test starts

def test_bucket_hydrology_heat_conservation():
    forcing_amount = surf_forcing * time_step.total_seconds()
    assert np.isclose(new_heat - old_heat, forcing_amount,
                      rtol=0, atol=1e-3)

def test_bucket_hydrology_moisture_conservation():
    forcing_amount = mois_forcing * time_step.total_seconds()
    assert np.isclose(new_moisture - old_moisture, forcing_amount,
                      rtol=0, atol=1e-3)

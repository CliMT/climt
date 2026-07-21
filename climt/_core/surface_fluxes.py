"""Shared bulk aerodynamic surface-flux formulae.

Matches the convention in BucketHydrology: sensible and latent fluxes use a
single bulk transfer coefficient times wind speed. This is a thin, tested
home for the formula so ocean/data-ocean and land paths agree.
"""
import numpy as np


def bulk_fluxes(wind_speed, surface_temperature, air_temperature,
                surface_specific_humidity, air_specific_humidity,
                air_density, bulk_coefficient=0.0011, latent_heat=2.5e6,
                beta=1.0):
    U = np.asarray(wind_speed)
    rho = np.asarray(air_density)
    potential_evap = bulk_coefficient * U * (
        np.asarray(surface_specific_humidity) - np.asarray(air_specific_humidity))
    evaporation_rate = beta * potential_evap
    latent = latent_heat * rho * evaporation_rate
    sensible = rho * bulk_coefficient * U * (
        np.asarray(surface_temperature) - np.asarray(air_temperature)) * 1004.0
    return {
        "sensible_heat_flux": sensible,
        "latent_heat_flux": latent,
        "evaporation_rate": evaporation_rate,
    }

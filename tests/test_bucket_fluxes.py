"""Dimensional correctness of BucketHydrology's bulk surface fluxes.

The sensible and latent heat fluxes must be W m^-2 and the evaporation rate a
liquid-water-equivalent depth rate (m s^-1), all built from one mass flux

    E = beta * rho * C * |v| * (q_s - q_a)      [kg m^-2 s^-1]

with ``rho`` the air density at the lowest model level. Before this was fixed
the density factors were missing: the sensible flux came out ~1200x too small,
the latent flux ~1.2x too small, and the evaporation rate ~830x too large.
"""
from datetime import timedelta

import numpy as np
from sympl import get_constant

from climt import BucketHydrology, get_default_state, get_grid


BULK_COEFFICIENT = 0.0011
LATENT_HEAT = 2260000.0


def _flux_state(soil_moisture=0.12, surface_humidity=0.02):
    """A windy, warm-surfaced column with a wet enough bucket that beta == 1."""
    comp = BucketHydrology(bulk_coefficient=BULK_COEFFICIENT,
                           specific_latent_heat_of_water=LATENT_HEAT)
    state = get_default_state([comp], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["soil_layer_thickness"].values[:] = 1.0
    state["northward_wind"].values[:] = 3.0
    state["eastward_wind"].values[:] = 4.0          # |v| = 5 m/s
    state["air_temperature"].values[:] = 290.0
    state["surface_temperature"].values[:] = 300.0  # dT = 10 K
    state["specific_humidity"].values[:] = 0.005
    state["surface_specific_humidity"].values[:] = surface_humidity
    state["lwe_thickness_of_soil_moisture_content"].values[:] = soil_moisture
    return comp, state


def _expected_terms(state):
    """(rho, |v|, dT, dq) at the lowest model level, from the state itself."""
    rd = get_constant("gas_constant_of_dry_air", "J kg^-1 K^-1")
    p = np.asarray(state["air_pressure"].values).ravel()[0]
    t_air = np.asarray(state["air_temperature"].values).ravel()[0]
    rho = p / (rd * t_air)
    wind = 5.0
    d_t = (np.asarray(state["surface_temperature"].values).ravel()[0] - t_air)
    d_q = (np.asarray(state["surface_specific_humidity"].values).ravel()[0]
           - np.asarray(state["specific_humidity"].values).ravel()[0])
    return rho, wind, d_t, d_q


def _scalar(field):
    return float(np.asarray(field.values).ravel()[0])


def test_air_pressure_is_an_input():
    """rho is derived from the ideal gas law, so pressure must be declared."""
    assert "air_pressure" in BucketHydrology().input_properties


def test_sensible_heat_flux_matches_the_bulk_aerodynamic_formula():
    comp, state = _flux_state()
    rho, wind, d_t, _ = _expected_terms(state)
    cp = get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                      "J kg^-1 K^-1")

    diagnostics, _ = comp(state, timedelta(seconds=1))

    expected = rho * cp * BULK_COEFFICIENT * wind * d_t
    assert np.isclose(_scalar(diagnostics["surface_upward_sensible_heat_flux"]),
                      expected, rtol=1e-10)


def test_sensible_heat_flux_has_a_physical_magnitude():
    """A 10 K land-air contrast at 5 m/s is tens of W m^-2, not hundredths."""
    comp, state = _flux_state()
    diagnostics, _ = comp(state, timedelta(seconds=1))
    assert 20.0 < _scalar(diagnostics["surface_upward_sensible_heat_flux"]) < 200.0


def test_latent_heat_flux_matches_the_bulk_aerodynamic_formula():
    comp, state = _flux_state()
    rho, wind, _, d_q = _expected_terms(state)

    diagnostics, _ = comp(state, timedelta(seconds=1))

    expected = LATENT_HEAT * rho * BULK_COEFFICIENT * wind * d_q
    assert np.isclose(_scalar(diagnostics["surface_upward_latent_heat_flux"]),
                      expected, rtol=1e-10)


def test_evaporation_rate_is_the_latent_flux_as_a_water_depth_rate():
    """evaporation_rate is declared m s^-1, so it is LH / (L * rho_water)."""
    comp, state = _flux_state()
    rho_water = get_constant("density_of_liquid_water", "kg m^-3")

    diagnostics, _ = comp(state, timedelta(seconds=1))

    expected = (_scalar(diagnostics["surface_upward_latent_heat_flux"])
                / (LATENT_HEAT * rho_water))
    assert np.isclose(_scalar(diagnostics["evaporation_rate"]), expected,
                      rtol=1e-10)


def test_evaporation_rate_has_a_physical_magnitude():
    """Millimetres per day, not metres per day.

    The unfixed code returned ~2400 mm/day here, which drains the 0.15 m
    bucket in about 90 minutes.
    """
    comp, state = _flux_state()
    diagnostics, _ = comp(state, timedelta(seconds=1))
    mm_per_day = _scalar(diagnostics["evaporation_rate"]) * 86400.0 * 1000.0
    assert 0.1 < mm_per_day < 50.0


def test_dry_bucket_shuts_off_evaporation_but_not_sensible_heat():
    """beta scales the moisture path only; the dry-surface limit."""
    comp, state = _flux_state(soil_moisture=0.0)
    diagnostics, _ = comp(state, timedelta(seconds=1))
    assert _scalar(diagnostics["evaporation_rate"]) == 0.0
    assert _scalar(diagnostics["surface_upward_latent_heat_flux"]) == 0.0
    assert _scalar(diagnostics["surface_upward_sensible_heat_flux"]) > 20.0


def test_fluxes_reverse_sign_over_a_surface_colder_than_the_air():
    comp, state = _flux_state()
    state["surface_temperature"].values[:] = 280.0     # 10 K colder than air
    state["surface_specific_humidity"].values[:] = 0.001   # drier than the air
    diagnostics, _ = comp(state, timedelta(seconds=1))
    assert _scalar(diagnostics["surface_upward_sensible_heat_flux"]) < 0.0
    assert _scalar(diagnostics["surface_upward_latent_heat_flux"]) < 0.0

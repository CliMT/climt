"""Surface-flux behaviour of climt.SimpleBoundaryLayer.

Covers the three ``surface_fluxes`` modes ('bulk', 'external', None), the
tridiagonal surface boundary condition that implements them, and the column
budgets they must close.
"""
from datetime import timedelta

import numpy as np
import pytest
from sympl import get_constant

import climt
from climt._components.simple_boundary_layer.component import _diffuse_profile


# ---------------------------------------------------------------- fixtures

def _test_column():
    """A small, well-conditioned 4-layer column for _diffuse_profile tests."""
    p_int = np.array([1.0e5, 9.0e4, 8.0e4, 7.0e4, 6.0e4])
    p = 0.5 * (p_int[:-1] + p_int[1:])
    rho = np.array([1.0, 0.9, 0.8])      # length num_layers - 1
    diff = np.array([10.0, 8.0, 6.0])    # length num_layers - 1
    return p, p_int, rho, diff


# -------------------------------------------- _diffuse_profile boundary terms

def test_diffuse_profile_zero_boundary_terms_reproduce_no_flux():
    """Both scalars zero must reproduce the pre-existing no-flux solve."""
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              0.0, 0.0)
    dp = p_int[:-1] - p_int[1:]
    assert np.isclose(np.sum(result * dp), np.sum(profile * dp), rtol=1e-12)


def test_diffuse_profile_does_not_mutate_profile():
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    before = profile.copy()
    _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665, 0.5, 140.0)
    assert np.array_equal(profile, before)


def test_diffuse_profile_bulk_term_matches_analytic_two_layer():
    """With diff == 0 the system is diagonal, so layer 0 has a closed form."""
    p, p_int, rho, _ = _test_column()
    diff = np.zeros(3)
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    beta = 0.25
    surface_value = 300.0
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              beta, beta * surface_value)
    expected0 = (profile[0] + beta * surface_value) / (1.0 + beta)
    assert np.isclose(result[0], expected0, rtol=1e-12)
    assert np.allclose(result[1:], profile[1:], rtol=1e-12)


def test_diffuse_profile_neumann_source_adds_exact_column_amount():
    """A pure Neumann source changes the mass-weighted integral by dp0 * S."""
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    source = 3.0
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              0.0, source)
    dp = p_int[:-1] - p_int[1:]
    change = np.sum(result * dp) - np.sum(profile * dp)
    assert np.isclose(change, dp[0] * source, rtol=1e-10)


# ------------------------------------------------------------------ helpers

def _constants():
    return (
        get_constant('gravitational_acceleration', 'm s^-2'),
        get_constant(
            'heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1'
        ),
        get_constant('latent_heat_of_condensation', 'J kg^-1'),
    )


def _column_state(component, nz=30):
    return climt.get_default_state(
        [component], grid_state=climt.get_grid(nx=None, ny=None, nz=nz)
    )


def _layer_mass(state):
    """dp for each mid-level of a single column (Pa)."""
    p_int = np.asarray(state['air_pressure_on_interface_levels'])
    return p_int[:-1] - p_int[1:]


def _column_budgets(state, dp):
    """(column enthalpy J/m2, column water kg/m2) for a single column."""
    g, cp, _ = _constants()
    enthalpy = cp * np.sum(np.asarray(state['air_temperature']) * dp) / g
    water = np.sum(np.asarray(state['specific_humidity']) * dp) / g
    return enthalpy, water


def _column_momentum(state, dp, name):
    """Column-integrated momentum of one wind component, kg m^-1 s^-1."""
    g, _, _ = _constants()
    return np.sum(np.asarray(state[name]) * dp) / g


def _scalar(field):
    return float(np.asarray(field).ravel()[0])


def _warm_wet_surface_state(component, nz=30):
    """320 K saturated-ish surface under a 250 K, bone-dry, windy atmosphere.

    This is the empirical case from the design doc that demonstrated the bug.
    """
    state = _column_state(component, nz=nz)
    np.asarray(state['air_temperature'])[:] = 250.0
    np.asarray(state['specific_humidity'])[:] = 0.0
    np.asarray(state['eastward_wind'])[:] = 10.0
    np.asarray(state['northward_wind'])[:] = 0.0
    np.asarray(state['surface_temperature'])[:] = 320.0
    np.asarray(state['surface_specific_humidity'])[:] = 0.02
    return state


# --------------------------------------------------------------------- API

def test_invalid_surface_fluxes_mode_raises():
    with pytest.raises(ValueError, match='surface_fluxes'):
        climt.SimpleBoundaryLayer(surface_fluxes='bogus')


def test_bulk_is_the_default_mode():
    assert climt.SimpleBoundaryLayer()._surface_fluxes == 'bulk'


def test_bulk_mode_diagnoses_surface_fluxes():
    c = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    assert set(c.diagnostic_properties) == {
        'northward_wind_stress',
        'eastward_wind_stress',
        'boundary_layer_height',
        'surface_upward_sensible_heat_flux',
        'surface_upward_latent_heat_flux',
    }
    assert (
        c.diagnostic_properties['surface_upward_sensible_heat_flux']['units']
        == 'W m^-2'
    )
    # bulk mode must NOT require the fluxes as inputs
    assert 'surface_upward_sensible_heat_flux' not in c.input_properties


def test_none_mode_has_no_flux_diagnostics():
    c = climt.SimpleBoundaryLayer(surface_fluxes=None)
    assert set(c.diagnostic_properties) == {
        'northward_wind_stress',
        'eastward_wind_stress',
        'boundary_layer_height',
    }
    assert 'surface_upward_sensible_heat_flux' not in c.input_properties


# ------------------------------------------------------------ bulk physics

def test_bulk_column_budget_closes_against_diagnosed_fluxes():
    """Column enthalpy/water change == diagnosed surface flux x dt, exactly."""
    _, _, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)

    enthalpy_before, water_before = _column_budgets(state, dp)
    timestep = timedelta(seconds=1200)
    diagnostics, new_state = component(state, timestep=timestep)
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    sh = _scalar(diagnostics['surface_upward_sensible_heat_flux'])
    lh = _scalar(diagnostics['surface_upward_latent_heat_flux'])
    dt = timestep.total_seconds()

    assert np.isclose(enthalpy_after - enthalpy_before, sh * dt, rtol=1e-10)
    assert np.isclose(water_after - water_before, lh * dt / lv, rtol=1e-10)


def test_bulk_momentum_budget_closes_against_the_stress_diagnostic():
    """Column momentum change == -wind stress x dt, exactly.

    The pre-reconciliation diagnostic (built from the pre-step *interface*
    wind) could not have passed this.
    """
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    np.asarray(state['northward_wind'])[:] = -4.0
    dp = _layer_mass(state)

    timestep = timedelta(seconds=1200)
    dt = timestep.total_seconds()
    before = {
        name: _column_momentum(state, dp, name)
        for name in ('eastward_wind', 'northward_wind')
    }
    diagnostics, new_state = component(state, timestep=timestep)

    for name, stress in (('eastward_wind', 'eastward_wind_stress'),
                         ('northward_wind', 'northward_wind_stress')):
        change = _column_momentum(new_state, dp, name) - before[name]
        assert np.isclose(
            change, -_scalar(diagnostics[stress]) * dt, rtol=1e-10
        ), name


def test_bulk_wind_stress_uses_the_post_solve_layer_zero_wind():
    """Guard the reconciliation against a revert to the old formula.

    kappa is recovered from the momentum budget, which the solver satisfies
    regardless of how the diagnostic is written -- so this is not circular.
    The old formula used the pre-step interface wind
    ``0.5 * (u0_old + u1_old)``; under shear that is a different number.
    """
    g, _, _ = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    # sheared profile: mid-level 0 is the surface, so this is weak wind at the
    # bottom and strong aloft, making u0 and the 0/1 interface clearly differ
    u = np.asarray(state['eastward_wind'])
    u[:] = np.linspace(2.0, 30.0, u.shape[0]).reshape(
        (-1,) + (1,) * (u.ndim - 1)
    )
    u_old_0, u_old_1 = _scalar(u[0]), _scalar(u[1])
    u_interface_old = 0.5 * (u_old_0 + u_old_1)

    dp = _layer_mass(state)
    before = _column_momentum(state, dp, 'eastward_wind')
    timestep = timedelta(seconds=1200)
    dt = timestep.total_seconds()
    diagnostics, new_state = component(state, timestep=timestep)

    tau = _scalar(diagnostics['eastward_wind_stress'])
    u_new_0 = _scalar(np.asarray(new_state['eastward_wind'])[0])
    momentum_lost = before - _column_momentum(new_state, dp, 'eastward_wind')
    kappa = momentum_lost / (u_new_0 * dt)

    assert tau > 0.0
    assert np.isclose(tau, kappa * u_new_0, rtol=1e-10)
    # the two candidate winds must be distinguishable, so that the assertion
    # below actually discriminates (which one is larger depends on how drag
    # and shear-driven diffusion compete -- do not assume an ordering)
    assert not np.isclose(u_interface_old, u_new_0, rtol=1e-2)
    # the old form would have reported a visibly different number
    assert not np.isclose(tau, kappa * u_interface_old, rtol=1e-3)


def test_none_mode_stress_is_diagnosed_but_not_applied():
    component = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)
    before = _column_momentum(state, dp, 'eastward_wind')

    diagnostics, new_state = component(state, timestep=timedelta(seconds=1200))

    # advisory: a real stress is reported ...
    assert _scalar(diagnostics['eastward_wind_stress']) > 0.0
    # ... but no momentum left the column
    after = _column_momentum(new_state, dp, 'eastward_wind')
    assert np.isclose(after, before, rtol=1e-12)


def test_bulk_mode_warms_moistens_and_decelerates():
    """The design doc's empirical case, with the assertions inverted."""
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)
    u_before = np.asarray(state['eastward_wind'])[0]

    for _ in range(20):
        _, new_state = component(state, timestep=timedelta(seconds=1200))
        state.update(new_state)

    enthalpy_after, water_after = _column_budgets(state, dp)
    assert enthalpy_after > enthalpy_before
    assert water_after > water_before
    assert np.asarray(state['air_temperature'])[0] > 250.0
    assert np.asarray(state['air_temperature'])[0] < 320.0
    assert abs(np.asarray(state['eastward_wind'])[0]) < abs(u_before)
    assert np.all(np.isfinite(np.asarray(state['air_temperature'])))


def test_bulk_mode_cools_column_over_a_cold_surface():
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _column_state(component)
    np.asarray(state['air_temperature'])[:] = 290.0
    np.asarray(state['surface_temperature'])[:] = 250.0
    np.asarray(state['eastward_wind'])[:] = 10.0
    dp = _layer_mass(state)
    enthalpy_before, _ = _column_budgets(state, dp)

    _, new_state = component(state, timestep=timedelta(seconds=1200))
    enthalpy_after, _ = _column_budgets(new_state, dp)
    assert enthalpy_after < enthalpy_before


def test_none_mode_conserves_over_a_warm_wet_surface():
    """The regression guard: opting out must restore exact conservation."""
    component = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)

    _, new_state = component(state, timestep=timedelta(seconds=1200))
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    assert np.isclose(enthalpy_after, enthalpy_before, rtol=1e-12)
    assert np.isclose(water_after, water_before, rtol=1e-12)
    assert np.asarray(new_state['air_temperature'])[0] == pytest.approx(250.0)


# -------------------------------------------------------- external physics

def test_external_mode_requires_the_flux_inputs():
    c = climt.SimpleBoundaryLayer(surface_fluxes='external')
    assert 'surface_upward_sensible_heat_flux' in c.input_properties
    assert 'surface_upward_latent_heat_flux' in c.input_properties
    assert (
        c.input_properties['surface_upward_latent_heat_flux']['units']
        == 'W m^-2'
    )
    # they are inputs, not diagnostics -- the two must not collide
    assert 'surface_upward_sensible_heat_flux' not in c.diagnostic_properties


def test_external_mode_applies_prescribed_fluxes_exactly():
    _, _, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _column_state(component)
    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = 50.0
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = 100.0

    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)
    timestep = timedelta(seconds=1200)
    _, new_state = component(state, timestep=timestep)
    enthalpy_after, water_after = _column_budgets(new_state, dp)
    dt = timestep.total_seconds()

    assert np.isclose(enthalpy_after - enthalpy_before, 50.0 * dt, rtol=1e-10)
    assert np.isclose(
        water_after - water_before, 100.0 * dt / lv, rtol=1e-10
    )


def test_external_mode_zero_fluxes_conserve():
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _column_state(component)
    np.asarray(state['specific_humidity'])[:] = 0.005
    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)

    _, new_state = component(state, timestep=timedelta(seconds=1200))
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    assert np.isclose(enthalpy_after, enthalpy_before, rtol=1e-12)
    assert np.isclose(water_after, water_before, rtol=1e-12)


def test_external_mode_still_applies_bulk_momentum_drag():
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _warm_wet_surface_state(component)
    # fluxes stay zero; only drag should act
    u_before = np.asarray(state['eastward_wind'])[0]
    _, new_state = component(state, timestep=timedelta(seconds=1200))
    assert abs(np.asarray(new_state['eastward_wind'])[0]) < abs(u_before)


def test_external_momentum_budget_closes_against_the_stress_diagnostic():
    """Momentum is bulk-internal in external mode, so its budget closes too."""
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _warm_wet_surface_state(component)
    np.asarray(state['northward_wind'])[:] = -4.0
    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = 50.0
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = 100.0
    dp = _layer_mass(state)

    timestep = timedelta(seconds=1200)
    dt = timestep.total_seconds()
    before = {
        name: _column_momentum(state, dp, name)
        for name in ('eastward_wind', 'northward_wind')
    }
    diagnostics, new_state = component(state, timestep=timestep)

    for name, stress in (('eastward_wind', 'eastward_wind_stress'),
                         ('northward_wind', 'northward_wind_stress')):
        change = _column_momentum(new_state, dp, name) - before[name]
        assert np.isclose(
            change, -_scalar(diagnostics[stress]) * dt, rtol=1e-10
        ), name


def test_bulk_and_external_agree_when_fed_the_bulk_fluxes():
    bulk = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    external = climt.SimpleBoundaryLayer(surface_fluxes='external')
    # build the state from the external component so it carries the flux
    # fields; bulk simply ignores the extras.
    state = _warm_wet_surface_state(external)

    diagnostics, bulk_new = bulk(state, timestep=timedelta(seconds=1200))

    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = _scalar(
        diagnostics['surface_upward_sensible_heat_flux']
    )
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = _scalar(
        diagnostics['surface_upward_latent_heat_flux']
    )
    external_diagnostics, external_new = external(
        state, timestep=timedelta(seconds=1200)
    )

    for name in ('air_temperature', 'specific_humidity',
                 'eastward_wind', 'northward_wind'):
        assert np.allclose(
            np.asarray(bulk_new[name]), np.asarray(external_new[name]),
            rtol=1e-10, atol=1e-12,
        ), name
    # momentum is bulk-internal in both modes, so the stresses match too
    for name in ('eastward_wind_stress', 'northward_wind_stress'):
        assert np.isclose(
            _scalar(diagnostics[name]),
            _scalar(external_diagnostics[name]),
            rtol=1e-10,
        ), name


def test_external_mode_applies_downward_fluxes_exactly():
    """Negative (downward) prescribed fluxes must be applied with equal exactness."""
    _, _, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _column_state(component)
    # start moist enough that a downward latent flux cannot drive q negative
    np.asarray(state['specific_humidity'])[:] = 0.005
    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = -30.0
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = -20.0

    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)
    timestep = timedelta(seconds=1200)
    _, new_state = component(state, timestep=timestep)
    enthalpy_after, water_after = _column_budgets(new_state, dp)
    dt = timestep.total_seconds()

    assert np.isclose(enthalpy_after - enthalpy_before, -30.0 * dt, rtol=1e-10)
    assert np.isclose(water_after - water_before, -20.0 * dt / lv, rtol=1e-10)
    # the budget above is only meaningful while the humidity stays physical
    assert np.all(np.asarray(new_state['specific_humidity']) > 0.0)


def test_bulk_flux_diagnostics_are_negative_over_a_cold_dry_surface():
    """Both diagnosed bulk fluxes point downward when the surface is cold and dry."""
    _, _, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _column_state(component)
    np.asarray(state['air_temperature'])[:] = 290.0
    np.asarray(state['specific_humidity'])[:] = 0.010
    np.asarray(state['eastward_wind'])[:] = 10.0
    np.asarray(state['northward_wind'])[:] = 0.0
    np.asarray(state['surface_temperature'])[:] = 250.0
    np.asarray(state['surface_specific_humidity'])[:] = 0.001

    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)
    timestep = timedelta(seconds=1200)
    diagnostics, new_state = component(state, timestep=timestep)
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    sh = _scalar(diagnostics['surface_upward_sensible_heat_flux'])
    lh = _scalar(diagnostics['surface_upward_latent_heat_flux'])
    dt = timestep.total_seconds()

    assert sh < 0.0
    assert lh < 0.0
    # a self-consistent but inverted sign convention must not survive this
    assert np.isclose(enthalpy_after - enthalpy_before, sh * dt, rtol=1e-10)
    assert np.isclose(water_after - water_before, lh * dt / lv, rtol=1e-10)


# ---------------------------------------------------------- multi-column

def _grid_state(component, nx=4, ny=3, nz=20):
    return climt.get_default_state(
        [component], grid_state=climt.get_grid(nx=nx, ny=ny, nz=nz)
    )


def _per_column_layer_mass(state):
    """dp with shape (nz, ncol); columns flattened as the component flattens them."""
    p_int = np.asarray(state['air_pressure_on_interface_levels'])
    p_int = p_int.reshape(p_int.shape[0], -1)
    return p_int[:-1] - p_int[1:]


def _per_column_budgets(state, dp, name):
    """Mass-weighted column integral of ``state[name]``, one entry per column."""
    g, _, _ = _constants()
    field = np.asarray(state[name])
    field = field.reshape(field.shape[0], -1)
    return np.sum(field * dp, axis=0) / g


def test_external_mode_multi_column_applies_per_column_fluxes():
    """Each column gets exactly its own prescribed flux, not column zero's."""
    _, cp, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _grid_state(component)
    np.asarray(state['specific_humidity'])[:] = 0.005

    sh_field = np.asarray(state['surface_upward_sensible_heat_flux'])
    lh_field = np.asarray(state['surface_upward_latent_heat_flux'])
    ncol = sh_field.size
    sh_col = np.linspace(-40.0, 80.0, ncol)
    lh_col = np.linspace(120.0, -30.0, ncol)
    sh_field[:] = sh_col.reshape(sh_field.shape)
    lh_field[:] = lh_col.reshape(lh_field.shape)
    # make the columns differ in kappa as well as in prescribed flux
    ts = np.asarray(state['surface_temperature'])
    ts[:] = np.linspace(275.0, 305.0, ncol).reshape(ts.shape)
    u = np.asarray(state['eastward_wind'])
    u[:] = np.linspace(2.0, 15.0, ncol).reshape((1,) + u.shape[1:])

    dp = _per_column_layer_mass(state)
    enthalpy_before = cp * _per_column_budgets(state, dp, 'air_temperature')
    water_before = _per_column_budgets(state, dp, 'specific_humidity')
    timestep = timedelta(seconds=1200)
    _, new_state = component(state, timestep=timestep)
    enthalpy_after = cp * _per_column_budgets(new_state, dp, 'air_temperature')
    water_after = _per_column_budgets(new_state, dp, 'specific_humidity')
    dt = timestep.total_seconds()

    assert np.allclose(enthalpy_after - enthalpy_before, sh_col * dt,
                       rtol=1e-10)
    assert np.allclose(water_after - water_before, lh_col * dt / lv, rtol=1e-10)
    assert np.all(np.asarray(new_state['specific_humidity']) > 0.0)


def test_none_mode_multi_column_conserves_every_column():
    """With surface exchange off, every column of a grid conserves exactly."""
    component = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = _grid_state(component)
    ts = np.asarray(state['surface_temperature'])
    ncol = ts.size
    ts[:] = np.linspace(260.0, 310.0, ncol).reshape(ts.shape)
    u = np.asarray(state['eastward_wind'])
    u[:] = np.linspace(2.0, 15.0, ncol).reshape((1,) + u.shape[1:])
    v = np.asarray(state['northward_wind'])
    v[:] = np.linspace(-8.0, 4.0, ncol).reshape((1,) + v.shape[1:])
    q = np.asarray(state['specific_humidity'])
    q[:] = np.linspace(0.002, 0.012, ncol).reshape((1,) + q.shape[1:])

    dp = _per_column_layer_mass(state)
    before = {
        name: _per_column_budgets(state, dp, name)
        for name in ('air_temperature', 'specific_humidity',
                     'northward_wind', 'eastward_wind')
    }
    _, new_state = component(state, timestep=timedelta(seconds=1200))

    for name, expected in before.items():
        after = _per_column_budgets(new_state, dp, name)
        assert np.allclose(after, expected, rtol=1e-12), name
        assert np.all(np.isfinite(np.asarray(new_state[name]))), name

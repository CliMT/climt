# -*- coding: utf-8 -*-
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax

jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import jax.numpy as jnp  # noqa: E402
from tests.jax_cases import emanuel_case, emanuel_loss  # noqa: E402


@pytest.mark.parametrize(
    "profile",
    [
        "unstable",
        pytest.param("moist", marks=pytest.mark.slow),
        pytest.param("stable", marks=pytest.mark.slow),
    ],
)
def test_emanuel_jax_differentiation_direct(profile):
    """Test that gradients through EmanuelConvectionPython.array_call are finite
    across convecting and non-convecting profiles, and non-zero when unstable."""
    component, state, timestep = emanuel_case(nlev=30, ncol=1, profile=profile, as_jax=True)
    grad_func = jax.grad(emanuel_loss(component, state, timestep))
    gradient = grad_func(state["air_temperature"])

    assert jnp.all(jnp.isfinite(gradient)), f"Gradient has non-finite entries for profile={profile!r}"
    if profile == "unstable":
        assert jnp.any(gradient != 0.0), "Gradient is zero on unstable profile."


def test_sat_vap_pressure_jax_grad_finite_at_extreme_cold():
    from climt._components.emanuel.kernel_jax import _sat_vap_pressure_jax

    temps = jnp.array([28.0, 29.65, 180.0, 273.15, 300.0], dtype=jnp.float64)
    grad = jax.grad(lambda t: jnp.sum(_sat_vap_pressure_jax(t, 0)))(temps)
    assert jnp.all(jnp.isfinite(grad)), f"Non-finite gradient in _sat_vap_pressure_jax: {grad}"

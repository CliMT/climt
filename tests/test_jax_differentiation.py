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
    """Test finite gradients through EmanuelConvectionPython.array_call."""
    component, state, timestep = emanuel_case(nlev=30, ncol=1, profile=profile, as_jax=True)
    loss_fn = emanuel_loss(component, state, timestep)
    gradient = jax.grad(loss_fn)(state["air_temperature"])

    assert jnp.all(jnp.isfinite(gradient)), f"Non-finite gradient for profile={profile}"
    if profile == "unstable":
        assert jnp.any(gradient != 0.0), "Gradient is zero for unstable profile."


def test_emanuel_jax_zero_humidity_gradient_is_finite():
    """Gradients w.r.t. T and q stay finite when specific_humidity is identically zero."""
    component, state, timestep = emanuel_case(nlev=30, ncol=1, profile="unstable", as_jax=True)
    T0 = state["air_temperature"]
    q0 = jnp.zeros_like(state["specific_humidity"])

    def loss(T, q):
        tendencies, _ = component.array_call(
            {**state, "air_temperature": T, "specific_humidity": q}, timestep
        )
        return jnp.sum(tendencies["air_temperature"] ** 2)

    gT, gq = jax.grad(loss, argnums=(0, 1))(T0, q0)
    assert jnp.all(jnp.isfinite(gT)), "Non-finite temperature gradient on q=0 column"
    assert jnp.all(jnp.isfinite(gq)), "Non-finite humidity gradient on q=0 column"

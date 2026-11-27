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

# -*- coding: utf-8 -*-
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax

jax.config.update("jax_enable_x64", True)

from datetime import timedelta  # noqa: E402
import jax.numpy as jnp  # noqa: E402
from climt import EmanuelConvectionPython, get_default_state, get_grid  # noqa: E402


@pytest.mark.slow
def test_emanuel_jax_differentiation_direct():
    """Test that we can compute gradients through EmanuelConvectionPython.array_call."""
    jax.config.update("jax_platform_name", "cpu")

    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    conv = EmanuelConvectionPython()
    timestep = timedelta(minutes=15)
    state_sympl = get_default_state([conv], grid_state=grid)

    # Create an unstable profile
    temp = jnp.linspace(310, 250, nlev).reshape(nlev, ncol)
    q = jnp.zeros((nlev, ncol)).at[0:10].set(0.015)

    jax_state = {
        "air_temperature": temp,
        "specific_humidity": q,
        "air_pressure": jnp.array(state_sympl["air_pressure"].values / 100.0).reshape(nlev, ncol),
        "air_pressure_on_interface_levels": jnp.array(
            state_sympl["air_pressure_on_interface_levels"].values / 100.0
        ).reshape(nlev + 1, ncol),
        "eastward_wind": jnp.zeros((nlev, ncol)),
        "northward_wind": jnp.zeros((nlev, ncol)),
    }

    def compute_loss(temp_values):
        current_state = {**jax_state, "air_temperature": temp_values}
        tendencies, _ = conv.array_call(current_state, timestep)
        return jnp.sum(jnp.square(tendencies["air_temperature"]))

    initial_temp = jax_state["air_temperature"]
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_temp)

    assert jnp.any(gradient != 0.0), "Gradient is zero on unstable profile."
    assert not jnp.any(jnp.isnan(gradient)), "Gradient contains NaNs"

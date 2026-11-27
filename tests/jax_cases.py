# -*- coding: utf-8 -*-
from datetime import timedelta

import numpy as np

from climt import EmanuelConvectionPython, get_default_state, get_grid
from climt._core.util import numpy_version_of


def emanuel_case(nlev=30, ncol=1, profile="unstable", as_jax=True):
    """Return (component, state, timestep): an EmanuelConvectionPython, the
    array_call state dictionary (pressures in hPa, arrays (nlev, ncol)), and
    timedelta(minutes=15). profile is 'moist', 'dry', 'unstable' (as in
    tests/test_emanuel_jax_parity.py) or 'stable' (isothermal 250 K,
    specific humidity 1e-6)."""
    component = EmanuelConvectionPython()
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    sympl_state = get_default_state([component], grid_state=grid)

    temp = sympl_state["air_temperature"].values
    q = np.zeros_like(sympl_state["specific_humidity"].values)

    if profile == "moist":
        temp[:] = 300.0
        q[0:5, :, :] = 0.020
    elif profile == "dry":
        temp[:] = 300.0
        q[:] = 1e-6
    elif profile == "unstable":
        for i in range(nlev):
            temp[i, :, :] = 310.0 - (i * 2.0)
        q[0:10, :, :] = 0.015
    elif profile == "stable":
        temp[:] = 250.0
        q[:] = 1e-6
    else:
        raise ValueError(f"Unknown profile: {profile!r}")

    sympl_state["specific_humidity"].values[:] = q
    raw_state = numpy_version_of(sympl_state)

    if as_jax:
        import jax.numpy as jnp

        xp = jnp
    else:
        xp = np

    state = {
        "air_temperature": xp.array(raw_state["air_temperature"].reshape(nlev, ncol)),
        "specific_humidity": xp.array(raw_state["specific_humidity"].reshape(nlev, ncol)),
        "air_pressure": xp.array(raw_state["air_pressure"].reshape(nlev, ncol) / 100.0),
        "air_pressure_on_interface_levels": xp.array(
            raw_state["air_pressure_on_interface_levels"].reshape(nlev + 1, ncol) / 100.0
        ),
        "eastward_wind": xp.zeros((nlev, ncol)),
        "northward_wind": xp.zeros((nlev, ncol)),
    }
    return component, state, timedelta(minutes=15)


def emanuel_forward(component, state, timestep):
    """Function of air_temperature returning the temperature tendency array."""

    def forward(air_temperature):
        tendencies, _ = component.array_call(
            {**state, "air_temperature": air_temperature}, timestep
        )
        return tendencies["air_temperature"]

    return forward


def emanuel_loss(component, state, timestep):
    """Function of air_temperature returning sum(tendency**2); for jax.grad."""
    import jax.numpy as jnp

    forward = emanuel_forward(component, state, timestep)

    def loss(air_temperature):
        tendency = forward(air_temperature)
        return jnp.sum(jnp.square(tendency))

    return loss

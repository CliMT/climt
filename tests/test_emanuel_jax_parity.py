# tests/test_emanuel_jax_parity.py
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
from datetime import timedelta  # noqa: E402
from climt import EmanuelConvectionPython, get_default_state, get_grid  # noqa: E402
from climt._core.util import numpy_version_of  # noqa: E402


def create_test_state(nlev, ncol, moisture_type="moist"):
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvectionPython()], grid_state=grid)
    state["air_temperature"].values[:] = 300.0
    q = np.zeros_like(state["specific_humidity"].values)
    if moisture_type == "moist":
        q[0:5, :, :] = 0.020
    elif moisture_type == "dry":
        q[:] = 1e-6
    elif moisture_type == "unstable":
        for i in range(nlev):
            state["air_temperature"].values[i, :, :] = 310.0 - (i * 2.0)
        q[0:10, :, :] = 0.015
    state["specific_humidity"].values[:] = q
    return state


def make_column_state(raw_state, nlev, ncol, as_jax=False):
    xp = jnp if as_jax else np
    return {
        "air_temperature": xp.array(raw_state["air_temperature"].reshape(nlev, ncol)),
        "specific_humidity": xp.array(raw_state["specific_humidity"].reshape(nlev, ncol)),
        "air_pressure": xp.array(raw_state["air_pressure"].reshape(nlev, ncol) / 100.0),
        "air_pressure_on_interface_levels": xp.array(
            raw_state["air_pressure_on_interface_levels"].reshape(nlev + 1, ncol) / 100.0
        ),
        "eastward_wind": xp.zeros((nlev, ncol)),
        "northward_wind": xp.zeros((nlev, ncol)),
    }


@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("moisture_type", ["moist", "dry", "unstable"])
def test_jax_vs_numpy_parity(ncol, moisture_type):
    nlev = 30
    state = create_test_state(nlev, ncol, moisture_type)
    timestep = timedelta(minutes=15)
    raw_state = numpy_version_of(state)
    conv = EmanuelConvectionPython()

    np_state = make_column_state(raw_state, nlev, ncol, as_jax=False)
    jax_state = make_column_state(raw_state, nlev, ncol, as_jax=True)

    tend_np, diag_np = conv.array_call(np_state, timestep)
    tend_jax, diag_jax = conv.array_call(jax_state, timestep)

    for key in tend_np:
        np.testing.assert_allclose(
            np.asarray(tend_jax[key]), tend_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Tendency '{key}' (ncol={ncol}, {moisture_type})",
        )
    for key in ["convective_precipitation_rate", "cloud_base_mass_flux",
                 "atmosphere_convective_available_potential_energy"]:
        np.testing.assert_allclose(
            np.asarray(diag_jax[key]), diag_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Diagnostic '{key}' (ncol={ncol}, {moisture_type})",
        )

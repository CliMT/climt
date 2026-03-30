# -*- coding: utf-8 -*-
"""Parity tests comparing EmanuelConvectionPythonV3 against V1 (the reference implementation).

V3 has minor algorithmic differences from V1 (vectorised kernel, NamedTuple params),
so tolerances are relaxed to atol=1e-3, rtol=1e-3.

V2 is intentionally NOT imported here to avoid triggering the V2 set_item Numba bug.
"""

from datetime import timedelta

import numpy as np
import pytest

from climt import EmanuelConvection, get_default_state, get_grid
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPythonV3
from climt._core.util import numpy_version_of


def create_test_state(nlev, ncol, moisture_type="moist"):
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
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


@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("moisture_type", ["moist", "dry", "unstable"])
def test_emanuel_v3_parity(ncol, moisture_type):
    """V3 tendencies and diagnostics must agree with V1 within atol=rtol=1e-3."""
    nlev = 30
    state = create_test_state(nlev, ncol, moisture_type)
    timestep = timedelta(minutes=15)

    raw_state = numpy_version_of(state)

    python_state = {
        "air_temperature": raw_state["air_temperature"].reshape(nlev, ncol),
        "specific_humidity": raw_state["specific_humidity"].reshape(nlev, ncol),
        "air_pressure": (raw_state["air_pressure"] / 100.0).reshape(nlev, ncol),
        "air_pressure_on_interface_levels": (
            raw_state["air_pressure_on_interface_levels"] / 100.0
        ).reshape(nlev + 1, ncol),
        "eastward_wind": raw_state["eastward_wind"].reshape(nlev, ncol),
        "northward_wind": raw_state["northward_wind"].reshape(nlev, ncol),
    }

    conv_v1 = EmanuelConvectionPython()
    conv_v3 = EmanuelConvectionPythonV3()

    tend_v1, diag_v1 = conv_v1.array_call(python_state, timestep)
    tend_v3, diag_v3 = conv_v3.array_call(python_state, timestep)

    for key in tend_v1:
        np.testing.assert_allclose(
            tend_v1[key],
            tend_v3[key],
            atol=1e-3,
            rtol=1e-3,
            err_msg=f"Tendency parity failed for '{key}' (ncol={ncol}, type={moisture_type})",
        )

    for key in diag_v1:
        np.testing.assert_allclose(
            diag_v1[key],
            diag_v3[key],
            atol=1e-3,
            rtol=1e-3,
            err_msg=f"Diagnostic parity failed for '{key}' (ncol={ncol}, type={moisture_type})",
        )

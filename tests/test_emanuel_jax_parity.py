# tests/test_emanuel_jax_parity.py
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import numpy as np  # noqa: E402
from tests.jax_cases import emanuel_case  # noqa: E402


@pytest.mark.parametrize("nlev", [30, 45])
@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("profile", ["moist", "dry", "unstable", "stable"])
def test_jax_vs_numpy_parity(nlev, ncol, profile):
    conv_np, np_state, timestep = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=False)
    conv_jax, jax_state, _ = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=True)

    tend_np, diag_np = conv_np.array_call(np_state, timestep)
    tend_jax, diag_jax = conv_jax.array_call(jax_state, timestep)

    for key in tend_np:
        np.testing.assert_allclose(
            np.asarray(tend_jax[key]), tend_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Tendency '{key}' (nlev={nlev}, ncol={ncol}, {profile})",
        )
    for key in [
        "convective_precipitation_rate",
        "cloud_base_mass_flux",
        "atmosphere_convective_available_potential_energy",
    ]:
        np.testing.assert_allclose(
            np.asarray(diag_jax[key]), diag_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Diagnostic '{key}' (nlev={nlev}, ncol={ncol}, {profile})",
        )

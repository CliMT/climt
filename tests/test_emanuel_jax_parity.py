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
    conv, np_state, timestep = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=False)
    _, jax_state, _ = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=True)

    tend_np, diag_np = conv.array_call(np_state, timestep)
    tend_jax, diag_jax = conv.array_call(jax_state, timestep)

    for key in tend_np:
        np.testing.assert_allclose(
            np.asarray(tend_jax[key]),
            tend_np[key],
            atol=1e-10,
            rtol=1e-10,
            err_msg=f"Tendency '{key}' (nlev={nlev}, ncol={ncol}, {profile})",
        )
    for key in diag_np:
        if key == "convective_state":
            assert diag_jax[key].dtype == np.int32
            np.testing.assert_array_equal(
                np.asarray(diag_jax[key]),
                diag_np[key],
                err_msg=f"Diagnostic '{key}' (nlev={nlev}, ncol={ncol}, {profile})",
            )
        else:
            np.testing.assert_allclose(
                np.asarray(diag_jax[key]),
                diag_np[key],
                atol=1e-10,
                rtol=1e-10,
                err_msg=f"Diagnostic '{key}' (nlev={nlev}, ncol={ncol}, {profile})",
            )


def test_emanuel_jax_exit5_preserves_outcape_and_iflag():
    """When exit1..exit4 pass (active parcel with TVP[ICB] > TV[ICB] - DTMAX)
    and positive CAPE exists aloft, yet DTMA <= 0 keeps CBMF == 0 from CBMFOLD == 0
    (exit5), NumPy returns IFLAG == 1 and OUTCAPE == CAPE > 0."""
    conv = __import__("climt").EmanuelConvectionPython(DTMAX=0.85, ALPHA=0.1)
    _, np_state, timestep = emanuel_case(nlev=30, ncol=1, profile="unstable", as_jax=False)
    # Cool only the sub-cloud/LCL transition slightly so TVP[ICB] > TV[ICB] - DTMAX
    # passes exit4 while DTMA <= 0 keeps CBMF == 0 (exit5), with unstable CAPE aloft.
    np_state["air_temperature"] = np_state["air_temperature"].copy()
    np_state["air_temperature"][1:4, 0] += 4.5
    np_state["specific_humidity"] = np_state["specific_humidity"].copy()
    np_state["specific_humidity"][0, 0] = 0.014
    np_state["specific_humidity"][1:10, 0] = 0.004

    import jax.numpy as jnp

    jax_state = {k: jnp.asarray(v) for k, v in np_state.items()}
    tend_np, diag_np = conv.array_call(np_state, timestep)
    tend_jax, diag_jax = conv.array_call(jax_state, timestep)

    assert diag_jax["convective_state"].dtype == np.int32
    np.testing.assert_array_equal(np.asarray(diag_jax["convective_state"]), diag_np["convective_state"])
    np.testing.assert_allclose(
        np.asarray(diag_jax["atmosphere_convective_available_potential_energy"]),
        diag_np["atmosphere_convective_available_potential_energy"],
        rtol=1e-10,
        atol=1e-10,
    )
    for key in tend_np:
        np.testing.assert_allclose(np.asarray(tend_jax[key]), tend_np[key], rtol=1e-10, atol=1e-10)

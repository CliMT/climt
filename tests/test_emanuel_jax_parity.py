# tests/test_emanuel_jax_parity.py
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402

from tests.jax_cases import emanuel_case  # noqa: E402


def _assert_parity(tend_jax, diag_jax, tend_np, diag_np, label):
    for key in tend_np:
        np.testing.assert_allclose(
            np.asarray(tend_jax[key]), tend_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Tendency '{key}' ({label})",
        )
    for key in diag_np:
        if key == "convective_state":
            assert diag_jax[key].dtype == np.int32
            np.testing.assert_array_equal(
                np.asarray(diag_jax[key]), diag_np[key], err_msg=f"Diagnostic '{key}' ({label})",
            )
        else:
            np.testing.assert_allclose(
                np.asarray(diag_jax[key]), diag_np[key], atol=1e-10, rtol=1e-10,
                err_msg=f"Diagnostic '{key}' ({label})",
            )


@pytest.mark.parametrize("nlev", [30, 45])
@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("profile", ["moist", "dry", "unstable", "stable"])
def test_jax_vs_numpy_parity(nlev, ncol, profile):
    conv_np, np_state, timestep = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=False)
    conv_jax, jax_state, _ = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=True)

    tend_np, diag_np = conv_np.array_call(np_state, timestep)
    tend_jax, diag_jax = conv_jax.array_call(jax_state, timestep)
    _assert_parity(tend_jax, diag_jax, tend_np, diag_np, f"nlev={nlev}, ncol={ncol}, {profile}")


@pytest.mark.parametrize("profile", ["moist", "unstable"])
def test_jax_vs_numpy_parity_with_winds_and_carried_mass_flux(profile):
    """Five steps with nonzero winds, feeding each step's cloud_base_mass_flux
    into the next. The single-call cases start from zero winds and zero mass
    flux, so they never compare momentum transport, the CBMFOLD > 0 branch of
    the closure, or (in most cases) any precipitation."""
    nlev, ncol = 30, 3
    conv_np, np_state, timestep = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=False)
    conv_jax, jax_state, _ = emanuel_case(nlev=nlev, ncol=ncol, profile=profile, as_jax=True)
    rng = np.random.default_rng(0)
    for name in ("eastward_wind", "northward_wind"):
        wind = rng.normal(0.0, 10.0, (nlev, ncol))
        np_state[name] = wind.copy()
        jax_state[name] = jnp.asarray(wind)

    for step in range(5):
        tend_np, diag_np = conv_np.array_call(dict(np_state), timestep)
        tend_jax, diag_jax = conv_jax.array_call(dict(jax_state), timestep)
        _assert_parity(tend_jax, diag_jax, tend_np, diag_np, f"{profile}, step {step}")
        cbmf = np.asarray(diag_np["cloud_base_mass_flux"]).copy()
        np_state["cloud_base_mass_flux"] = cbmf
        jax_state["cloud_base_mass_flux"] = jnp.asarray(cbmf)

    assert np.any(cbmf > 0.0), "the case never built up a cloud-base mass flux"


# Found by random search: exit 5 (IFLAG == 1, CBMF == 0 from CBMFOLD == 0)
# is rare, and rarer still with nonzero CAPE (1 in ~8000 random profiles).
_EXIT5_TEMPERATURE_OFFSET = [
    1.03, 3.58, 5.91, 5.45, 4.47, 5.26, 6.45, 5.13, 3.9, 2.14, 1.29, 1.39, -0.86, -1.59, -2.3,
    -1.25, -3.07, -3.27, -0.94, -0.69, 0.03, -0.99, 1.56, 0.21, -0.31, -0.39, -0.04, -0.27, -1.57, -1.27,
]


def test_emanuel_jax_exit5_preserves_outcape_and_iflag():
    """Exits 1-4 pass but DTMA <= 0 keeps CBMF == 0 with CBMFOLD == 0 (exit
    5). Numba returns IFLAG == 1 and the CAPE it computed before exiting, here
    nonzero, so gating OUTCAPE or IFLAG on the final active flag shows up."""
    conv, np_state, timestep = emanuel_case(nlev=30, ncol=1, profile="unstable", as_jax=False)
    # Set on _params directly: EmanuelConvectionPython(DTMAX=...) raises
    # TypeError, because __init__ forwards its tuning keywords to sympl.
    conv._params = conv._params._replace(DTMAX=2.2)
    np_state["air_temperature"] = np_state["air_temperature"] + np.array(_EXIT5_TEMPERATURE_OFFSET)[:, None]
    q = np.zeros_like(np_state["specific_humidity"])
    q[:12, 0] = 0.0198 * np.exp(-np.arange(12) / 9.24)
    np_state["specific_humidity"] = q
    jax_state = {k: jnp.asarray(v) for k, v in np_state.items()}

    tend_np, diag_np = conv.array_call(np_state, timestep)
    tend_jax, diag_jax = conv.array_call(jax_state, timestep)

    assert diag_np["convective_state"][0] == 1
    assert diag_np["cloud_base_mass_flux"][0] == 0.0
    assert diag_np["atmosphere_convective_available_potential_energy"][0] != 0.0
    _assert_parity(tend_jax, diag_jax, tend_np, diag_np, "exit 5")

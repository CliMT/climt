import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax

import jax.numpy as jnp       # noqa: E402
import numpy as np            # noqa: E402

import climt                  # noqa: E402
from climt._core.jax_support import is_jax_array   # noqa: E402
from tests.test_jax_isolation import run            # noqa: E402


def test_is_jax_array():
    assert is_jax_array(jnp.zeros(3)) is True
    assert is_jax_array(np.zeros(3)) is False
    assert jax.jit(lambda x: jnp.asarray(float(is_jax_array(x))))(jnp.zeros(3)) == 1.0   # tracers count


def test_x64_off_is_refused_with_the_fix():
    code = (
        "import jax, jax.numpy as jnp, climt\n"
        "from datetime import timedelta\n"
        "from climt._core.jax_support import JaxPrecisionError\n"
        "assert not jax.config.jax_enable_x64\n"
        "state = {k: jnp.ones((30, 1)) for k in ('air_temperature', 'specific_humidity', 'air_pressure',"
        " 'eastward_wind', 'northward_wind')}\n"
        "state['air_pressure_on_interface_levels'] = jnp.ones((31, 1))\n"
        "try:\n"
        "    climt.EmanuelConvectionPython().array_call(state, timedelta(minutes=15))\n"
        "except JaxPrecisionError as error:\n"
        "    assert 'jax.config.update(\"jax_enable_x64\", True)' in str(error)\n"
        "else:\n"
        "    raise SystemExit('no error')\n"
    )
    result = run(code)
    assert result.returncode == 0, result.stderr


def test_cork_lw_x64_off_is_refused_with_the_fix():
    code = (
        "import jax, jax.numpy as jnp, climt\n"
        "from climt._core.jax_support import JaxPrecisionError\n"
        "assert not jax.config.jax_enable_x64\n"
        "comp = climt.CorkLongwaveRadiation(optics='correlated_k', table='earth_low_res_lw')\n"
        "try:\n"
        "    comp.array_call({'T': jnp.ones((20, 2))})\n"
        "except JaxPrecisionError as error:\n"
        "    assert 'jax.config.update(\"jax_enable_x64\", True)' in str(error)\n"
        "else:\n"
        "    raise SystemExit('no error')\n"
    )
    result = run(code)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("kwargs, word", [
    ({"optics": "correlated_k", "table": "earth_low_res_lw", "diagnostics_level": 1}, "diagnostics_level"),
    ({"optics": "parmentier"}, "optics"),
])
def test_cork_lw_refuses_uncovered_configurations(kwargs, word):
    jax.config.update("jax_enable_x64", True)
    component = climt.CorkLongwaveRadiation(**kwargs)
    with pytest.raises(NotImplementedError, match=word) as info:
        component.array_call({"T": jnp.ones((20, 2))})
    assert "CorkLongwaveRadiation" in str(info.value)


def test_cork_lw_refuses_non_7d_or_non_premixed_table():
    jax.config.update("jax_enable_x64", True)
    component = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    component._premixed_bg = False
    with pytest.raises(NotImplementedError, match="pre-mixed") as info:
        component.array_call({"T": jnp.ones((20, 2))})
    assert "CorkLongwaveRadiation" in str(info.value)

    component2 = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    component2._table = dict(component2._table, k_coefficients=np.zeros((1, 2, 2, 2, 2)))
    with pytest.raises(NotImplementedError, match="7-dimensional") as info2:
        component2.array_call({"T": jnp.ones((20, 2))})
    assert "CorkLongwaveRadiation" in str(info2.value)


def test_has_jax_kernel_flags():
    flagged = {cls.__name__ for cls in climt._components.__all__ if getattr(cls, "has_jax_kernel", False)}
    assert flagged == {"EmanuelConvectionPython", "CorkLongwaveRadiation"}


def test_dispatch_triggers_on_non_temperature_jax_input():
    jax.config.update("jax_enable_x64", True)
    from tests.jax_cases import emanuel_case
    from tests.test_cork_lw_jax_integration import _make_state

    # 1. Uncovered CORK config with NumPy T and JAX T_surf must still raise NotImplementedError
    parmentier = climt.CorkLongwaveRadiation(optics="parmentier")
    with pytest.raises(NotImplementedError, match="optics") as info:
        parmentier.array_call({"T": np.ones((10, 1)), "T_surf": jnp.ones((1,))})
    assert "CorkLongwaveRadiation" in str(info.value)

    # 2. Supported CORK config: jax.grad w.r.t. T_surf with NumPy T
    comp, s_np = _make_state(20, 2, as_jax=False)

    def olr_wrt_tsurf(T_surf):
        _, diag = comp.array_call({**s_np, "T_surf": T_surf})
        return diag["upwelling_longwave_flux_in_air"][-1].sum()

    g_ts = jax.grad(olr_wrt_tsurf)(jnp.asarray(s_np["T_surf"]))
    assert jnp.all(jnp.isfinite(g_ts)) and float(g_ts.sum()) > 0.0

    # 3. Emanuel: jax.grad w.r.t. specific_humidity with NumPy air_temperature
    conv, em_np, timestep = emanuel_case(nlev=30, ncol=1, profile="unstable", as_jax=False)

    def em_loss_wrt_q(q):
        tend, _ = conv.array_call({**em_np, "specific_humidity": q}, timestep)
        return jnp.sum(tend["air_temperature"] ** 2)

    g_q = jax.grad(em_loss_wrt_q)(jnp.asarray(em_np["specific_humidity"]))
    assert jnp.all(jnp.isfinite(g_q)) and jnp.any(g_q != 0.0)


def test_float32_arrays_with_x64_enabled_are_refused():
    from datetime import timedelta
    from climt._core.jax_support import JaxPrecisionError

    jax.config.update("jax_enable_x64", True)
    state_em = {
        k: jnp.ones((30, 1), dtype=jnp.float32)
        for k in ("air_temperature", "specific_humidity", "air_pressure", "eastward_wind", "northward_wind")
    }
    state_em["air_pressure_on_interface_levels"] = jnp.ones((31, 1), dtype=jnp.float32)
    with pytest.raises(JaxPrecisionError, match="float32") as exc_em:
        climt.EmanuelConvectionPython().array_call(state_em, timedelta(minutes=15))
    assert 'jax.config.update("jax_enable_x64", True)' in str(exc_em.value)

    comp_lw = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    with pytest.raises(JaxPrecisionError, match="float32"):
        comp_lw.array_call({"T": jnp.ones((20, 2), dtype=jnp.float32)})


@pytest.mark.parametrize("change, word", [
    ({"overlap_method": np.array("esft")}, "additive"),
    ({"k_coefficients": np.zeros((2, 2, 2, 2, 2, 2, 2))}, "single-gas"),
    ({"co2_vmr_grid": None}, "CO2"),
])
def test_cork_lw_refuses_7d_tables_the_jax_optics_do_not_cover(change, word):
    jax.config.update("jax_enable_x64", True)
    component = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    table = dict(component._table, **change)
    if change.get("co2_vmr_grid", ...) is None:
        del table["co2_vmr_grid"]
        component._has_co2_axis = False
    component._table = table
    with pytest.raises(NotImplementedError, match=word) as info:
        component.array_call({"T": jnp.ones((20, 2))})
    assert "CorkLongwaveRadiation" in str(info.value)

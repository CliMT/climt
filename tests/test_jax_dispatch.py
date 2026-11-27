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

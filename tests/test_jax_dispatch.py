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


@pytest.mark.parametrize("kwargs, word", [
    ({"optics": "correlated_k", "table": "earth_low_res_lw", "diagnostics_level": 1}, "diagnostics_level"),
    ({"optics": "parmentier"}, "optics"),
    ({"optics": "correlated_k", "table": "single_band_gray_lw"}, "7-dimensional"),
])
def test_cork_lw_refuses_uncovered_configurations(kwargs, word):
    jax.config.update("jax_enable_x64", True)
    component = climt.CorkLongwaveRadiation(**kwargs)
    with pytest.raises(NotImplementedError, match=word) as info:
        component.array_call({"T": jnp.ones((20, 2))})
    assert "CorkLongwaveRadiation" in str(info.value)


def test_cork_lw_x64_off_is_refused():
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


def test_cork_lw_lazy_table_inside_jit_does_not_leak_tracers():
    jax.config.update("jax_enable_x64", True)
    comp = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    assert comp._jax_table is None
    nband = comp.num_longwave_bands
    nlev, ncol = 10, 2
    state = {
        "T": jnp.full((nlev, ncol), 280.0),
        "p": jnp.linspace(95000.0, 5000.0, nlev)[:, None] * jnp.ones((1, ncol)),
        "p_int": jnp.linspace(100000.0, 0.0, nlev + 1)[:, None] * jnp.ones((1, ncol)),
        "T_surf": jnp.full((ncol,), 290.0),
        "emissivity": jnp.ones((nband, ncol)),
        "h2o": jnp.full((nlev, ncol), 0.005),
        "co2": jnp.full((nlev, ncol), 400e-6),
        "tau_cloud_lw": jnp.zeros((nlev, ncol, nband)),
    }
    tend_jit = jax.jit(lambda T: comp.array_call({**state, "T": T})[0]["T"])(state["T"])
    tend_eager, _ = comp.array_call(state)
    np.testing.assert_allclose(np.asarray(tend_jit), np.asarray(tend_eager["T"]), rtol=1e-12, atol=1e-14)


def test_has_jax_kernel_flags():
    flagged = {cls.__name__ for cls in climt._components.__all__ if getattr(cls, "has_jax_kernel", False)}
    assert flagged == {"EmanuelConvectionPython", "CorkLongwaveRadiation"}


def test_float32_arrays_with_x64_enabled_are_refused():
    from datetime import timedelta
    from climt._core.jax_support import JaxPrecisionError

    jax.config.update("jax_enable_x64", True)
    state_em = {
        k: jnp.ones((30, 1), dtype=jnp.float32)
        for k in ("air_temperature", "specific_humidity", "air_pressure", "eastward_wind", "northward_wind")
    }
    state_em["air_pressure_on_interface_levels"] = jnp.ones((31, 1), dtype=jnp.float32)
    with pytest.raises(JaxPrecisionError) as exc_em:
        climt.EmanuelConvectionPython().array_call(state_em, timedelta(minutes=15))
    assert 'jax.config.update("jax_enable_x64", True)' in str(exc_em.value)

    comp_lw = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    nband = comp_lw.num_longwave_bands
    state_lw = {
        "T": jnp.full((10, 2), 280.0, dtype=jnp.float32),
        "p": jnp.full((10, 2), 50000.0, dtype=jnp.float32),
        "p_int": jnp.linspace(100000.0, 0.0, 11, dtype=jnp.float32)[:, None] * jnp.ones((1, 2), dtype=jnp.float32),
        "T_surf": jnp.full((2,), 290.0, dtype=jnp.float32),
        "emissivity": jnp.ones((nband, 2), dtype=jnp.float32),
        "h2o": jnp.full((10, 2), 0.005, dtype=jnp.float32),
        "co2": jnp.full((10, 2), 400e-6, dtype=jnp.float32),
        "tau_cloud_lw": jnp.zeros((10, 2, nband), dtype=jnp.float32),
    }
    with pytest.raises(JaxPrecisionError) as exc_lw:
        comp_lw.array_call(state_lw)
    assert 'jax.config.update("jax_enable_x64", True)' in str(exc_lw.value)


def test_cork_lw_refuses_esft_or_multigas_7d_table():
    jax.config.update("jax_enable_x64", True)
    comp = climt.CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    comp._table = dict(comp._table, overlap_method=np.array("esft"))
    with pytest.raises(NotImplementedError, match="additive") as info:
        comp.array_call({"T": jnp.ones((20, 2), dtype=jnp.float64)})
    assert "CorkLongwaveRadiation" in str(info.value)

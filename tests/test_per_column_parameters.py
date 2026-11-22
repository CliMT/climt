"""Per-column parameters: with_parameters and the per-column numba kernels."""
import copy
import sys
import types
from datetime import timedelta

import numpy as np
import pytest

import climt
from climt import get_default_state, get_grid
from climt._components.emanuel import pure_python_v3
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPython
from climt._core.parameters import PARALLEL_MIN_COLUMNS, per_column, with_parameters

NCOL = 6
TIMESTEP = timedelta(minutes=15)


def _make_state(component, ncol, nlev=30):
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([component], grid_state=grid)
    t = state["air_temperature"].values
    q = state["specific_humidity"].values
    for j in range(ncol):
        for k in range(nlev):
            t[k, :, j] = 305.0 + 0.8 * j - 2.0 * k
        q[0:10, :, j] = 0.012 + 0.001 * j
    state["eastward_wind"].values[:] = 5.0 + np.arange(ncol)[np.newaxis, np.newaxis, :]
    return state


def _take_column(value, i):
    """Column i of a state quantity (sympl DataArray or UnytStateContainer)
    as a numpy array keeping the lon axis, or None if it has no lon axis."""
    dims = list(getattr(value, "dims", ()))
    if "lon" not in dims:
        return None
    return np.take(np.asarray(value.values), [i], axis=dims.index("lon"))


def _column(component, state, i):
    """A one-column state holding column i of state."""
    nlev = state["air_temperature"].shape[list(state["air_temperature"].dims).index("mid_levels")]
    out = _make_state(component, 1, nlev)
    for key, value in state.items():
        column = _take_column(value, i)
        if column is not None and key in out:
            out[key].values[...] = column
    return out


def _call(component, state):
    result = component(state, TIMESTEP)
    merged = {}
    for part in result:
        merged.update(part)
    return merged


def _component(name):
    return getattr(climt, name)() if name != "EmanuelConvectionPython" \
        else EmanuelConvectionPython()


@pytest.mark.parametrize("name,parameter,values", [
    ("EmanuelConvectionPython", "ALPHA", np.linspace(0.1, 0.4, 6)),
    ("SimplePhysics", "base_momentum_drag_coefficient", np.linspace(5e-4, 1.2e-3, 6)),
])
def test_per_column_equals_columns_alone(name, parameter, values):
    component = _component(name)
    state = _make_state(component, NCOL)
    together = _call(with_parameters(component, {parameter: values}, NCOL), state)
    varied = False
    for i in range(NCOL):
        alone_component = copy.copy(component)
        alone_component._params = component._params._replace(**{parameter: values[i]})
        alone = _call(alone_component, _column(component, state, i))
        for key, expected in alone.items():
            got = _take_column(together[key], i)
            if got is not None:
                expected = np.asarray(expected.values)
                np.testing.assert_allclose(
                    got, expected.reshape(got.shape),
                    rtol=1e-14, atol=0, err_msg="{} column {}".format(key, i))
                varied = varied or np.any(expected != 0)
    assert varied


def test_undeclared_parameter_is_refused():
    with pytest.raises(ValueError) as info:
        with_parameters(EmanuelConvectionPython(), {"IPBL": np.zeros(6)}, 6)
    message = str(info.value)
    for word in ("EmanuelConvectionPython", "IPBL", "ALPHA"):
        assert word in message


def test_original_component_untouched():
    component = EmanuelConvectionPython()
    component._params = component._params._replace(ALPHA=0.2)
    before = component._params
    copied = with_parameters(component, {"ALPHA": np.linspace(0.1, 0.4, 6)}, 6)
    assert component._params is before
    assert component._params.ALPHA == 0.2
    assert copied._params.ALPHA.shape == (6,)
    assert copied is not component


def test_simple_physics_copy_does_not_reuse_cache():
    component = climt.SimplePhysics()
    state = _make_state(component, NCOL)
    _call(component, state)
    values = np.linspace(5e-4, 1.2e-3, NCOL)
    copied = with_parameters(component, {"base_momentum_drag_coefficient": values}, NCOL)
    np.testing.assert_array_equal(
        copied._per_column_params(NCOL).base_momentum_drag_coefficient, values)
    assert np.all(component._per_column_params(NCOL).base_momentum_drag_coefficient
                  == component._params.base_momentum_drag_coefficient)


def test_emanuel_serial_and_parallel_agree():
    ncol = 64
    assert ncol >= PARALLEL_MIN_COLUMNS
    nlev = 30
    rng = np.random.default_rng(0)
    p = np.linspace(1000.0, 30.0, nlev)[:, None] * np.ones((1, ncol))
    ph = np.linspace(1005.0, 20.0, nlev + 1)[:, None] * np.ones((1, ncol))
    t = 305.0 - 2.0 * np.arange(nlev)[:, None] + rng.uniform(0, 3, (1, ncol))
    q = np.zeros((nlev, ncol))
    q[:10] = 0.012 + rng.uniform(0, 0.004, (1, ncol))
    u = np.ones((nlev, ncol))
    v = np.zeros((nlev, ncol))
    comp = EmanuelConvectionPython()
    qs = pure_python_v3.compute_qs(t, p, comp._cond, comp._params.RD)
    params = per_column(comp._params._replace(ALPHA=np.linspace(0.1, 0.4, ncol)), ncol)
    tra = np.zeros((nlev, 1, ncol))
    args = (t, q, qs, u, v, p, ph, nlev, nlev - 3, 0, 900.0, np.zeros(ncol), tra,
            params, comp._cond)
    serial = pure_python_v3._convect_serial(*args)
    parallel = pure_python_v3._convect_parallel(*args)
    for a, b in zip(serial, parallel):
        np.testing.assert_allclose(a, b, rtol=0, atol=0)


@pytest.mark.parametrize("ncol, expected", [
    (PARALLEL_MIN_COLUMNS - 1, "serial"), (PARALLEL_MIN_COLUMNS, "parallel")])
def test_emanuel_picks_variant_by_column_count(monkeypatch, ncol, expected):
    called = []
    for variant in ("serial", "parallel"):
        monkeypatch.setattr(pure_python_v3, "_convect_" + variant,
                            lambda *args, _v=variant: called.append(_v))
    t = np.ones((3, ncol))
    pure_python_v3._numpy_vectorized_convect(
        t, t, t, t, t, t, np.ones((4, ncol)), 3, 0, 0, 1.0, np.ones(ncol),
        np.zeros((3, 1, ncol)), None, None)
    assert called == [expected]


def test_int_fields_stay_int64_under_postponed_annotations():
    module = types.ModuleType("_postponed_params_mod")
    sys.modules[module.__name__] = module
    try:
        exec(
            "from __future__ import annotations\n"
            "from typing import NamedTuple\n"
            "class P(NamedTuple):\n"
            "    n: int\n"
            "    x: float\n",
            module.__dict__)
        result = per_column(module.P(n=3, x=1.5), 4)
    finally:
        del sys.modules[module.__name__]
    assert result.n.dtype == np.int64
    assert result.x.dtype == np.float64

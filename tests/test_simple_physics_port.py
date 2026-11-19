"""Parity of the numba SimplePhysics kernel with the Fortran reference.

The reference outputs in ``tests/fixtures/fortran_reference`` were produced
by the Fortran-backed ``SimplePhysics.array_call`` (see
``scripts/generate_fortran_reference.py``). Here the kernel is called on the
same inputs, with parameters built the way the component builds them.
"""

import inspect
import json
import os
import subprocess
import sys
from datetime import timedelta

import numpy as np
import pytest
from sympl import get_constant

import climt
from climt import SimplePhysics
from climt._components.simple_physics import kernel as sp_kernel
from climt._components.simple_physics.kernel import (
    SimplePhysicsFlags, SimplePhysicsParams, simple_physics_kernel)
from climt._core.parameters import PARALLEL_MIN_COLUMNS, per_column

FIXTURES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "fixtures", "fortran_reference")

with open(os.path.join(FIXTURES, "cases.json")) as _f:
    CASES = json.load(_f)["simple_physics"]
SIMPLE_PHYSICS_CASES = sorted(CASES)

# Kernel output order -> fixture quantity name.
OUTPUT_NAMES = (
    "air_temperature",
    "eastward_wind",
    "northward_wind",
    "specific_humidity",
    "stratiform_precipitation_rate",
    "surface_upward_sensible_heat_flux",
    "surface_upward_latent_heat_flux",
)

# Physical constants as SimplePhysics reads them: params field -> units.
CONSTANT_UNITS = {
    "gravitational_acceleration": "m/s^2",
    "heat_capacity_of_dry_air_at_constant_pressure": "J/kg/degK",
    "gas_constant_of_dry_air": "J/kg/degK",
    "gas_constant_of_vapor_phase": "J/kg/degK",
    "latent_heat_of_condensation": "J/kg",
    "planetary_radius": "m",
    "planetary_rotation_rate": "s^-1",
    "density_of_liquid_water": "kg/m^3",
}


def _load_npz():
    return np.load(os.path.join(FIXTURES, "simple_physics.npz"))


def constructor_kwargs(case):
    """The case's kwargs on top of SimplePhysics's constructor defaults."""
    sig = inspect.signature(climt.SimplePhysics.__init__)
    kwargs = {name: p.default for name, p in sig.parameters.items()
              if p.default is not inspect.Parameter.empty}
    kwargs.update(CASES[case]["kwargs"])
    return kwargs


def build_params_and_flags(case, **overrides):
    kwargs = constructor_kwargs(case)
    fields = {}
    for name in SimplePhysicsParams._fields:
        if name in CONSTANT_UNITS:
            fields[name] = get_constant(name, CONSTANT_UNITS[name])
        else:
            fields[name] = kwargs[name]
    fields.update(overrides)
    flags = SimplePhysicsFlags(
        **{name: bool(kwargs[name]) for name in SimplePhysicsFlags._fields})
    return SimplePhysicsParams(**fields), flags


def fixture_inputs(case):
    """Kernel positional inputs, in the kernel's (nlev, ncol)/(ncol,) layout."""
    data = _load_npz()

    def get(quantity):
        return np.array(data["{}/in/{}".format(case, quantity)], dtype=np.float64)

    t = get("air_temperature")
    nlev, ncol = t.shape
    arrays = (
        get("eastward_wind"), get("northward_wind"), t,
        get("specific_humidity"), get("air_pressure"),
        get("air_pressure_on_interface_levels"),
        get("surface_air_pressure"), get("surface_temperature"),
        get("surface_specific_humidity"), get("latitude"),
    )
    for a in arrays[:5]:
        assert a.shape == (nlev, ncol)
    assert arrays[5].shape == (nlev + 1, ncol)
    for a in arrays[6:]:
        assert a.shape == (ncol,)
    return arrays


def fixture_outputs(case):
    data = _load_npz()
    return {name: data["{}/out/{}".format(case, name)] for name in OUTPUT_NAMES}


def call_kernel(func, arrays, case, params, flags):
    ncol = arrays[2].shape[1]
    before = [a.copy() for a in arrays]
    out = func(*arrays, CASES[case]["timestep_seconds"],
               per_column(params, ncol), flags)
    for a, b in zip(arrays, before):
        np.testing.assert_array_equal(a, b)  # inputs unmodified
    out = dict(zip(OUTPUT_NAMES, out))
    # SimplePhysics clips negative latent heat flux after the kernel.
    lhf = out["surface_upward_latent_heat_flux"].copy()
    lhf[lhf < 0] = 0
    out["surface_upward_latent_heat_flux"] = lhf
    return out


def run_kernel(case):
    params, flags = build_params_and_flags(case)
    return call_kernel(simple_physics_kernel, fixture_inputs(case), case,
                       params, flags)


@pytest.mark.parametrize("case", SIMPLE_PHYSICS_CASES)
def test_kernel_matches_fortran(case):
    out = run_kernel(case)
    for name, expected in fixture_outputs(case).items():
        np.testing.assert_allclose(out[name], expected, rtol=1e-10,
                                   atol=1e-12, err_msg=name)


def test_all_eleven_cases_present():
    assert len(SIMPLE_PHYSICS_CASES) == 11


def test_per_column_parameters_match_columns_run_alone():
    case = "default"
    arrays = fixture_inputs(case)
    ncol = arrays[2].shape[1]
    cd0 = np.linspace(0.0005, 0.0012, ncol)
    params, flags = build_params_and_flags(
        case, base_momentum_drag_coefficient=cd0)
    together = call_kernel(simple_physics_kernel, arrays, case, params, flags)
    for i in range(ncol):
        column = tuple(
            np.ascontiguousarray(a[..., i:i + 1]) for a in arrays)
        params_i = params._replace(base_momentum_drag_coefficient=cd0[i])
        alone = call_kernel(simple_physics_kernel, column, case, params_i,
                            flags)
        for name in OUTPUT_NAMES:
            np.testing.assert_allclose(
                together[name][..., i:i + 1], alone[name], rtol=1e-14,
                atol=0, err_msg="{} column {}".format(name, i))


def test_parallel_and_serial_variants_agree():
    case = "default"
    reps = 8
    arrays = tuple(np.tile(a, reps) for a in fixture_inputs(case))
    ncol = arrays[2].shape[1]
    assert ncol == 64 and ncol >= PARALLEL_MIN_COLUMNS
    params, flags = build_params_and_flags(case)
    dispatched = call_kernel(simple_physics_kernel, arrays, case, params,
                             flags)
    serial = call_kernel(sp_kernel._simple_physics_serial, arrays, case,
                         params, flags)
    for name in OUTPUT_NAMES:
        np.testing.assert_allclose(dispatched[name], serial[name], rtol=0,
                                   atol=0, err_msg=name)


@pytest.mark.parametrize("ncol, expected", [
    (PARALLEL_MIN_COLUMNS - 1, "serial"), (PARALLEL_MIN_COLUMNS, "parallel")])
def test_kernel_picks_variant_by_column_count(monkeypatch, ncol, expected):
    called = []
    for variant in ("serial", "parallel"):
        monkeypatch.setattr(
            sp_kernel, "_simple_physics_" + variant,
            lambda *args, _v=variant: called.append(_v))
    nlev = 3
    col2, col1 = np.ones((nlev, ncol)), np.ones(ncol)
    params, flags = build_params_and_flags("default")
    simple_physics_kernel(col2, col2, col2, col2, col2,
                          np.ones((nlev + 1, ncol)), col1, col1, col1, col1,
                          600.0, per_column(params, ncol), flags)
    assert called == [expected]


# --- The SimplePhysics component on the kernel ---------------------------

def component_inputs(case):
    """The raw dict the reference passed to ``SimplePhysics.array_call``."""
    data = _load_npz()
    prefix = "{}/in/".format(case)
    return {key[len(prefix):]: np.array(data[key]) for key in data.files
            if key.startswith(prefix)}


def component_timestep(case):
    return timedelta(seconds=CASES[case]["timestep_seconds"])


def run_component(component, raw, case="default"):
    diagnostics, new_state = component.array_call(raw, component_timestep(case))
    out = dict(new_state)
    out.update(diagnostics)
    return out


@pytest.mark.parametrize("case", SIMPLE_PHYSICS_CASES)
def test_component_matches_fortran(case):
    component = SimplePhysics(**CASES[case]["kwargs"])
    out = run_component(component, component_inputs(case), case)
    for name, expected in fixture_outputs(case).items():
        np.testing.assert_allclose(out[name], expected, rtol=1e-10,
                                   atol=1e-12, err_msg=name)


def test_alternating_instances_keep_their_parameters():
    a = SimplePhysics()
    b = SimplePhysics(base_momentum_drag_coefficient=0.0014)
    raw = component_inputs("default")
    first = run_component(a, raw)
    other = run_component(b, raw)
    again = run_component(a, raw)
    assert any(not np.array_equal(first[k], other[k]) for k in first)
    for name in first:
        np.testing.assert_allclose(again[name], first[name], rtol=0, atol=0,
                                   err_msg=name)


def test_non_contiguous_input_matches():
    raw = component_inputs("default")
    fortran_ordered = {
        k: (np.asfortranarray(v) if v.ndim == 2 else v)
        for k, v in raw.items()}
    assert not fortran_ordered["air_temperature"].flags["C_CONTIGUOUS"]
    expected = run_component(SimplePhysics(), raw)
    actual = run_component(SimplePhysics(), fortran_ordered)
    for name in expected:
        np.testing.assert_allclose(actual[name], expected[name], rtol=0,
                                   atol=0, err_msg=name)


def test_runs_without_climt_fortran():
    code = (
        "import sys\n"
        "sys.modules['climt_fortran'] = None\n"
        "from datetime import timedelta\n"
        "import climt\n"
        "component = climt.SimplePhysics()\n"
        "state = climt.get_default_state([component])\n"
        "component(state, timedelta(minutes=10))\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().endswith("ok")


_NO_JIT_RUN = """
import numpy as np
import sympl
import climt

sympl.set_backend(climt.UnytBackend())
component = climt.SimplePhysics()
state = climt.get_default_state(
    [component], grid_state=climt.get_grid(nx=1, ny=1, nz=28))
diagnostics, new_state = component(state, climt.UnytTimeDelta(minutes=10))
print(float(np.asarray(new_state['air_temperature']).ravel()[0]))
"""


def test_component_runs_with_the_jit_disabled():
    """Runs the Pyodide-shaped path: no numba JIT, UnytBackend, one step.

    numba reads NUMBA_DISABLE_JIT at import time, so this has to happen in a
    fresh interpreter.
    """
    result = subprocess.run(
        [sys.executable, "-c", _NO_JIT_RUN],
        env={**os.environ, "NUMBA_DISABLE_JIT": "1"},
        capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, (
        "SimplePhysics failed with the JIT disabled -- it will fail the same "
        "way in Pyodide, which has no numba at all:\n" + result.stderr)
    assert np.isfinite(float(result.stdout.strip().splitlines()[-1]))


def test_constants_read_once(monkeypatch):
    component = SimplePhysics()
    monkeypatch.setattr(
        "climt._components.simple_physics.component.get_constant",
        lambda *a, **k: pytest.fail("get_constant called in array_call"))
    run_component(component, component_inputs("default"))


def test_per_column_parameters_built_once_per_column_count(monkeypatch):
    from climt._components.simple_physics import component as sp_component
    calls = []

    def counting_per_column(params, ncol):
        calls.append(ncol)
        return per_column(params, ncol)

    monkeypatch.setattr(sp_component, "per_column", counting_per_column)
    component = SimplePhysics()
    raw = component_inputs("default")
    ncol = raw["air_temperature"].shape[1]
    first = run_component(component, raw)
    again = run_component(component, raw)
    assert calls == [ncol]
    for name in first:
        np.testing.assert_array_equal(again[name], first[name], err_msg=name)

    # A different column count rebuilds them, and so does going back.
    one = {k: np.ascontiguousarray(v[..., :1]) for k, v in raw.items()}
    single = run_component(component, one)
    assert calls == [ncol, 1]
    for name in first:
        np.testing.assert_allclose(single[name], first[name][..., :1],
                                   rtol=1e-14, atol=0, err_msg=name)
    run_component(component, raw)
    assert calls == [ncol, 1, ncol]


# --- Internal surface temperature ------------------------------------------

def _run_with_latitude(case, latitude):
    component = SimplePhysics(**CASES[case]["kwargs"])
    raw = component_inputs(case)
    raw["latitude"] = np.asarray(latitude, dtype=np.float64)
    return run_component(component, raw, case)


def test_cyclone_sst_does_not_depend_on_latitude():
    # Test 5-1 uses a constant SST of 302.15 K.
    raw = component_inputs("cyclone")
    ncol = raw["latitude"].shape[0]
    a = _run_with_latitude("cyclone", np.zeros(ncol))
    b = _run_with_latitude("cyclone", np.full(ncol, 45.0))
    for name in a:
        np.testing.assert_array_equal(a[name], b[name], err_msg=name)


def test_internal_sst_depends_on_latitude():
    # Without the cyclone flag, the internal SST is the latitude-dependent
    # SST of the moist baroclinic wave (test 4-3), symmetric about the
    # equator.
    raw = component_inputs("internal_sst")
    ncol = raw["latitude"].shape[0]
    equator = _run_with_latitude("internal_sst", np.zeros(ncol))
    north = _run_with_latitude("internal_sst", np.full(ncol, 60.0))
    south = _run_with_latitude("internal_sst", np.full(ncol, -60.0))
    flux = "surface_upward_sensible_heat_flux"
    assert np.all(equator[flux] > north[flux])
    for name in north:
        np.testing.assert_allclose(north[name], south[name], rtol=1e-14,
                                   atol=0, err_msg=name)

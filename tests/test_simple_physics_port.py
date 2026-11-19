"""Parity of the numba SimplePhysics kernel with the Fortran reference.

The reference outputs in ``tests/fixtures/fortran_reference`` were produced
by the Fortran-backed ``SimplePhysics.array_call`` (see
``scripts/generate_fortran_reference.py``). Here the kernel is called on the
same inputs, with parameters built the way the component builds them.
"""

import inspect
import json
import os

import numpy as np
import pytest
from sympl import get_constant

import climt
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

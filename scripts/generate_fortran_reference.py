"""Generate Fortran reference outputs for SimplePhysics and DCMIP ICs.

These fixtures are the parity oracle for the numba ports of both components.
They are ONLY valid when generated from the Fortran-backed components, i.e.
before the numba port replaces them. Once SimplePhysics and
DcmipInitialConditions run numba, re-running this script would compare the
port with itself and would NOT reproduce a Fortran reference; do not
regenerate the fixtures after that point.

The script refuses to run unless both components are still Fortran-backed:
the pre-port component modules hold the loaded extension (``phys`` in
``simple_physics/component.py``, ``_dcmip`` in ``dcmip/component.py``). From
the numba port on (after commit b7c3302e) it always refuses; to add a case,
check out b7c3302e or earlier, with climt-fortran installed, and run it there.
It also refuses unless ``fortran_status()['state'] == 'available'``.

Writes, under ``tests/fixtures/fortran_reference/``:

* ``cases.json``: constructor kwargs (and timestep) for every case.
* ``simple_physics.npz`` and ``dcmip.npz``: keys ``"<case>/in/<quantity>"``
  (raw input arrays passed to ``array_call``, in the component's dims) and
  ``"<case>/out/<quantity>"`` (every diagnostic and new-state quantity).

Run from the worktree root:
  python scripts/generate_fortran_reference.py
"""

import json
import os
import sys
from datetime import timedelta

import numpy as np
from sympl import get_numpy_arrays_with_properties

import climt
from climt._core.fortran import fortran_status

OUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir,
    "tests", "fixtures", "fortran_reference")
SEED = 20261119
TIMESTEP_SECONDS = 600.0

# Saturation specific humidity as in the Fortran large-scale condensation.
EPSILON = 0.622
E0 = 610.78
T0 = 273.16
LATVAP = 2.501e6
RH2O = 461.5

SIMPLE_PHYSICS_CASES = {
    "default": {},
    "lsc_off": {"large_scale_condensation": False},
    "pbl_off": {"boundary_layer": False},
    "fluxes_off": {"surface_fluxes": False},
    "internal_sst": {"use_external_surface_temperature": False},
    "cyclone": {"simulate_cyclone": True,
                "use_external_surface_temperature": False},
    "external_qsurf": {"use_external_surface_specific_humidity": True},
    "drag_params": {
        "top_of_boundary_layer": 1.5 * 85000.0,
        "boundary_layer_influence_height": 1.5 * 20000.0,
        "drag_coefficient_heat_fluxes": 1.5 * 0.0011,
        "base_momentum_drag_coefficient": 1.5 * 0.0007,
        "wind_dependent_momentum_drag_coefficient": 1.5 * 0.000065,
        "maximum_momentum_drag_coefficient": 1.5 * 0.002,
    },
    "lsc_supersaturated": {},
    "strong_wind": {},
    "cold_surface": {},
}

DCMIP_CASES = {
    "bw_dry": {"moist": False},
    "bw_moist": {"moist": True},
    "bw_unperturbed": {"moist": False, "add_perturbation": False},
    "tc": {"condition_type": "tropical_cyclone", "moist": True},
}


def qsat(t, p):
    return EPSILON * E0 / p * np.exp(-LATVAP / RH2O * (1.0 / t - 1.0 / T0))


def raw_inputs(component, grid):
    state = climt.get_default_state([component], grid_state=grid)
    raw = get_numpy_arrays_with_properties(state, component.input_properties)
    return {k: np.array(v, dtype=np.float64, copy=True)
            for k, v in raw.items()}


def simple_physics_state(base, case, rng):
    """Give every column a distinct state; then apply case-specific changes.

    Arrays have shape (levels, columns); level 0 is the surface (highest
    pressure).
    """
    s = {k: v.copy() for k, v in base.items()}
    p = s["air_pressure"]
    nz, ncol = p.shape
    ps = s["surface_air_pressure"]
    # Surface air temperature 220-305 K across columns, dry-adiabatic-ish
    # lapse with a 190 K floor plus small noise.
    t_sfc = np.linspace(220.0, 305.0, ncol) + rng.uniform(-1, 1, ncol)
    t = t_sfc[None, :] * (p / ps[None, :]) ** (287.0 * 0.0065 / 9.81)
    t = np.maximum(t, 190.0) + rng.uniform(-0.5, 0.5, t.shape)
    s["air_temperature"] = t
    # Relative humidity 0-1.2 of saturation, distinct per column.
    rh = rng.uniform(0.0, 1.2, (nz, ncol))
    s["specific_humidity"] = rh * qsat(t, p)
    s["eastward_wind"] = rng.uniform(-10, 10, (nz, ncol))
    s["northward_wind"] = rng.uniform(-10, 10, (nz, ncol))
    s["surface_temperature"] = t_sfc + rng.uniform(-2, 2, ncol)
    s["surface_specific_humidity"] = rng.uniform(0.5, 1.0, ncol) * qsat(
        s["surface_temperature"], ps)
    s["latitude"] = np.linspace(-60, 60, ncol) + rng.uniform(-2, 2, ncol)

    if case == "lsc_supersaturated":
        factor = rng.uniform(1.1, 1.3, (11, ncol))
        s["specific_humidity"][5:16] = factor * qsat(t[5:16], p[5:16])
    elif case == "strong_wind":
        s["eastward_wind"][0] = np.linspace(15.0, 30.0, ncol)
    elif case == "cold_surface":
        s["surface_temperature"] = np.linspace(250.0, 275.0, ncol)
    return s


def generate_simple_physics(rng):
    grid = climt.get_grid(nx=4, ny=2, nz=28)
    base = raw_inputs(climt.SimplePhysics(), grid)
    assert base["air_pressure"][0, 0] > base["air_pressure"][-1, 0]
    arrays, cases = {}, {}
    for case, kwargs in SIMPLE_PHYSICS_CASES.items():
        component = climt.SimplePhysics(**kwargs)
        state = simple_physics_state(base, case, rng)
        inputs = {k: v.copy() for k, v in state.items()}
        diagnostics, new_state = component.array_call(
            state, timedelta(seconds=TIMESTEP_SECONDS))
        for k, v in inputs.items():
            arrays["{}/in/{}".format(case, k)] = v
        for k, v in list(diagnostics.items()) + list(new_state.items()):
            arrays["{}/out/{}".format(case, k)] = np.asarray(v)
        cases[case] = {"kwargs": kwargs, "timestep_seconds": TIMESTEP_SECONDS}
    return arrays, cases


def generate_dcmip():
    grid = climt.get_grid(nx=16, ny=8, nz=10)
    arrays, cases = {}, {}
    for case, kwargs in DCMIP_CASES.items():
        component = climt.DcmipInitialConditions(**kwargs)
        state = raw_inputs(component, grid)
        inputs = {k: v.copy() for k, v in state.items()}
        diagnostics = component.array_call(state)
        for k, v in inputs.items():
            arrays["{}/in/{}".format(case, k)] = v
        for k, v in diagnostics.items():
            arrays["{}/out/{}".format(case, k)] = np.asarray(v)
        cases[case] = {"kwargs": kwargs}
    return arrays, cases


def check_fortran_backed():
    """Exit unless SimplePhysics and DcmipInitialConditions call the Fortran.

    ``fortran_status()`` alone no longer tells: after the numba port it
    reports ``available`` without ``_simple_physics`` or ``_dcmip``, and the
    components then run numba, so the "reference" would be the port itself.
    """
    from climt._components.dcmip import component as dcmip_component
    from climt._components.simple_physics import component as sp_component
    numba_backed = [
        name for name, module, attribute in (
            ("SimplePhysics", sp_component, "phys"),
            ("DcmipInitialConditions", dcmip_component, "_dcmip"))
        if getattr(module, attribute, None) is None]
    if numba_backed:
        sys.exit(
            "{} not Fortran-backed in this checkout (numba port), so the "
            "output would not be a Fortran reference. Regenerate the "
            "fixtures from a pre-port commit (b7c3302e or earlier) with "
            "climt-fortran installed.".format(
                " and ".join(numba_backed)
                + (" are" if len(numba_backed) > 1 else " is")))


def main():
    check_fortran_backed()
    state = fortran_status()["state"]
    if state != "available":
        sys.exit("Fortran extensions not available (state={!r}); the "
                 "reference must come from the Fortran.".format(state))
    rng = np.random.default_rng(SEED)
    sp_arrays, sp_cases = generate_simple_physics(rng)
    dc_arrays, dc_cases = generate_dcmip()
    for name, arrays in (("simple_physics", sp_arrays), ("dcmip", dc_arrays)):
        for k, v in arrays.items():
            if not np.all(np.isfinite(v)):
                sys.exit("non-finite values in {} {}".format(name, k))
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez(os.path.join(OUT_DIR, "simple_physics.npz"), **sp_arrays)
    np.savez(os.path.join(OUT_DIR, "dcmip.npz"), **dc_arrays)
    with open(os.path.join(OUT_DIR, "cases.json"), "w") as f:
        json.dump({"simple_physics": sp_cases, "dcmip": dc_cases}, f,
                  indent=2, sort_keys=True)
        f.write("\n")


if __name__ == "__main__":
    main()

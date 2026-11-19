"""Time the numba SimplePhysics kernel against the Fortran, one column.

One column of fixture case ``default``. The Fortran is called through
``climt_fortran._simple_physics`` (``init_simple_physics`` and
``set_physical_constants`` once, then ``get_new_state`` timed); the numba
kernel is timed after one warm-up call, with its parameters built once.
2,000 calls each, best of 5. Run with one numba thread:

  NUMBA_NUM_THREADS=1 python benchmarks/benchmark_simple_physics_port.py

Prints the timings and writes them, with the machine, numba version and
commit, to ``benchmarks/results/simple_physics_port.txt``.
"""

import json
import os
import platform
import subprocess
import time

import numba
import numpy as np
from climt_fortran import _simple_physics as fortran
from sympl import get_constant

from climt._components.simple_physics.kernel import (
    SimplePhysicsFlags, SimplePhysicsParams, simple_physics_kernel)
from climt._core.parameters import per_column

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
FIXTURES = os.path.join(ROOT, "tests", "fixtures", "fortran_reference")
RESULT = os.path.join(ROOT, "benchmarks", "results", "simple_physics_port.txt")
CASE = "default"
CALLS = 2000
REPEATS = 5

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
# SimplePhysics constructor defaults (case "default" has no kwargs).
DEFAULTS = {
    "top_of_boundary_layer": 85000.0,
    "boundary_layer_influence_height": 20000.0,
    "drag_coefficient_heat_fluxes": 0.0011,
    "base_momentum_drag_coefficient": 0.0007,
    "wind_dependent_momentum_drag_coefficient": 0.000065,
    "maximum_momentum_drag_coefficient": 0.002,
}
FLAGS = SimplePhysicsFlags(
    simulate_cyclone=False, large_scale_condensation=True,
    boundary_layer=True, surface_fluxes=True,
    use_external_surface_temperature=True,
    use_external_surface_specific_humidity=False)


def one_column_inputs():
    data = np.load(os.path.join(FIXTURES, "simple_physics.npz"))

    def get(quantity):
        return np.ascontiguousarray(
            data["{}/in/{}".format(CASE, quantity)][..., :1], dtype=np.float64)

    return {q: get(q) for q in (
        "eastward_wind", "northward_wind", "air_temperature",
        "specific_humidity", "air_pressure",
        "air_pressure_on_interface_levels", "surface_air_pressure",
        "surface_temperature", "surface_specific_humidity", "latitude")}


def best_of(func):
    func()
    best = float("inf")
    for _ in range(REPEATS):
        start = time.perf_counter()
        for _ in range(CALLS):
            func()
        best = min(best, (time.perf_counter() - start) / CALLS)
    return best * 1e6


def main():
    with open(os.path.join(FIXTURES, "cases.json")) as f:
        dt = json.load(f)["simple_physics"][CASE]["timestep_seconds"]
    x = one_column_inputs()
    fields = {name: get_constant(name, units)
              for name, units in CONSTANT_UNITS.items()}
    fields.update(DEFAULTS)
    params = SimplePhysicsParams(**fields)

    fortran.init_simple_physics(*(int(flag) for flag in FLAGS))
    fortran.set_physical_constants(
        params.gravitational_acceleration,
        params.heat_capacity_of_dry_air_at_constant_pressure,
        params.gas_constant_of_dry_air,
        params.latent_heat_of_condensation,
        params.gas_constant_of_vapor_phase,
        params.planetary_radius,
        params.planetary_rotation_rate,
        params.density_of_liquid_water,
        params.top_of_boundary_layer,
        params.boundary_layer_influence_height,
        params.drag_coefficient_heat_fluxes,
        params.base_momentum_drag_coefficient,
        params.wind_dependent_momentum_drag_coefficient,
        params.maximum_momentum_drag_coefficient)

    def call_fortran():
        return fortran.get_new_state(
            x["eastward_wind"], x["northward_wind"], x["air_temperature"],
            x["air_pressure"], x["air_pressure_on_interface_levels"],
            x["specific_humidity"], x["surface_air_pressure"],
            x["surface_temperature"], x["surface_specific_humidity"],
            x["latitude"], dt)

    column_params = per_column(params, 1)

    def call_numba():
        return simple_physics_kernel(
            x["eastward_wind"], x["northward_wind"], x["air_temperature"],
            x["specific_humidity"], x["air_pressure"],
            x["air_pressure_on_interface_levels"], x["surface_air_pressure"],
            x["surface_temperature"], x["surface_specific_humidity"],
            x["latitude"], dt, column_params, FLAGS)

    for a, b in zip(call_fortran(), call_numba()):
        np.testing.assert_allclose(b, a, rtol=1e-10, atol=1e-12)

    t_fortran = best_of(call_fortran)
    t_numba = best_of(call_numba)
    nlev = x["air_temperature"].shape[0]
    lines = [
        "SimplePhysics, 1 column x {} levels, {} calls, best of {}".format(
            nlev, CALLS, REPEATS),
        "Fortran (get_new_state): {:8.2f} us/call".format(t_fortran),
        "numba (kernel):          {:8.2f} us/call".format(t_numba),
        "ratio numba / Fortran:   {:8.2f} (bound 1.5)".format(
            t_numba / t_fortran),
    ]
    print("\n".join(lines))
    commit = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
        capture_output=True, text=True).stdout.strip()
    meta = [
        "machine: {} ({})".format(platform.platform(), platform.machine()),
        "numba: {}, NUMBA_NUM_THREADS={}".format(
            numba.__version__, numba.config.NUMBA_NUM_THREADS),
        "commit: {}".format(commit),
    ]
    with open(RESULT, "w") as f:
        f.write("\n".join(lines + meta) + "\n")


if __name__ == "__main__":
    main()

import time
from datetime import timedelta

import numpy as np
import sympl
import unyt
from sympl._core.backend import DataArrayBackend

import climt
from climt import (
    BergerSolarInsolation,
    BucketHydrology,
    DryConvectiveAdjustment,
    EmanuelConvection,
    Frierson06LongwaveOpticalDepth,
    GrayLongwaveRadiation,
    GridScaleCondensation,
    HeldSuarez,
    IceSheet,
    Instellation,
    RRTMGLongwave,
    RRTMGShortwave,
    SimplePhysics,
    SlabSurface,
    UnytBackend,
    UnytStateContainer,
    get_default_state,
    get_grid,
)

# List of components to benchmark
COMPONENTS = [
    BergerSolarInsolation,
    BucketHydrology,
    DryConvectiveAdjustment,
    EmanuelConvection,
    Frierson06LongwaveOpticalDepth,
    GrayLongwaveRadiation,
    GridScaleCondensation,
    HeldSuarez,
    IceSheet,
    Instellation,
    RRTMGLongwave,
    RRTMGShortwave,
    SimplePhysics,
    SlabSurface,
]


def get_timestep_quantity(timestep, backend_obj):
    if isinstance(backend_obj, UnytBackend):
        from climt import UnytTimeDelta

        return UnytTimeDelta(timestep)
    elif isinstance(backend_obj, DataArrayBackend):
        from sympl import DataArray

        return DataArray(timestep.total_seconds(), attrs={"units": "s"})
    return timestep


def run_simulation(component, state, timestep, backend_obj, iterations):
    is_diagnostic = isinstance(component, sympl.DiagnosticComponent)
    is_implicit = isinstance(component, sympl.ImplicitTendencyComponent)
    is_stepper = isinstance(component, sympl.Stepper)

    dt_quant = get_timestep_quantity(timestep, backend_obj)

    start_time = time.perf_counter()

    for i in range(iterations):
        if is_diagnostic:
            diagnostics = component(state)
            state.update(diagnostics)
        elif is_stepper:
            diagnostics, new_state = component(state, timestep)
            state.update(diagnostics)
            state.update(new_state)
        elif is_implicit:
            tendencies, diagnostics = component(state, timestep)
            for name, tendency in tendencies.items():
                if name in state:
                    state[name] = state[name] + tendency * dt_quant
            state.update(diagnostics)
        else:
            # Explicit TendencyComponent
            tendencies, diagnostics = component(state)
            for name, tendency in tendencies.items():
                if name in state:
                    state[name] = state[name] + tendency * dt_quant
            state.update(diagnostics)

        state["time"] += timestep

    end_time = time.perf_counter()
    return end_time - start_time


def benchmark_component(component_cls, iterations=100):

    print(f"Benchmarking {component_cls.__name__}...", end="", flush=True)

    # Setup

    timestep = timedelta(minutes=20)

    try:
        # Initialize component (some might need specific args, assume defaults work for now)

        # Note: Radiation might need mcica=False to be stable/faster/deterministic?

        # But let's use defaults.

        component = component_cls()

    except Exception as e:
        print(f" FAILED initialization: {e}")

        return None, None

    # Benchmark DataArrayBackend

    sympl.set_backend(DataArrayBackend())

    grid_da = get_grid(nx=32, ny=32, nz=20)

    state_da = get_default_state([component], grid_state=grid_da)

    try:
        time_da = run_simulation(
            component, state_da, timestep, DataArrayBackend(), iterations
        )

    except Exception as e:
        print(f" FAILED DataArray: {e}")

        time_da = None

    # Benchmark UnytBackend

    sympl.set_backend(UnytBackend())

    grid_unyt = get_grid(nx=32, ny=32, nz=20)

    # Re-get state to ensure it uses the new backend

    state_unyt = get_default_state([component], grid_state=grid_unyt)

    try:
        time_unyt = run_simulation(
            component, state_unyt, timestep, UnytBackend(), iterations
        )

    except Exception as e:
        print(f" FAILED Unyt: {e}")

        time_unyt = None

    print(" Done.")
    return time_da, time_unyt


def main():
    iterations = 50
    results = []

    print(f"Starting benchmark suite with {iterations} iterations per component.")
    print(
        f"{'Component':<35} | {'DataArray (s)':<15} | {'Unyt (s)':<15} | {'Speedup':<10}"
    )
    print("-" * 85)

    for comp_cls in COMPONENTS:
        # Special handling for slow components if needed
        current_iters = iterations
        if comp_cls in [RRTMGLongwave, RRTMGShortwave]:
            current_iters = 10  # Radiation is slow

        time_da, time_unyt = benchmark_component(comp_cls, iterations=current_iters)

        if time_da is not None and time_unyt is not None:
            speedup = time_da / time_unyt
            results.append((comp_cls.__name__, time_da, time_unyt, speedup))
            print(
                f"{comp_cls.__name__:<35} | {time_da:<15.4f} | {time_unyt:<15.4f} | {speedup:<10.2f}x"
            )
        else:
            results.append((comp_cls.__name__, "N/A", "N/A", "N/A"))
            print(f"{comp_cls.__name__:<35} | {'N/A':<15} | {'N/A':<15} | {'N/A':<10}")


if __name__ == "__main__":
    main()

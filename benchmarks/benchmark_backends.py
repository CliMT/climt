from datetime import timedelta

import sympl
from pyinstrument import Profiler
from sympl._core.backend import DataArrayBackend

import climt
from climt import EmanuelConvection, UnytBackend, get_default_state, get_grid


def get_timestep_quantity(timestep, backend_obj):
    if isinstance(backend_obj, UnytBackend):
        from climt import UnytTimeDelta

        return UnytTimeDelta(timestep)
    elif isinstance(backend_obj, DataArrayBackend):
        from sympl import DataArray

        return DataArray(timestep.total_seconds(), attrs={"units": "s"})
    return timestep


def run_simulation(backend_obj, iterations=100, component_class=EmanuelConvection):
    print(f"Initializing {component_class.__name__}...")
    component = component_class()

    print("Getting default state...")
    grid = get_grid(nx=32, ny=32, nz=20)
    state = get_default_state([component], grid_state=grid)

    print(f"Running {iterations} iterations...")
    timestep = timedelta(minutes=10)
    dt_quant = get_timestep_quantity(timestep, backend_obj)

    for i in range(iterations):
        tendencies, diagnostics = component(state, timestep)
        for name, tendency in tendencies.items():
            if name in state:
                state[name] = state[name] + tendency * dt_quant
        state.update(diagnostics)
        state["time"] += timestep


def benchmark(backend_name, backend_obj, iterations=100):
    print(f"\n{'=' * 40}")
    print(f"Benchmarking: {backend_name}")
    print(f"{'=' * 40}")

    sympl.set_backend(backend_obj)

    profiler = Profiler()
    profiler.start()

    try:
        run_simulation(backend_obj, iterations=iterations)
    except Exception as e:
        print(f"Simulation failed with {backend_name}: {e}")
        import traceback

        traceback.print_exc()
    finally:
        profiler.stop()

    print(profiler.output_text(unicode=True, color=True, show_all=False))


if __name__ == "__main__":
    # Benchmark DataArrayBackend (Default)
    benchmark("DataArrayBackend (Default)", DataArrayBackend(), iterations=5000)

    # Benchmark UnytBackend
    benchmark("UnytBackend", UnytBackend(), iterations=5000)

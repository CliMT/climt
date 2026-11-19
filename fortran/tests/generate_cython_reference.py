"""Record the inputs and outputs of the Cython climt-fortran calls.

The reference (``fixtures/cython_reference.npz``) is the yardstick that later
builds of climt-fortran must reproduce exactly. Each call is captured by
running the climt component that makes it, with the extension module swapped
for a recorder, so the arguments are exactly the ones climt passes by default.

Keys, for each call name in ``CALLS``:

* ``<call>/init/<function>/<arg>``: the set-up calls the component makes
  before the call (``set_constants``, ``initialise_rrtm_radiation[_mcica]``,
  ``init_emanuel_convection``), in the order given by ``<call>/init_order``.
* ``<call>/function``: name of the extension function called.
* ``<call>/in/<arg>``: every argument, as passed (scalars as 0-d arrays).
* ``<call>/out/<arg>``: every array argument after the call.

Run: ``python fortran/tests/generate_cython_reference.py [output.npz]``.
"""

import datetime
import inspect
import os
import sys
from contextlib import contextmanager

import numpy as np

SEED = 20261119
NUM_COLUMNS = 4
NUM_LEVELS = 30
CLOUD_LEVELS = slice(10, 15)
TIMESTEP = datetime.timedelta(minutes=10)
CALLS = ("emanuel", "rrtmg_lw", "rrtmg_lw_mcica", "rrtmg_sw", "rrtmg_sw_mcica")
DEFAULT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "fixtures", "cython_reference.npz"
)


class _Recorder:
    """Stands in for an extension module, recording every call into ``calls``."""

    def __init__(self, module, calls):
        self._module = module
        self._calls = calls

    def __getattr__(self, name):
        function = getattr(self._module, name)
        names = list(inspect.signature(function).parameters)

        def record(*args):
            before = [np.array(a, order="K", copy=True) for a in args]
            result = function(*args)
            after = {
                n: np.array(a, order="K", copy=True)
                for n, a in zip(names, args)
                if isinstance(a, np.ndarray)
            }
            self._calls.append((name, dict(zip(names, before)), after))
            return result

        return record


@contextmanager
def recording(component_module, attribute):
    """Replace ``component_module.<attribute>`` with a recorder; yield the call list."""
    calls = []
    original = getattr(component_module, attribute)
    setattr(component_module, attribute, _Recorder(original, calls))
    try:
        yield calls
    finally:
        setattr(component_module, attribute, original)


def _make_state(component, num_columns, rng, clouds):
    """A default climt state for ``component`` with fixed pseudo-random fields."""
    import climt

    state = climt.get_default_state(
        [component], grid_state=climt.get_grid(nx=num_columns, nz=NUM_LEVELS)
    )
    shape = state["air_temperature"].shape  # (levels, 1, columns)
    surface_shape = shape[1:]

    p = state["air_pressure"].values
    ps = state["surface_air_pressure"].values
    surface_temperature = 300.0 + rng.uniform(-2.0, 2.0, surface_shape)
    temperature = np.maximum(surface_temperature * (p / ps) ** 0.19, 200.0)
    temperature += rng.normal(0.0, 0.5, shape)
    if "surface_temperature" in state:
        state["surface_temperature"].values[:] = surface_temperature
    state["air_temperature"].values[:] = temperature
    state["specific_humidity"].values[:] = (
        0.018 * (p / ps) ** 3 * rng.uniform(0.8, 1.0, shape)
    )
    if "zenith_angle" in state:
        state["zenith_angle"].values[:] = rng.uniform(0.0, 1.2, surface_shape)
    if "eastward_wind" in state:
        state["eastward_wind"].values[:] = rng.normal(5.0, 3.0, shape)
        state["northward_wind"].values[:] = rng.normal(0.0, 3.0, shape)
    if clouds:
        cloud_shape = state["cloud_area_fraction_in_atmosphere_layer"][
            CLOUD_LEVELS
        ].shape
        state["cloud_area_fraction_in_atmosphere_layer"][CLOUD_LEVELS] = (
            rng.uniform(0.3, 0.9, cloud_shape)
        )
        state["mass_content_of_cloud_liquid_water_in_atmosphere_layer"][
            CLOUD_LEVELS
        ] = rng.uniform(0.01, 0.1, cloud_shape)
        state["mass_content_of_cloud_ice_in_atmosphere_layer"][CLOUD_LEVELS] = (
            rng.uniform(0.0, 0.02, cloud_shape)
        )
    return state


def record_call(name, num_columns=NUM_COLUMNS):
    """Run the component behind call ``name``; return (init calls, main call).

    Each record is ``(function name, inputs by argument, array outputs by argument)``.
    """
    import climt
    from climt._components.emanuel import component as emanuel_module
    from climt._components.rrtmg.lw import component as lw_module
    from climt._components.rrtmg.sw import component as sw_module

    rng = np.random.default_rng(SEED)
    mcica = name.endswith("_mcica")
    if name == "emanuel":
        module, attribute = emanuel_module, "_emanuel_convection"
        factory = climt.EmanuelConvection
    elif name.startswith("rrtmg_lw"):
        module, attribute = lw_module, "_rrtmg_lw"
        factory = lambda: climt.RRTMGLongwave(mcica=mcica)  # noqa: E731
    else:
        module, attribute = sw_module, "_rrtmg_sw"
        factory = lambda: climt.RRTMGShortwave(mcica=mcica)  # noqa: E731

    with recording(module, attribute) as calls:
        component = factory()
        state = _make_state(component, num_columns, rng, clouds=mcica)
        # The MCICA components draw permute_seed from numpy's global generator.
        np.random.seed(SEED)
        if name == "emanuel":
            component(state, TIMESTEP)
        else:
            component(state)
    return calls[:-1], calls[-1]


def generate(path=DEFAULT_PATH):
    """Record every call in ``CALLS`` and save the reference to ``path``."""
    arrays = {}
    for name in CALLS:
        init_calls, (function, inputs, outputs) = record_call(name)
        arrays[name + "/init_order"] = np.array([c[0] for c in init_calls])
        for init_function, init_inputs, _ in init_calls:
            for arg, value in init_inputs.items():
                arrays["%s/init/%s/%s" % (name, init_function, arg)] = value
        arrays[name + "/function"] = np.array(function)
        for arg, value in inputs.items():
            arrays["%s/in/%s" % (name, arg)] = value
        for arg, value in outputs.items():
            arrays["%s/out/%s" % (name, arg)] = value
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    np.savez_compressed(path, **arrays)
    print("Wrote %s with calls: %s" % (path, ", ".join(CALLS)))


if __name__ == "__main__":
    generate(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATH)

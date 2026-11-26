"""An ensemble of GMD radiative-convective columns run as one compiled model.

Each column is one scenario: surface temperature from 295 K to 305 K, and the
Emanuel scheme's ALPHA from 0.1 to 0.3. ``climt.build_column_model`` compiles
the components once and steps all columns on NumPy arrays.
"""
import time
from datetime import timedelta

import numpy as np

import climt


def gmd_column():
    """(state, model kwargs, convection) for the column of gmd_radiative_convective_unyt.py."""
    convection = climt.EmanuelConvectionPython()
    tendency_components = [convection, climt.RRTMGLongwave(), climt.RRTMGShortwave(),
                           climt.SlabSurface()]
    stepper = climt.SimplePhysics()
    state = climt.get_default_state([stepper] + tendency_components)
    values = {
        "air_temperature": 270.,
        "surface_albedo_for_direct_shortwave": 0.5,
        "surface_albedo_for_direct_near_infrared": 0.5,
        "surface_albedo_for_diffuse_shortwave": 0.5,
        "zenith_angle": np.pi / 2.5,
        "surface_temperature": 300.,
        "ocean_mixed_layer_thickness": 5.,
        "area_type": "sea",
        "eastward_wind": 3.,
    }
    for name, value in values.items():
        state[name].values[:] = value
    kwargs = dict(timestep=timedelta(minutes=5), tendency_components=tendency_components,
                  steppers=[stepper], prescribed=["eastward_wind"])
    return state, kwargs, convection


def main(n=64, n_steps=2000):
    """Run ``n`` scenarios for ``n_steps`` steps of 5 minutes; return the RunResult."""
    state, kwargs, convection = gmd_column()
    # Surface temperature rises across the columns; ALPHA cycles through n // 8 values.
    surface_temperature = np.linspace(295., 305., n)
    alpha = np.resize(np.linspace(0.1, 0.3, max(1, n // 8)), n)
    scenarios = climt.Scenarios(
        n, inputs={"surface_temperature": surface_temperature},
        parameters={(convection, "ALPHA"): alpha})
    model = climt.build_column_model(state, scenarios=scenarios, **kwargs)

    start = time.perf_counter()
    result = model.run(n_steps, output_every=max(1, n_steps // 20),
                       outputs=["air_temperature", "surface_temperature"])
    elapsed = time.perf_counter() - start

    final = result.arrays["surface_temperature"][-1]
    print("{} columns, {} steps in {:.1f} s ({} threads)".format(
        n, n_steps, elapsed, model.threads))
    print("final surface temperature: {:.2f} K to {:.2f} K".format(final.min(), final.max()))
    return result


if __name__ == "__main__":
    main()

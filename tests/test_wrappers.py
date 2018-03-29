import climt
import numpy as np
from datetime import timedelta


def test_scaling_prognostic_tendency():

    radiation = climt.RRTMGLongwave()
    radiation = radiation.scaled_version(
        tendency_scale_factors=dict(air_temperature=0))

    state = climt.get_default_state([radiation])

    tendency, diagnostics = radiation(state)

    assert np.all(tendency['air_temperature'].values == 0)


def test_scaling_implicit_output():

    simple_physics = climt.SimplePhysics()

    simple_physics = simple_physics.scaled_version(
        output_scale_factors=dict(air_temperature=0))

    state = climt.get_default_state([simple_physics])

    diagnostics, output = simple_physics(state, timedelta(hours=1))

    assert np.all(output['air_temperature'].values == 0)


def test_scaling_diagnostic_output():

    dcmip = climt.DcmipInitialConditions()
    dcmip = dcmip.scaled_version(
        diagnostic_scale_factors=dict(eastward_wind=0))

    state = climt.get_default_state([dcmip])

    diagnostics = dcmip(state)

    assert np.all(diagnostics['eastward_wind'].values == 0)

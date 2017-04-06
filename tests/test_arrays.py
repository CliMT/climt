from climt import (
    create_state_dict_for, get_numpy_arrays_from_state,
    get_default_state, RRTMGLongwave
)

import pytest
import numpy as np

from .test_classes import (
    MockPrognostic,
    MockPrognosticWithExtraQuantities,
    MockPrognosticWithExtraQuantitiesNotDefined,
    MockPrognosticWithAllAttributes
)


def test_get_diagnostics():

    dummy = MockPrognosticWithAllAttributes()
    state = get_default_state([dummy])

    diag = create_state_dict_for(dummy, 'diagnostics', state)

    for quantity in dummy.diagnostics.keys():
        assert quantity in diag


def test_get_diagnostics_with_real_component():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy])

    diag = create_state_dict_for(dummy, 'diagnostics', state)

    for quantity in dummy.diagnostics.keys():
        assert quantity in diag


def test_get_inputs_with_real_component():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy])

    diag = create_state_dict_for(dummy, 'inputs', state)

    for quantity in dummy.inputs.keys():
        assert quantity in diag
        assert diag[quantity].dims == state[quantity].dims


def test_get_tendencies():

    dummy = MockPrognosticWithAllAttributes()
    state = get_default_state([dummy])

    diag = create_state_dict_for(dummy, 'tendencies', state)

    for quantity in dummy.tendencies.keys():
        assert quantity in diag


def test_get_outputs():

    dummy = MockPrognosticWithAllAttributes()
    state = get_default_state([dummy])

    diag = create_state_dict_for(dummy, 'outputs', state)

    for quantity in dummy.outputs.keys():
        assert quantity in diag


def test_get_diagnostics_with_real_component_with_2d_coordinates():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy],
                              x=dict(label='shore', values=np.random.randn(2, 2), units='km'),
                              y=dict(label='latitude', values=np.random.randn(2, 2), units='degrees east'))

    diag = create_state_dict_for(dummy, 'diagnostics', state)

    for quantity in dummy.diagnostics.keys():
        assert quantity in diag
        assert 'shore' in diag[quantity].coords
        assert 'latitude' in diag[quantity].coords
        assert diag[quantity].coords['shore'].ndim == 2
        assert diag[quantity].coords['latitude'].ndim == 2


def test_extracting_arrays_from_real_component():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy])

    arrays = get_numpy_arrays_from_state(dummy, 'inputs', state)

    for quantity in dummy.inputs.keys():
        units = dummy.inputs[quantity]
        state_values = state[quantity].to_units(units).values
        assert np.all(arrays[quantity] == state_values)


def test_inputs_is_not_dict():

    dummy = MockPrognostic()
    state = get_default_state([dummy])
    dummy.inputs = ('air_temperature', 'oxygen_mixing_ratio')

    with pytest.raises(NotImplementedError) as excinfo:
        get_numpy_arrays_from_state(dummy, 'inputs', state)
    assert 'with a dict-like' in str(excinfo.value)


def test_c_memory_layout():

    dummy = MockPrognostic()
    state = get_default_state([dummy], x=dict(label='longitude',
                                              values=np.arange(10),
                                              units='degrees_east'))

    get_numpy_arrays_from_state(dummy, 'inputs', state, memory_layout='c')


def test_wrong_memory_layout():

    dummy = MockPrognostic()
    state = get_default_state([dummy])

    with pytest.raises(ValueError) as excinfo:
        get_numpy_arrays_from_state(dummy, 'inputs', state, memory_layout='abcd')
    assert 'memory_layout' in str(excinfo.value)


def test_unknown_quantity_in_component():

    dummy = MockPrognostic()
    state = get_default_state([dummy])
    state.pop('air_temperature')

    with pytest.raises(IndexError) as excinfo:
        get_numpy_arrays_from_state(dummy, 'inputs', state, memory_layout='c')
    assert 'does not contain' in str(excinfo.value)


def test_get_class_defined_quantity():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])
    input_arrays = get_numpy_arrays_from_state(dummy, 'inputs', state)

    assert np.all(
        state['sigma_on_interface_levels'].values == input_arrays['sigma_on_interface_levels'])


def test_get_undefined_array():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])
    dummy = MockPrognosticWithExtraQuantitiesNotDefined()

    with pytest.raises(IndexError) as excinfo:
        get_numpy_arrays_from_state(dummy, 'inputs', state)
    assert 'not described' in str(excinfo.value)


def test_unknown_attribute():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])

    with pytest.raises(IndexError) as excinfo:
        get_numpy_arrays_from_state(dummy, 'random', state)
    assert 'no attribute called' in str(excinfo.value)

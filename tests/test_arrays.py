from climt import (
    get_default_state, RRTMGLongwave
)

import pytest
import numpy as np

from .test_classes import (
    MockPrognostic,
    MockPrognosticWithExtraQuantities,
    MockPrognosticWithExtraQuantitiesNotDefined,
    MockPrognosticWithAllAttributes,
    MockImplicitWithAllAttributes
)


def test_get_diagnostics():

    dummy = MockPrognosticWithAllAttributes()
    state = get_default_state([dummy])

    diag = dummy.create_state_dict_for('_climt_diagnostics', state)

    for quantity in dummy.diagnostics:
        assert quantity in diag


def test_get_diagnostics_with_real_component():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy])

    diag = dummy.create_state_dict_for('_climt_diagnostics', state)

    for quantity in dummy._climt_diagnostics:
        assert quantity in diag


def test_get_inputs_with_real_component():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy])

    diag = dummy.create_state_dict_for('_climt_inputs', state)

    for quantity in dummy._climt_inputs:
        assert quantity in diag
        assert diag[quantity].dims == state[quantity].dims


def test_get_tendencies():

    dummy = MockPrognosticWithAllAttributes()
    state = get_default_state([dummy])

    diag = dummy.create_state_dict_for('_climt_tendencies', state)

    for quantity in dummy.tendencies:
        assert quantity in diag


def test_get_outputs():

    dummy = MockImplicitWithAllAttributes()
    state = get_default_state([dummy])

    diag = dummy.create_state_dict_for('_climt_outputs', state)

    for quantity in dummy.outputs:
        assert quantity in diag


def test_get_diagnostics_with_real_component_with_2d_coordinates():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy],
                              x=dict(label='shore', values=np.random.randn(2, 2), units='km'),
                              y=dict(label='latitude', values=np.random.randn(2, 2), units='degrees east'))

    diag = dummy.create_state_dict_for('_climt_diagnostics', state)

    for quantity in dummy.diagnostics:
        assert quantity in diag
        assert 'shore' in diag[quantity].coords
        assert 'latitude' in diag[quantity].coords
        assert diag[quantity].coords['shore'].ndim == 2
        assert diag[quantity].coords['latitude'].ndim == 2


def test_extracting_arrays_from_real_component():
    dummy = RRTMGLongwave()
    state = get_default_state([dummy])

    arrays = dummy.get_numpy_arrays_from_state('_climt_inputs', state)

    for quantity in dummy.inputs:
        units = dummy._climt_inputs[quantity]
        state_values = state[quantity].to_units(units).values
        assert np.all(arrays[quantity] == state_values)


def test_inputs_is_not_dict():

    dummy = MockPrognostic()
    state = get_default_state([dummy])
    dummy._climt_inputs = ('air_temperature', 'oxygen_mixing_ratio')

    with pytest.raises(NotImplementedError) as excinfo:
        dummy.get_numpy_arrays_from_state('_climt_inputs', state)
    assert 'with a dict-like' in str(excinfo.value)


def test_c_memory_layout():

    dummy = MockPrognostic()
    state = get_default_state([dummy], x=dict(label='longitude',
                                              values=np.arange(10),
                                              units='degrees_east'))

    dummy.get_numpy_arrays_from_state('_climt_inputs', state, memory_layout='c')


def test_wrong_memory_layout():

    dummy = MockPrognostic()
    state = get_default_state([dummy])

    with pytest.raises(ValueError) as excinfo:
        dummy.get_numpy_arrays_from_state('_climt_inputs', state, memory_layout='abcd')
    assert 'memory_layout' in str(excinfo.value)


def test_unknown_quantity_in_component():

    dummy = MockPrognostic()
    state = get_default_state([dummy])
    state.pop('air_temperature')

    with pytest.raises(IndexError) as excinfo:
        dummy.get_numpy_arrays_from_state('_climt_inputs', state, memory_layout='c')
    assert 'does not contain' in str(excinfo.value)


def test_get_class_defined_quantity():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])
    input_arrays = dummy.get_numpy_arrays_from_state('_climt_inputs', state)

    assert np.all(
        state['sigma_on_interface_levels'].values == input_arrays['sigma_on_interface_levels'])


def test_get_undefined_array():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])
    dummy = MockPrognosticWithExtraQuantitiesNotDefined()

    with pytest.raises(IndexError) as excinfo:
        dummy.get_numpy_arrays_from_state('_climt_inputs', state)
    assert 'not described' in str(excinfo.value)


def test_unknown_attribute():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])

    with pytest.raises(IndexError) as excinfo:
        dummy.get_numpy_arrays_from_state('random', state)
    assert 'no attribute called' in str(excinfo.value)

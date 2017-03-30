from climt import (
    get_default_state)

from sympl import SharedKeyException
import numpy as np
import pytest
from .test_classes import (MockPrognostic,
                           MockPrognosticWithExtraDimensions,
                           MockPrognosticWithExtraQuantities)


def test_no_components():

    with pytest.raises(ValueError) as excinfo:
        get_default_state([])
    assert 'at least one' in str(excinfo.value)


def test_input_state_has_overlapping_keys():

    dummy = MockPrognostic()
    with pytest.raises(SharedKeyException):
        get_default_state([dummy], input_state={'air_temperature': 0})


def test_basic_case_for_two_inputs():

    dummy = MockPrognostic()
    state = get_default_state([dummy])

    required_quantities = list(dummy.inputs.keys())
    required_quantities.extend(['longitude',
                               'latitude',
                                'mid_levels'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['oxygen_mixing_ratio'].dims == ('longitude', 'latitude', 'mid_levels')


def test_case_for_x_dim_defined():

    dummy = MockPrognostic()
    state = get_default_state([dummy], x=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_east'))

    required_quantities = list(dummy.inputs.keys())
    required_quantities.extend(['along_shore',
                               'latitude',
                                'mid_levels'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['air_temperature'].dims == ('along_shore', 'latitude', 'mid_levels')


def test_case_for_y_dim_defined():

    dummy = MockPrognostic()
    state = get_default_state([dummy], y=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_north'))

    required_quantities = list(dummy.inputs.keys())
    required_quantities.extend(['longitude',
                               'along_shore',
                                'mid_levels'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['air_temperature'].dims == ('longitude', 'along_shore', 'mid_levels')


def test_case_for_z_dim_defined():

    dummy = MockPrognostic()
    state = get_default_state([dummy], z=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_east'))

    required_quantities = list(dummy.inputs.keys())
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['air_temperature'].dims == ('longitude', 'latitude', 'mid_levels')


def test_with_extra_dimensions():

    dummy = MockPrognosticWithExtraDimensions()
    state = get_default_state([dummy], z=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_east'))

    required_quantities = list(dummy.inputs.keys())
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore',
                                'some_other_dimension'])

    for quantity in state.keys():
        assert quantity in required_quantities


def test_with_extra_quantities():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy], z=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_east'))

    required_quantities = list(dummy.inputs.keys())
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore',
                                'some_quantity'])

    for quantity in state.keys():
        assert quantity in required_quantities


if __name__ == '__main__':
    pytest.main([__file__])

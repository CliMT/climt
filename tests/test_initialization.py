from climt import (
    get_default_state)

from sympl import SharedKeyError
import numpy as np
import pytest
from .test_classes import (MockPrognostic,
                           MockPrognosticWithExtraDimensions,
                           MockPrognosticWithExtraQuantities,
                           MockPrognosticWithMalformedExtraQuantities,
                           MockPrognosticWithExtraDimensionsAndSigmaLevels,
                           MockPrognosticWithExtraDimensionsIn2d)


def test_no_components():

    with pytest.raises(ValueError) as excinfo:
        get_default_state([])
    assert 'at least one' in str(excinfo.value)


def test_input_state_has_overlapping_keys():

    dummy = MockPrognostic()
    with pytest.raises(SharedKeyError):
        get_default_state([dummy], initial_state={'air_temperature': 0})


def test_only_mid_levels_specified():

    dummy = MockPrognostic()
    with pytest.raises(ValueError):
        get_default_state([dummy], mid_levels=dict(
            label='vert_coord', values=np.linspace(0, 10, 10),
            units='kilometer'))


def test_only_int_levels_specified():

    dummy = MockPrognostic()
    with pytest.raises(ValueError):
        get_default_state([dummy], interface_levels=dict(
            label='vert_coord', values=np.linspace(0, 10, 10),
            units='kilometer'))


def test_mid_level_2d():

    dummy = MockPrognostic()
    with pytest.raises(ValueError):
        get_default_state([dummy], mid_levels=dict(
            label='vert_coord', values=np.random.randn(4, 4),
            units='kilometer'), interface_levels=dict(
                label='vert_coord', values=np.linspace(0, 10, 10),
                units='kilometer'))


def test_int_level_2d():

    dummy = MockPrognostic()
    with pytest.raises(ValueError):
        get_default_state([dummy],
                          interface_levels=dict(
                              label='vert_coord',
                              values=np.random.randn(4, 4),
                              units='kilometer'),
                          mid_levels=dict(
                              label='vert_coord',
                              values=np.linspace(0, 10, 10),
                              units='kilometer'))


def test_mid_int_levels_not_consistent():

    dummy = MockPrognostic()
    with pytest.raises(ValueError):
        get_default_state([dummy], interface_levels=dict(
            label='vert_coord', values=np.linspace(0, 10, 10),
            units='kilometer'), mid_levels=dict(
            label='vert_coord', values=np.linspace(0, 10, 11),
            units='kilometer'))


def test_basic_case_for_two_inputs():

    dummy = MockPrognostic()
    state = get_default_state([dummy])

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['longitude',
                               'latitude',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['mole_fraction_of_oxygen_in_air'].dims == ('longitude', 'latitude', 'mid_levels')


def test_case_for_x_dim_defined():

    dummy = MockPrognostic()
    state = get_default_state([dummy], x=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_east'))

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['along_shore',
                               'latitude',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['air_temperature'].dims == ('along_shore', 'latitude', 'mid_levels')


def test_case_for_y_dim_defined():

    dummy = MockPrognostic()
    state = get_default_state([dummy], y=dict(
        label='along_shore', values=np.linspace(0, 10, 10),
        units='degrees_north'))

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['longitude',
                               'along_shore',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['air_temperature'].dims == ('longitude', 'along_shore', 'mid_levels')


def test_case_for_z_dim_defined():

    dummy = MockPrognostic()
    state = get_default_state([dummy],
                              mid_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 10, 10),
                                  units='degrees_east'),
                              interface_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 11, 11),
                                  units='degrees_east'))

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities

    assert state['air_temperature'].dims == ('longitude', 'latitude', 'mid_levels')


def test_with_extra_dimensions():

    dummy = MockPrognosticWithExtraDimensions()
    state = get_default_state([dummy],
                              mid_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 10, 10),
                                  units='degrees_east'),
                              interface_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 11, 11),
                                  units='degrees_east'))

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore',
                                'some_other_dimension',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities


def test_with_extra_dimensions_and_sigma_levels():

    dummy = MockPrognosticWithExtraDimensionsAndSigmaLevels()
    state = get_default_state([dummy],
                              mid_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 10, 10),
                                  units='degrees_east'),
                              interface_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 11, 11),
                                  units='degrees_east'))

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore',
                                'some_other_dimension',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities


def test_with_extra_quantities():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy],
                              mid_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 10, 10),
                                  units='degrees_east'),
                              interface_levels=dict(
                                  label='along_shore',
                                  values=np.linspace(0, 11, 11),
                                  units='degrees_east'))

    required_quantities = list(dummy.inputs)
    required_quantities.extend(['longitude',
                               'latitude',
                                'along_shore',
                                'some_quantity',
                                'mid_levels',
                                'x',
                                'y',
                                'interface_levels', 'time'])

    for quantity in state.keys():
        assert quantity in required_quantities


def test_with_malformed_extra_quantity():

    dummy = MockPrognosticWithMalformedExtraQuantities()

    with pytest.raises(ValueError) as excinfo:
        get_default_state([dummy])
    assert 'Malformed' in str(excinfo.value)


def test_different_dimension_units():

    dummy = MockPrognostic()
    with pytest.raises(ValueError) as excinfo:
        get_default_state([dummy], y=dict(
            label='along_shore', values=np.ones((2, 2)),
            units='degrees_north'))
    assert 'must have the same shape' in str(excinfo.value)


def test_basic_2d_coordinates():

    dummy = MockPrognosticWithExtraQuantities()
    random_x_values = np.random.randn(3, 4)
    random_y_values = np.random.randn(3, 4)
    state = get_default_state([dummy],
                              x=dict(label='shore', values=random_x_values, units='km'),
                              y=dict(label='latitude', values=random_y_values, units='degrees east'))

    assert state['x'].values.shape[0] == 3
    assert state['x'].values.shape[1] == 4

    assert state['y'].values.shape[0] == 3
    assert state['y'].values.shape[1] == 4

    assert np.all(state['x'].values == random_x_values)
    assert np.all(state['y'].values == random_y_values)

    for quantity in dummy.inputs:
        assert 'shore' in state[quantity].coords
        assert 'latitude' in state[quantity].coords
        assert state[quantity].coords['shore'].units is 'km'
        assert state[quantity].coords['latitude'].units is 'degrees east'


def test_2d_coordinates_wrong_shape():

    dummy = MockPrognosticWithExtraQuantities()

    with pytest.raises(ValueError) as excinfo:
        get_default_state([dummy],
                          x=dict(label='shore', values=np.random.randn(3, 4), units='km'),
                          y=dict(label='latitude', values=np.random.randn(3, 3), units='degrees east'))
    assert '2d coordinates, they' in str(excinfo.value)


def test_2d_coordinates_in_extra_dimensions():

    dummy = MockPrognosticWithExtraDimensionsIn2d()

    with pytest.raises(NotImplementedError) as excinfo:
        get_default_state([dummy])
    assert 'not yet supported' in str(excinfo.value)


if __name__ == '__main__':
    pytest.main([__file__])

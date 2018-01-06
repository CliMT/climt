from climt import (
    get_constant, constant_library)

from copy import deepcopy
import pytest

sample_constants = dict(
    constant_one=dict(
        type='oceanographic_constants',
        value=10.0,
        units='J kg^-1 m^-2'),
    constant_two=dict(
        type='miscellaneous_constants',
        value=1.002,
        units='J kg^-1'),
    constant_three=dict(
        type='stellar_constants',
        value=1.024e6,
        units='kg^-1'))


def test_adding_one_constant():

    constant = {'constant_one':
                sample_constants['constant_one']}

    constant_library.add_constants_to_library(constant)

    assert 'constant_one' in constant_library.current_constants.keys()
    assert 'constant_one' in constant_library.constant_categories['oceanographic_constants']


def test_resetting_library():

    constant = {'constant_two':
                sample_constants['constant_two']}

    constant_library.add_constants_to_library(constant)

    constant_library.reset_constant_library()

    for name in ['constant_one', 'constant_two', 'constant_three']:
        assert name not in constant_library.current_constants.keys()


def test_adding_more_constants():

    constant_library.add_constants_to_library(sample_constants)

    assert 'constant_one' in constant_library.current_constants.keys()
    assert 'constant_one' in constant_library.constant_categories['oceanographic_constants']
    assert 'constant_two' in constant_library.current_constants.keys()
    assert 'constant_two' in constant_library.constant_categories['miscellaneous_constants']
    assert 'constant_three' in constant_library.current_constants.keys()
    assert 'constant_three' in constant_library.constant_categories['stellar_constants']


def test_list_constants():

    constant_library.list_available_constants()


def test_get_non_existent_constant():

    constant_library.reset_constant_library()
    with pytest.raises(IndexError) as excinfo:
        get_constant('constant_one')

    assert 'not present' in str(excinfo.value)


def test_modify_constant():

    constant = {'gravitational_acceleration': {'value': 11.0,
                                               'units': 'km s^-2',
                                               'type': 'physical_constants'}}

    constant_library.modify_constants_in_library(constant)

    gravity = get_constant('gravitational_acceleration')

    assert gravity.values.item() == 11.0
    assert gravity.units == 'km s^-2'

    constant_library.reset_constant_library()

    gravity = get_constant('gravitational_acceleration')

    assert gravity.values.item() == 9.80665
    assert gravity.units == 'm s^-2'


def test_adding_existing_constant():

    constant = {'gravitational_acceleration': {'value': 11.0,
                                               'units': 'km s^-2',
                                               'type': 'physical_constants'}}

    with pytest.raises(IndexError) as excinfo:
        constant_library.add_constants_to_library(constant)

    assert 'already present' in str(excinfo.value)


def test_adding_constant_with_wrong_type():

    constant = {'constant_one': deepcopy(sample_constants['constant_one'])}
    constant['constant_one']['type'] = 'xyz'

    with pytest.raises(IndexError) as excinfo:
        constant_library.add_constants_to_library(constant)

    assert 'valid category' in str(excinfo.value)


def test_modify_constant_which_does_not_exist():

    constant = {'constant_one': deepcopy(sample_constants['constant_one'])}
    constant['constant_one']['type'] = 'xyz'

    with pytest.raises(IndexError) as excinfo:
        constant_library.modify_constants_in_library(constant)

    assert 'does not exist' in str(excinfo.value)

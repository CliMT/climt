from climt import (
    mass_to_volume_mixing_ratio,
    get_input_arrays_from_state,
    get_default_state)

from sympl import DataArray
import numpy as np
import pytest
from .test_classes import MockPrognostic, MockPrognosticWithExtraQuantities


def test_mol_weight_not_passed():

    with pytest.raises(ValueError) as excinfo:
        mass_to_volume_mixing_ratio(None)

    assert 'molecular weight' in str(excinfo.value)


def test_for_365ppm():

    mass_mixing_ratio_co2 = 605e-6
    molecular_weight_co2 = 48

    expected_vol_mixing_ratio = 365.0670833333333e-6

    co2 = DataArray(
        mass_mixing_ratio_co2*np.ones((1, 1, 1)),
        dims=['longitude', 'latitude', 'mid_levels'],
        attrs=dict(units='g/g'))

    vol_mixing_ratio = mass_to_volume_mixing_ratio(
        co2, molecular_weight=molecular_weight_co2)

    assert np.all(np.isclose(
        vol_mixing_ratio.values, expected_vol_mixing_ratio))

    assert vol_mixing_ratio.dims == co2.dims


def test_for_g_per_kg():

    mass_mixing_ratio_co2 = 605e-3
    molecular_weight_co2 = 48

    expected_vol_mixing_ratio = 365.0670833333333e-6

    co2 = DataArray(
        mass_mixing_ratio_co2*np.ones((1, 1, 1)),
        dims=['longitude', 'latitude', 'mid_levels'],
        attrs=dict(units='g/kg'))

    vol_mixing_ratio = mass_to_volume_mixing_ratio(
        co2, molecular_weight=molecular_weight_co2)

    assert np.all(np.isclose(
        vol_mixing_ratio.values, expected_vol_mixing_ratio))


def test_inputs_is_not_dict():

    dummy = MockPrognostic()
    state = get_default_state([dummy])
    dummy.inputs = ('air_temperature', 'oxygen_mixing_ratio')

    with pytest.raises(NotImplementedError) as excinfo:
        get_input_arrays_from_state(dummy, state)
    assert 'with dict-like' in str(excinfo.value)


def test_c_memory_layout():

    dummy = MockPrognostic()
    state = get_default_state([dummy], x=dict(label='longitude',
                                              values=np.arange(10),
                                              units='degrees_east'))

    get_input_arrays_from_state(dummy, state, memory_layout='c')


def test_wrong_memory_layout():

    dummy = MockPrognostic()
    state = get_default_state([dummy])

    with pytest.raises(ValueError) as excinfo:
        get_input_arrays_from_state(dummy, state, memory_layout='abcd')
    assert 'memory_layout' in str(excinfo.value)


def test_unknown_quantity_in_component():

    dummy = MockPrognostic()
    state = get_default_state([dummy])
    state.pop('air_temperature')

    with pytest.raises(IndexError) as excinfo:
        get_input_arrays_from_state(dummy, state, memory_layout='c')
    assert 'does not contain' in str(excinfo.value)


def test_get_class_defined_quantity():

    dummy = MockPrognosticWithExtraQuantities()
    state = get_default_state([dummy])
    input_arrays = get_input_arrays_from_state(dummy, state)

    assert np.all(
        state['sigma_on_interface_levels'].values == input_arrays['sigma_on_interface_levels'])


if __name__ == '__main__':
    pytest.main([__file__])

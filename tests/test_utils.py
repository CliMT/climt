from climt import (
    mass_to_volume_mixing_ratio,
    get_interface_values,
    calculate_q_sat, numpy_version_of,
    RRTMGShortwave, get_default_state,
    set_constants_from_dict,
    list_available_constants
)

from sympl import DataArray, get_constant
import numpy as np
import pytest

sample_constants = dict(
    constant_one=dict(
        value=10.0,
        units='J kg^-1 m^-2'),
    constant_two=dict(
        value=1.002,
        units='J kg^-1'),
    constant_three=dict(
        value=1.024e6,
        units='kg^-1'))


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
        vol_mixing_ratio, expected_vol_mixing_ratio))


def test_for_g_per_kg():

    mass_mixing_ratio_co2 = 605e-3
    molecular_weight_co2 = 48

    expected_vol_mixing_ratio = [[[365.0670833333333e-6]]]

    co2 = DataArray(
        mass_mixing_ratio_co2*np.ones((1, 1, 1)),
        dims=['longitude', 'latitude', 'mid_levels'],
        attrs=dict(units='g/kg'))

    vol_mixing_ratio = mass_to_volume_mixing_ratio(
        co2.to_units('g/g').values, molecular_weight=molecular_weight_co2)

    print(vol_mixing_ratio, expected_vol_mixing_ratio)
    assert np.all(np.isclose(
        vol_mixing_ratio, expected_vol_mixing_ratio))


def test_interface_levels():

    mid_level_values = np.ones((10, 1))
    surface_value = np.ones((1, 1))

    pressure_mid_level = np.linspace(0.995, 0.001, 10)[:, None]
    surface_pressure = 1.

    pressure_interface_level = np.zeros((11, 1))

    pressure_interface_level[1:-1, :] = (
        pressure_mid_level[1:, :] + pressure_mid_level[:-1, :])/2.

    pressure_interface_level[0, :] = surface_pressure
    pressure_interface_level[-1, :] = 0.0005

    interface_values = get_interface_values(
        mid_level_values, surface_value,
        pressure_mid_level, pressure_interface_level)

    assert np.all(interface_values == np.ones((11, 1)))


def test_qsat():

    surf_temp = 290*np.random.randn(10, 10)
    surf_temp[surf_temp < 260] = 260
    surf_temp[surf_temp > 280] = 280

    qsat_at_280 = 0.0062856560708380816
    qsat_at_260 = 0.0012408979354134003

    surf_press = 1e5*np.ones(surf_temp.shape)

    Rd = 287.058  # J/kg/K
    Rv = 461.5  # J/kg/K

    qsat = calculate_q_sat(surf_temp, surf_press, Rd, Rv)

    assert np.all(qsat[surf_temp == 280] == qsat_at_280)
    assert np.all(qsat[surf_temp == 260] == qsat_at_260)


def test_get_numpy_version_with_numpy_array_in_state():

    component = RRTMGShortwave()
    state = get_default_state([component])

    state['test_values'] = np.arange(100)

    raw_array = numpy_version_of(state)

    assert np.all(raw_array['test_values'] == np.arange(100))


def test_adding_constant():

    set_constants_from_dict(sample_constants)

    for constant in sample_constants.keys():

        constant_value = get_constant(constant,
                                      sample_constants[constant]['units'])

        assert constant_value == sample_constants[constant]['value']


def test_print_constant():

    list_available_constants()


if __name__ == '__main__':
    pytest.main([__file__])

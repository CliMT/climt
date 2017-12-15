from sympl import default_constants, DataArray
from copy import deepcopy

climt_constants = deepcopy(default_constants)

constant_categories = dict(
    planetary_constants=('gravitational_acceleration',
                         'planetary_radius',
                         'planetary_rotation_rate',
                         'seconds_per_day'),

    physical_constants=('stefan_boltzmann_constant',
                        'avogadro_constant',
                        'speed_of_light',
                        'boltzmann_constant',
                        'loschmidt_constant',
                        'universal_gas_constant',
                        'planck_constant'),

    condensible_constants=('density_of_liquid_phase',
                           'heat_capacity_of_liquid_phase',
                           'heat_capacity_of_vapor_phase',
                           'specific_enthalpy_of_vapor_phase',
                           'gas_constant_of_vapor_phase',
                           'latent_heat_of_condensation',
                           'latent_heat_of_fusion',
                           'density_of_solid_phase_as_ice',
                           'density_of_solid_phase_as_snow',
                           'heat_capacity_of_solid_phase_as_ice',
                           'heat_capacity_of_solid_phase_as_snow',
                           'thermal_conductivity_of_solid_phase_as_ice',
                           'thermal_conductivity_of_solid_phase_as_snow',
                           'thermal_conductivity_of_liquid_phase',
                           'freezing_temperature_of_liquid_phase',
                           ),

    atmospheric_constants=('heat_capacity_of_dry_air_at_constant_pressure',
                           'gas_constant_of_dry_air',
                           'thermal_conductivity_of_dry_air',
                           'reference_pressure'),

    stellar_constants=('stellar_irradiance',),

    oceanographic_constants=(),

    chemical_constants=('heat_capacity_of_water_vapor_at_constant_pressure',
                        'density_of_liquid_water',
                        'gas_constant_of_water_vapor',
                        'latent_heat_of_vaporization_of_water',
                        'heat_capacity_of_liquid_water',
                        'latent_heat_of_fusion_of_water',),

    miscellaneous_constants=()
)


def _list_constants(constant_list):

    for quantity in constant_list:
        print("{}: {} {}".format(quantity,
                                 climt_constants[quantity].values,
                                 climt_constants[quantity].units))
        print()


def list_available_constants():

    print('Physical Constants')
    print('===========================')
    _list_constants(constant_categories['physical_constants'])
    print()
    print('Chemical Constants')
    print('===========================')
    _list_constants(constant_categories["chemical_constants"])
    print()
    print('Planetary Constants')
    print('===========================')
    _list_constants(constant_categories["planetary_constants"])
    print()
    print('Stellar Constants')
    print('===========================')
    _list_constants(constant_categories["stellar_constants"])
    print()
    print('Atmospheric Constants')
    print('===========================')
    _list_constants(constant_categories["atmospheric_constants"])
    print()
    print('Condensible Constants')
    print('===========================')
    _list_constants(constant_categories["condensible_constants"])
    print()
    print('Oceanographic Constants')
    print('===========================')
    _list_constants(constant_categories["oceanographic_constants"])
    print()
    print('Miscellaneous Constants')
    print('===========================')
    _list_constants(constant_categories["miscellaneous_constants"])
    print()


def _init_constant_library():

    climt_constants['gas_constant_of_vapor_phase'] = deepcopy(
        climt_constants['gas_constant_of_water_vapor'])

    climt_constants['heat_capacity_of_vapor_phase'] = deepcopy(
        climt_constants['heat_capacity_of_water_vapor_at_constant_pressure'])

    climt_constants['specific_enthalpy_of_vapor_phase'] = deepcopy(
        climt_constants['specific_enthalpy_of_water_vapor'])

    climt_constants['stellar_irradiance'] = deepcopy(
        climt_constants['solar_constant'])

    climt_constants['density_of_liquid_phase'] = deepcopy(
        climt_constants['density_of_liquid_water'])

    climt_constants['heat_capacity_of_liquid_phase'] = deepcopy(
        climt_constants['heat_capacity_of_liquid_water'])

    climt_constants['latent_heat_of_condensation'] = deepcopy(
        climt_constants['latent_heat_of_vaporization_of_water'])

    climt_constants['freezing_temperature_of_liquid_phase'] = deepcopy(
        climt_constants['freezing_temperature_of_liquid_water'])

    climt_constants['density_of_solid_phase_as_ice'] = deepcopy(
        climt_constants['density_of_ice'])

    climt_constants['density_of_solid_phase_as_snow'] = deepcopy(
        climt_constants['density_of_snow'])

    climt_constants['heat_capacity_of_solid_phase_as_ice'] = deepcopy(
        climt_constants['heat_capacity_of_ice'])

    climt_constants['heat_capacity_of_solid_phase_as_snow'] = deepcopy(
        climt_constants['heat_capacity_of_snow'])

    climt_constants['latent_heat_of_fusion'] = deepcopy(
        climt_constants['latent_heat_of_fusion_of_water'])

    climt_constants['freezing_temperature_of_liquid_phase'] = deepcopy(
        climt_constants['freezing_temperature_of_liquid_water'])

    climt_constants['thermal_conductivity_of_solid_phase_as_ice'] = deepcopy(
        climt_constants['thermal_conductivity_of_ice'])

    climt_constants['thermal_conductivity_of_solid_phase_as_snow'] = deepcopy(
        climt_constants['thermal_conductivity_of_snow'])

    climt_constants['thermal_conductivity_of_liquid_phase'] = deepcopy(
        climt_constants['thermal_conductivity_of_liquid_water'])


def add_constant_to_library(name, value, units):
    """
    Add new constants to the library.

    Args:

        name (string):
            name of the constant.

        value (object):
            The value assigned.

        units (string):
            The units, e.g, m/s, J/kg.
    """

    climt_constants[name] = DataArray(value, attrs={'units': units})


def modify_constant_in_library(name, value, units):
    """
    Modify existing constants in the library.

    Args:

        name (string):
            name of the constant.

        value (object):
            The value assigned.

        units (string):
            The units, e.g, m/s, J/kg.

    Raises:

        IndexError: If the constant does not exist in the library.


    """

    try:
        assert name in climt_constants.keys()
    except AssertionError:
        IndexError("{} not present in the library!".format(name))

    climt_constants[name] = DataArray(value, attrs={'units': units})


def get_constant(name):
    """
    Retrieve constant from the library.

    Args:

        name(string):
            The name of the constant to be retrieved.

    Returns:
        constant(DataArray):
            The constant in a DataArray.

    Raises:

        IndexError: If the constant does not exist in the library.
    """

    try:
        assert name in climt_constants.keys()
    except AssertionError:
        IndexError("{} not present in the library!".format(name))

    return climt_constants[name]

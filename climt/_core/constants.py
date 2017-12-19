from sympl import default_constants, DataArray
from copy import deepcopy

constant_name_map = {

    'gas_constant_of_vapor_phase': 'gas_constant_of_water_vapor',
    'heat_capacity_of_vapor_phase': 'heat_capacity_of_water_vapor_at_constant_pressure',
    'specific_enthalpy_of_vapor_phase': 'specific_enthalpy_of_water_vapor',
    'stellar_irradiance': 'solar_constant',
    'density_of_liquid_phase': 'density_of_liquid_water',
    'heat_capacity_of_liquid_phase': 'heat_capacity_of_liquid_water',
    'latent_heat_of_condensation': 'latent_heat_of_vaporization_of_water',
    'freezing_temperature_of_liquid_phase': 'freezing_temperature_of_liquid_water',
    'density_of_solid_phase_as_ice': 'density_of_ice',
    'density_of_solid_phase_as_snow': 'density_of_snow',
    'heat_capacity_of_solid_phase_as_ice': 'heat_capacity_of_ice',
    'heat_capacity_of_solid_phase_as_snow': 'heat_capacity_of_snow',
    'latent_heat_of_fusion': 'latent_heat_of_fusion_of_water',
    'freezing_temperature_of_liquid_phase': 'freezing_temperature_of_liquid_water',
    'thermal_conductivity_of_solid_phase_as_ice': 'thermal_conductivity_of_ice',
    'thermal_conductivity_of_solid_phase_as_snow': 'thermal_conductivity_of_snow',
    'thermal_conductivity_of_liquid_phase': 'thermal_conductivity_of_liquid_water'
}

default_categories = dict(
    planetary_constants=set(('gravitational_acceleration',
                             'planetary_radius',
                             'planetary_rotation_rate',
                             'seconds_per_day')),

    physical_constants=set(('stefan_boltzmann_constant',
                            'avogadro_constant',
                            'speed_of_light',
                            'boltzmann_constant',
                            'loschmidt_constant',
                            'universal_gas_constant',
                            'planck_constant')),

    condensible_constants=set(('density_of_liquid_phase',
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
                               'freezing_temperature_of_liquid_phase',)),

    atmospheric_constants=set(('heat_capacity_of_dry_air_at_constant_pressure',
                               'gas_constant_of_dry_air',
                               'thermal_conductivity_of_dry_air',
                               'reference_air_pressure')),

    stellar_constants=set(('stellar_irradiance',)),

    oceanographic_constants=set(()),

    chemical_constants=set(('heat_capacity_of_water_vapor_at_constant_pressure',
                            'density_of_liquid_water',
                            'gas_constant_of_water_vapor',
                            'latent_heat_of_vaporization_of_water',
                            'heat_capacity_of_liquid_water',
                            'latent_heat_of_fusion_of_water',)),

    miscellaneous_constants=set(())
)


class ConstantLibrary(object):
    """
    Maintain constant values and categories. Currently available
    constants are accessible as :code:`constant_library.current_constants`.
    """

    def __init__(self):

        self.init_constant_library()

    def init_constant_library(self):
        """
        Initialise the constant library from sympl's :code:`default_constants`.
        """

        self.current_constants = deepcopy(default_constants)

        self.constant_categories = deepcopy(default_categories)

        for climt_name in constant_name_map.keys():
            self.current_constants[climt_name] = deepcopy(
                default_constants[constant_name_map[climt_name]])

    def reset_constant_library(self):
        """
        Remove all additions to constant library that have happened
        during runtime.
        """

        self.init_constant_library()

    def _list_constants(self, constant_list):

        for quantity in constant_list:
            print("{}: {} {}".format(quantity,
                                     self.current_constants[quantity].values,
                                     self.current_constants[quantity].units))
            print()

    def list_available_constants(self):
        """
        List all constants currently available in the library, along
        with the category to which they belong.

        """

        print('Physical Constants')
        print('===========================')
        self._list_constants(self.constant_categories['physical_constants'])
        print()
        print('Chemical Constants')
        print('===========================')
        self._list_constants(self.constant_categories["chemical_constants"])
        print()
        print('Planetary Constants')
        print('===========================')
        self._list_constants(self.constant_categories["planetary_constants"])
        print()
        print('Stellar Constants')
        print('===========================')
        self._list_constants(self.constant_categories["stellar_constants"])
        print()
        print('Atmospheric Constants')
        print('===========================')
        self._list_constants(self.constant_categories["atmospheric_constants"])
        print()
        print('Condensible Constants')
        print('===========================')
        self._list_constants(self.constant_categories["condensible_constants"])
        print()
        print('Oceanographic Constants')
        print('===========================')
        self._list_constants(self.constant_categories["oceanographic_constants"])
        print()
        print('Miscellaneous Constants')
        print('===========================')
        self._list_constants(self.constant_categories["miscellaneous_constants"])
        print()

    def add_constants_to_library(self, constant_descriptions):
        """
        Add new constants to the library.

        Args:

           constant_descriptions (dict):
               Dictionary containing the description of the constants.
               The key should be the name of the constant, and the value
               should be a dictionary containing the following keys:

               * type (string):
                   Type of constant. Is one of:

                       * :code:`'physical_constants'`,
                       * :code:`'chemical_constants'`,
                       * :code:`'planetary_constants'`,
                       * :code:`'atmospheric_constants'`,
                       * :code:`'oceanographic_constants'`,
                       * :code:`'stellar_constants'`,
                       * :code:`'condensible_constants'`,
                       * :code:`'miscellaneous_constants'`,
                       * value (float):
                           The value assigned.
               * units (string):
                   The units of the value, e.g, m/s, J/kg.

        Raises:

           IndexError:
               If the constant exists in the library or if the type
               of constant is not valid.
        """

        for name in constant_descriptions.keys():

            if name in self.current_constants.keys():
                raise IndexError('{} already present in library.'.format(name))
            value = constant_descriptions[name]['value']
            units = constant_descriptions[name]['units']
            type = constant_descriptions[name]['type']

            if type not in self.constant_categories.keys():
                raise IndexError('{} is not a valid category.'.format(type))

            self.constant_categories[type] = self.constant_categories[type].union((name,))

            self.current_constants[name] = DataArray(value, attrs={'units': units})

    def modify_constants_in_library(self, constant_descriptions):
        """
        Modify existing constants in the library.

        Args:

           constant_descriptions (dict):
               Dictionary containing the description of the constants.
               The key should be the name of the constant, and the value
               should be a dictionary containing the following keys:

               * value (float):
                   The value assigned.
               * units (string):
                   The units of the value, e.g, m/s, J/kg.

        Raises:

            IndexError: If the constant does not exist in the library.

        """

        for name in constant_descriptions.keys():

            if name not in self.current_constants.keys():
                raise IndexError('{} does not exist in library.'.format(
                    name))
            value = constant_descriptions[name]['value']
            units = constant_descriptions[name]['units']

            self.current_constants[name] = DataArray(value, attrs={'units': units})


constant_library = ConstantLibrary()


def determine_constant_type(name):
    """
    Determine constant type using constant_categories.

    Args:

        name (string):
            the name of the constant.

    Returns:
        constant_type (string):
            the type of the constant.
    """

    for constant_type in constant_library.constant_categories.keys():
        if name in constant_library.constant_categories[constant_type]:
            return constant_type

    raise IndexError("{} not present in any category".format(name))


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

    if name not in constant_library.current_constants.keys():
        raise IndexError("{} not present in the library!".format(name))

    return constant_library.current_constants[name]


def remove_constant(name):
    """
    Delete constant from the library.

        Args:

            name(string):
                The name of the constant to be deleted.

        Raises:

            IndexError: If the constant does not exist in the library.
    """

    if name not in constant_library.current_constants.keys():
        raise IndexError("{} not present in the library!".format(name))

    constant_library.current_constants.pop(name)
    category = determine_constant_type(name)
    constant_library.constant_categories[category].remove(name)

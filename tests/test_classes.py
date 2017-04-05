from sympl import Prognostic
import numpy as np


class MockPrognostic(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
    }

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithExtraDimensions(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
    }

    extra_dimensions = {'some_other_dimension': np.arange(10)}

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithExtraDimensionsAndSigmaLevels(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
        'sigma_levels': ''
    }

    extra_dimensions = {'some_other_dimension': np.arange(10)}

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithExtraQuantities(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
        'some_quantity': 'dimensionless',
        'sigma_on_interface_levels': 'dimensionless'
    }

    quantity_descriptions = {
        'some_quantity': {
            'dims': ['x', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'init_value': 1.
        }
    }

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithExtraQuantitiesNotDefined(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
        'some_quantity': 'dimensionless',
        'sigma_on_interface_levels': 'dimensionless'
    }

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithMalformedExtraQuantities(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
        'some_quantity': 'dimensionless',
    }

    quantity_descriptions = {
        'some_quantity': {
            'dims': ['x', 'y', 'mid_levels'],
            'units': 'dimensionless',
        }
    }

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithAllAttributes(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
        'some_quantity': 'dimensionless',
    }

    tendencies = {
        'air_temperature': 'degK s^-1',
        'some_quantity': 's^-1'
    }

    diagnostics = {
        'air_pressure': 'Pa',
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': 'kg m^-2'
    }

    outputs = {
        'air_temperature': 'degK',
        'mole_fraction_of_oxygen_in_air': 'millimole/mole',
    }

    quantity_descriptions = {
        'some_quantity': {
            'dims': ['x', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'init_value': 1.
        }
    }

    def __init__(self):
        return

    def __call__(self, state):
        return

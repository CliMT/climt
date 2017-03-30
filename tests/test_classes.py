from sympl import Prognostic
import numpy as np


class MockPrognostic(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'oxygen_mixing_ratio': 'millimole/mole',
    }

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithExtraDimensions(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'oxygen_mixing_ratio': 'millimole/mole',
    }

    extra_dimensions = {'some_other_dimension': np.arange(10)}

    def __init__(self):
        return

    def __call__(self, state):
        return


class MockPrognosticWithExtraQuantities(Prognostic):

    inputs = {
        'air_temperature': 'degK',
        'oxygen_mixing_ratio': 'millimole/mole',
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

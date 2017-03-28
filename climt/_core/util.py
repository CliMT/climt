from pint import UnitRegistry
from sympl import DataArray, get_numpy_array
from .initialization import _quantity_descriptions
import numpy as np

pint_units = UnitRegistry()

def mass_to_volume_mixing_ratio(
    mass_mixing_ratio,
    molecular_weight=None,
    molecular_weight_air=28.964):
    """
    Converts from mass mixing ratio (mass per unit mass) to volume
    mixing ratio (volume per unit volume)

    Args:

    mass_mixing_ratio (DataArray): The quantity to be transformed. It
        must be a DataArray with the units attribute correctly set (i.e, to
        something like 'g/kg' or 'kg/kg').

    molecular_weight (float): The molecular weight of the gas in grams/mole.

    molecular_weight_air (float,optional): The molecular weight of dry air.
        If it is not provided, the value for dry air on earth (28.964 g/mol)
        is used.

    Returns:

    volume_mixing_ratio (DataArray): The volume mixing ratio of the gas.
    """

    if molecular_weight is None:
        raise ValueError('The molecular weight must be provided')

    mass_mixing_ratio_with_units = mass_mixing_ratio.values*pint_units(mass_mixing_ratio.units)

    dims = mass_mixing_ratio.dims

    volume_mixing_ratio = mass_mixing_ratio_with_units*molecular_weight_air/molecular_weight

    volume_mixing_ratio = DataArray(volume_mixing_ratio.to_base_units(), 
                                    dims=dims, attrs={'units': str(volume_mixing_ratio.units)})

    return volume_mixing_ratio

def get_input_arrays_from_state(component, state, memory_layout='fortran'):
    """
    Helper function to get a list of numpy arrays from the state dictionary.
    Returns arrays with dimensions (x,y,z) in the same order as specified in
    ``component`` (first priority) or ``_quantity_descriptions`` in
    ``initialization.py`` (second priority).

    Args:

    component (Prognostic, Diagnostic, Implicit, TimeStepper): The component for
        which the input arrays are required.

    state (dict): The state dictionary.

    memory_layout (string): String which is either 'fortran' or 'c'. This specifies
        the memory layout which the component expects the arrays in. If the arrays
        in ``state`` are not in this memory layout and/or not memory aligned,
        a copy will be made.

    Returns:

    array_dict (dict): dictionary whose values are numpy arrays corresponding
        to the input quantities specified in ``component``. The returned arrays will
        be in the units specified in ``component.inputs``.
    """

    if not hasattr(component, 'inputs'):
        raise IndexError(
            'The component does not have the "inputs" attribute. \
            This is quite unusual! Contact the developers.')

    if not isinstance(component.inputs, dict):
        raise NotImplementedError(
            'This method will only work with components with dict-like "inputs" attribute')

    if memory_layout not in ['fortran', 'c']:
        raise ValueError(
            'memory_layout can be either fortran or c')

    array_dict = {}

    for input_quantity in component.inputs.keys():

        dims = get_dimensions_for(component, input_quantity)
        units = component.inputs[input_quantity]

        input_array = get_array_from_state(state,
                                           input_quantity,
                                           units,
                                           dims)

        if memory_layout is 'fortran' and not input_array.flags['FARRAY']:
            input_array = np.asfortranarray(input_array)
        elif memory_layout is 'c' and not input_array.flags['CARRAY']:
            input_array = np.ascontiguousarray(input_array)

        array_dict[input_quantity] = input_array

    return array_dict

def get_dimensions_for(component, quantity_name):

    if hasattr(component, 'quantity_descriptions'):
        if quantity_name in component.quantity_descriptions:
            return component.quantity_descriptions[quantity_name]['dims']

    if quantity_name in _quantity_descriptions:
        return _quantity_descriptions[quantity_name]['dims']

    #Should never come here.
    raise IndexError(
        '{} not described either by the component or by CliMT!'.format(quantity_name))

def get_array_from_state(state,
                         quantity_name,
                         quantity_units,
                         quantity_dims):

    if quantity_name not in state:
        raise IndexError(
            'The input state does not contain {}'.format(quantity_name))

    return get_numpy_array(
        state[quantity_name].to_units(quantity_units), quantity_dims)

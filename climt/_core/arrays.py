from pint import UnitRegistry
from sympl import DataArray, get_numpy_array, jit
from .initialization import _quantity_descriptions
import numpy as np


def get_numpy_arrays_from_state(component, attribute, state, memory_layout='fortran'):
    """

    Extract required inputs as numpy arrays from state.

    Returns arrays with dimensions (x,y,z) in the same order as specified in
    :code:`component` (first priority) or :code:`_quantity_descriptions` in
    :code:`initialization.py` (second priority).

    Args:

        component (Prognostic, Diagnostic, Implicit, TimeStepper):
            The component for which the input arrays are required.
            
        attribute (basestring):
            The attribute (:code:`inputs`, :code:`tendencies`, :code:`outputs` or
            :code:`diagnostics`) for which to return numpy arrays.

        state (dict):
            The state dictionary.

        memory_layout (string, optional):
            String which is either :code:`'fortran'` or :code:`'c'`. This specifies
            the memory layout which the component expects the arrays in. If the arrays
            in :code:`state` are not in this memory layout and/or not memory aligned,
            a copy will be made.

    Returns:

        array_dict (dict):
            dictionary whose values are numpy arrays corresponding
            to the input quantities specified in :code:`component`. The returned arrays will
            be in the units specified in :code:`attribute`.

    Raises:

        NotImplementedError:
            If the component's :code:`inputs` attribute is not a dictionary.

        ValueError:
            If the :code:`memory_layout` argument is neither "fortran" nor "c".

    """


    if not hasattr(component, attribute):
        raise IndexError(
            'Component has no attribute called {}'.format(attribute))

    quantities_to_extract = getattr(component, attribute)

    if not isinstance(quantities_to_extract, dict):
        raise NotImplementedError(
            'This method will only work with components with dict-like "inputs" attribute')

    if memory_layout not in ['fortran', 'c']:
        raise ValueError(
            'memory_layout can be either fortran or c')

    array_dict = {}

    for quantity in quantities_to_extract.keys():

        dims = get_dimensions_for(component, quantity)
        units = quantities_to_extract[quantity]

        new_array = get_array_from_state(state,
                                           quantity,
                                           units,
                                           dims)

        if memory_layout is 'fortran' and not new_array.flags['FARRAY']:
            new_array = np.asfortranarray(new_array)
        elif memory_layout is 'c' and not new_array.flags['CARRAY']:
            new_array = np.ascontiguousarray(new_array)

        array_dict[quantity] = new_array

    return array_dict


def get_dimensions_for(component, quantity_name):

    if hasattr(component, 'quantity_descriptions'):
        if quantity_name in component.quantity_descriptions:
            return component.quantity_descriptions[quantity_name]['dims']

    if quantity_name in _quantity_descriptions:
        return _quantity_descriptions[quantity_name]['dims']

    # Should never come here.
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

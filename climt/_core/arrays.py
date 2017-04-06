from sympl import get_numpy_array
from .initialization import quantity_descriptions, get_default_values
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

        attribute (string):
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

    quantities_to_extract = check_if_sane_and_return_attribute(component, attribute)

    if memory_layout not in ['fortran', 'c']:
        raise ValueError(
            'memory_layout can be either fortran or c')

    array_dict = {}

    for quantity in quantities_to_extract.keys():

        dims = get_dimensions_for(component, quantity)
        units = quantities_to_extract[quantity]

        new_array = get_array_from_state(
            state,
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

    if quantity_name in quantity_descriptions:
        return quantity_descriptions[quantity_name]['dims']

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


def create_state_dict_for(component, attribute, state):
    """
    Create dictionaries to return to caller.

    Use quantities in :code:`component.attribute` to create a
    dictionary of DataArrays which is returned by the component
    to the caller.

    Args:

        component (Prognostic, Implicit, Diagnostic, ImplicitPrognostic, TimeStepper):
            component for which the output dictionary is required.

        attribute (basestring):
            The attribute of the component which should be used to create the
            dictionary. Typically, one of :code:`inputs`, :code:`tendencies`, :code:`outputs` or
            :code:`diagnostics`.

        state (dict):
            The state dictionary that was passed in to the component

    Returns:

        output_dict (dict):
            The dictionary whose keys are labels from :code:`component.attribute` and
            values are the appropriate DataArrays.

    """

    quantities_to_extract = check_if_sane_and_return_attribute(component, attribute)
    # quantities_to_extract is a dictionary whose keys are quantity names
    # and values are units produced by the code. Note that if this is a
    # tendency term, the final units returned to the caller must be in per second,
    # since the TimeStepper requires quantities in per seconds.

    output_state = {}
    for quantity in quantities_to_extract.keys():
        description = quantity_descriptions.copy()

        if hasattr(component, 'quantity_descriptions'):
            if quantity in component.quantity_descriptions:
                description[quantity] = component.quantity_descriptions[quantity]

        additional_dimensions = {}
        for dimension in description[quantity]['dims']:
            if dimension not in ['x', 'y', 'mid_levels', 'interface_levels']:
                additional_dimensions[dimension] = state[dimension]

        # Set the units according to component's description, not global
        # description
        description[quantity]['units'] = quantities_to_extract[quantity]

        using_2d_coordinates = False
        x_coord = state['x']
        y_coord = state['y']
        z_coord = state['z']

        if x_coord.ndim == 2:
            assert y_coord.ndim == 2
            using_2d_coordinates = True
            x_coord = state['logical_x_coordinate']
            y_coord = state['logical_y_coordinate']

        quantity_data_array = get_default_values(quantity,
                                                 x_coord, y_coord, z_coord,
                                                 description,
                                                 additional_dimensions)

        if using_2d_coordinates:
            physical_x = state['x']
            physical_y = state['y']

            quantity_data_array.coords[physical_x.label] = (
                physical_x.dims, physical_x.values)

            quantity_data_array.coords[physical_y.label] = (
                physical_y.dims, physical_y.values)

        output_state[quantity] = quantity_data_array

    return output_state


def check_if_sane_and_return_attribute(component, attribute):
    """ Check if attribute exists and is a dict"""

    if not hasattr(component, attribute):
        raise IndexError(
            'Component has no attribute called {}'.format(attribute))

    quantities_to_extract = getattr(component, attribute)

    if not isinstance(quantities_to_extract, dict):
        raise NotImplementedError(
            'This method will only work with components with a dict-like {} attribute'.format(attribute))

    return quantities_to_extract

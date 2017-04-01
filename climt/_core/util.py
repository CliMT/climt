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
    g/g or g/kg to mole/mole

    Converts from mass mixing ratio (mass per unit mass) to volume
    mixing ratio (volume per unit volume)

    Args:

        mass_mixing_ratio (DataArray):
            The quantity to be transformed. It
            must be a DataArray with the units attribute correctly set (i.e, to
            something like 'g/kg' or 'kg/kg').

        molecular_weight (float):
            The molecular weight of the gas in grams/mole.

        molecular_weight_air (float,optional):
            The molecular weight of dry air.
            If it is not provided, the value for dry air on earth (28.964 g/mol)
            is used.

    Returns:

        volume_mixing_ratio (DataArray):
            The volume mixing ratio of the gas.

    Raises:

        ValueError:
            if the molecular weight is not provided.

    """

    if molecular_weight is None:
        raise ValueError('The molecular weight must be provided')

    mass_mixing_ratio_with_units = mass_mixing_ratio.values*pint_units(mass_mixing_ratio.units)

    dims = mass_mixing_ratio.dims

    volume_mixing_ratio = mass_mixing_ratio_with_units*molecular_weight_air/molecular_weight

    volume_mixing_ratio = volume_mixing_ratio.to_base_units()
    volume_mixing_ratio = DataArray(volume_mixing_ratio,
                                    dims=dims,
                                    attrs={'units': str(volume_mixing_ratio.units)})

    return volume_mixing_ratio


def get_input_arrays_from_state(component, state, memory_layout='fortran'):
    """

    Extract required inputs as numpy arrays from state.

    Returns arrays with dimensions (x,y,z) in the same order as specified in
    :code:`component` (first priority) or :code:`_quantity_descriptions` in
    :code:`initialization.py` (second priority).

    Args:

        component (Prognostic, Diagnostic, Implicit, TimeStepper):
            The component for which the input arrays are required.

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
            be in the units specified in :code:`component.inputs`.

    Raises:

        NotImplementedError:
            If the component's :code:`inputs` attribute is not a dictionary.

        ValueError:
            If the :code:`memory_layout` argument is neither "fortran" nor "c".

    """

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


def get_interface_values(
        mid_level_values,
        surface_value,
        mid_level_pressure,
        interface_level_pressure):
    """
    Calculate interface values given mid-level values.

    Given 3D values of a quantity on model mid levels (cell centers) and the 2D surface
    value, return the 3D values of that quantity on model full levels (cell interfaces).
    If the z-dimension of :code:`mid_level_values` is length K, the returned array will have a
    z-dimension of length K+1.

    Routine borrowed from CESM (radiation.F90 in rrtmg folder)

    Args:

        mid_level_values (array):
            The values of the quantity on mid-levels.

        surface_value (array):
            The value of the quantity at the surface. Must be in the
            same units as `mid_level_values`

        mid_level_pressure (array):
            Pressure values on mid-levels. Can be in any units.

        interface_level_pressure (array):
            Pressure values on interface levels. Must be in
            in the same units as `mid_level_pressure`.

    Returns:

        interface_values (array):
            values of the quantity on mid-levels.

    """

    interface_values = np.zeros(
        (mid_level_values.shape[0], mid_level_values.shape[1], mid_level_values.shape[2]+1),
        dtype=np.double, order='F')

    log_mid_p = np.log(mid_level_pressure)
    log_int_p = np.log(interface_level_pressure)

    interp_weight = (log_int_p[:, :, 1:-1] - log_mid_p[:, :, 1::])/(log_mid_p[:, :, :-1] -
                                                                    log_mid_p[:, :, 1::])

    interface_values[:, :, 1:-1] = \
        mid_level_values[:, :, 1::] - interp_weight*(mid_level_values[:, :, 1::] -
                                                     mid_level_values[:, :, 0:-1])

    interface_values[:, :, 0] = surface_value[:, :]
    interface_values[:, :, -1] = mid_level_values[:, :, -1]

    return interface_values

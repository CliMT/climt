from .exceptions import SharedKeyException, InvalidStateException
from .constants import default_constants
try:
    from numba import jit
except ImportError:
    # define a function with the same call signature as jit that does nothing
    def jit(signature_or_function=None, **kwargs):
        if signature_or_function is None:
            return lambda x: x
        else:
            return signature_or_function


vertical_dimension_names = [
    'height', 'z', 'alt', 'pres', 'pressure', 'air_pressure', 'altitude']
for name in tuple(vertical_dimension_names):
    vertical_dimension_names.append(name + '_on_full_levels')
    vertical_dimension_names.append(name + '_on_half_levels')
vertical_dimension_names.extend([
    'lev', 'levels', 'half_levels', 'mid_level', 'mid_levels',
    'interface_level', 'interface_levels'])

x_dimension_names = ('x', 'lon', 'longitude')

y_dimension_names = ('y', 'lat', 'latitude')

horizontal_dimension_names = x_dimension_names + y_dimension_names


def combine_dimensions_in_3d(*args):
    """Returns a tuple of 3 dimensions (x, y, z), with names corresponding to
    dimension names from the DataArray objects given by *args when present.
    """
    dims = [None, None, None]
    for value in args:
        for dim in value.dims:
            for i, dim_list in enumerate(
                    [x_dimension_names, y_dimension_names,
                     vertical_dimension_names]):
                if (dim in dim_list) and (dim not in ('x', 'y', 'z')):
                    if dims[i] is None:
                        dims[i] = dim
                    elif dims[i] != dim:
                        raise ValueError(
                            'Multiple dimensions along same x/y/z coordinate')
    for i, default in enumerate(['x', 'y', 'z']):
        if dims[i] is None:
            dims[i] = default
    return dims


def combine_dimensions_in_2d(*args):
    """Returns a tuple of 2 dimensions from (x, y, z), with names corresponding to
    dimension names from the DataArray objects given by *args when present.
    """
    dims = [None, None]
    for value in args:
        for dim in value.dims:
            for i, dim_list in enumerate(
                    [x_dimension_names, y_dimension_names]):
                if dim in dim_list:
                    if dims[i] is None:
                        dims[i] = dim
                    elif dims[i] != dim:
                        raise ValueError(
                            'Multiple dimensions along same x/y coordinate')
    for i, default in enumerate(['x', 'y']):
        if dims[i] is None:
            dims[i] = default
    return dims


def set_prognostic_update_frequency(prognostic_class, update_timedelta):
    """
    Alters a prognostic class so that when it is called, it only computes its
    output once for every period of length update_timedelta. In between these
    calls, the cached output from the last computation will be returned.

    Note that the *class* itself must be updated, not an *instance* of that
    class.

    Once modified, the class requires that the 'time' quantity is set on
    states it receives, and that it is a datetime or timedelta object.

    Example:
        This how the function should be used on a Prognostic class GrayLongwaveRadiation.

        >>> from climt import GrayLongwaveRadiation
        >>> from datetime import timedelta
        >>> set_prognostic_update_frequency(GrayLongwaveRadiation, timedelta(hours=1))
        >>> prognostic = GrayLongwaveRadiation()

    Args:
        prognostic_class (type): A Prognostic class (not an instance).
        update_timedelta (timedelta): The amount that state['time'] must differ
            from when output was cached before new output is
            computed.

    Returns:
        prognostic (Prognostic): The input Prognostic, altered in-place
    """
    prognostic_class._update_timedelta = update_timedelta
    prognostic_class._last_update_time = None
    original_call = prognostic_class.__call__

    def __call__(self, state):
        if (self._last_update_time is None or
                state['time'] >= self._last_update_time + self._update_timedelta):
            self._cached_output = original_call(self, state)
            self._last_update_time = state['time']
        return self._cached_output

    prognostic_class.__call__ = __call__
    return prognostic_class


def replace_none_with_default(constant_name, value):
    """If value is None, returns the default constant for the constant name.
    Otherwise, returns value."""
    if value is None:
        return default_constants[constant_name]
    else:
        return value


def add_dicts_inplace(dict1, dict2):
    """
    Takes two dictionaries. For any keys in both dictionary, it adds the value
    in dict2 to the value in dict1. This is done in-place if the values are
    array-like, to avoid data copying. None is returned.
    """
    for key in dict2.keys():
        if key not in dict1:
            dict1[key] = dict2[key]
        else:
            try:
                # works for array-like objects, in-place
                dict1[key][:] += dict2[key][:]
            except TypeError:
                dict1[key] += dict2[key]
    return  # not returning anything emphasizes that this is in-place


def ensure_no_shared_keys(dict1, dict2):
    """
    Raises SharedKeyException if there exists a key present in both
    dictionaries.
    """
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    if len(shared_keys) > 0:
        raise SharedKeyException(
            'unexpected shared keys: {}'.format(shared_keys))


def ensure_horizontal_only(value, name):
    """Raises InvalidStateException if the DataArray given by value contains
    any non-horizontal dimensions."""
    if len(value.dims) > 2:
        raise InvalidStateException(
            '{} is not horizontal, it has more than 2 dimensions'.format(name))
    for dim in value.dims:
        if dim not in horizontal_dimension_names:
            raise InvalidStateException(
                '{} has non-horizontal dimension {}'.format(name, dim)
            )


def get_3d_numpy_array(data_array):
    """Takes in a DataArray, and returns a (x, y, z) 3-dimensional numpy
    array from that DataArray."""
    indices = [None, None, None]
    dimensions = []
    for i, dimension_names in zip(
            range(3), [x_dimension_names, y_dimension_names,
                       vertical_dimension_names]):
        dims = set(data_array.dims).intersection(dimension_names)
        if len(dims) == 1:
            dim = dims.pop()
            dimensions.append(dim)
            indices[i] = slice(0, len(data_array.coords[dim]))
        elif len(dims) > 1:
            raise ValueError(
                'DataArray has multiple dimensions for a single direction')
    if len(dimensions) < len(data_array.dims):
        raise ValueError(
            'Was not able to classify all dimensions as x/y/z: {}'.format(
                data_array.dims))
    # Transpose correctly orders existing dimensions, indices creates new ones
    return data_array.transpose(*dimensions).values[indices]


def get_2d_numpy_array(data_array):
    """Takes in a DataArray, and returns a (x, y) 2-dimensional numpy
    array from that DataArray."""
    indices = [None, None]
    dimensions = []
    for i, dimension_names in zip(
            range(2), [x_dimension_names, y_dimension_names]):
        dims = set(data_array.dims).intersection(dimension_names)
        if len(dims) == 1:
            dim = dims.pop()
            dimensions.append(dim)
            indices[i] = slice(0, len(data_array.coords[dim]))
        elif len(dims) > 1:
            raise ValueError(
                'DataArray has multiple dimensions for a single direction')
    if len(dimensions) < len(data_array.dims):
        raise ValueError(
            'Was not able to classify all dimensions as x/y: {}'.format(
                data_array.dims))
    # Transpose correctly orders existing dimensions, indices creates new ones
    return data_array.transpose(*dimensions).values[indices]

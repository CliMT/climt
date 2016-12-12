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
vertical_dimension_names.extend(['lev', 'levels', 'half_levels'])

horizontal_dimension_names = (
    'x', 'y', 'lon', 'lat', 'longitude', 'latitude',
)


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


def ensure_shared_coordinates(**args):
    """Raises InvalidStateException if the DataArrays given do not share the
    same coordinates."""
    reference_dims = args[0].dims
    reference_coords = args[0].coords
    for value in args:
        assert value.dims == reference_dims
        for i in range(len(value.dims)):
            assert (
                value.coords[value.dims[i]] ==
                reference_coords[reference_dims[i]])


def ensure_third_dim_is_vertical(value, name):
    """Raises InvalidStateException if the DataArray given by value does not
    contain a third dimension which is vertical."""
    if len(value.dims) < 3:
        raise InvalidStateException(
            '{} has fewer than 3 dimensions'.format(name))
    elif value.dims[2] not in vertical_dimension_names:
        raise InvalidStateException(
            '{} does not have a vertical coordinate as its '
            'third dimension'.format(name))


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


def ensure_number_of_dims(value, num_dims, name):
    """Raises InvalidStateException if the DataArray given by value does not
    have precisely num_dims dimensions."""
    if len(value.dims) != num_dims:
        raise InvalidStateException(
            '{} must have {} dimensions'.format(name, num_dims))


def ensure_3d(value, data_dim_if_1d):
    """If the given array is 3D, it is returned. If it is 1D, that array is
    returned broadcast to 3D with the original dim becoming dimension
    data_dim_if_1d. data_dim_if_1d should be an integer from 1 to 3."""
    if len(value.shape) == 3:
        return value
    elif len(value.shape) == 1:
        if data_dim_if_1d == 1:
            return value[:, None, None]
        elif data_dim_if_1d == 2:
            return value[None, :, None]
        elif data_dim_if_1d == 3:
            return value[None, None, :]
        else:
            raise ValueError(
                'data_dim_if_1d should be an integer from 1 to 3, '
                'got {}'.format(data_dim_if_1d))
    else:
        raise ValueError(
            'value should be a 1D or 3D array, instead got {}D'.format(
                len(value.shape)))

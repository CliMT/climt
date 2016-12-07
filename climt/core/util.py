from .exceptions import SharedKeyException, InvalidStateException

vertical_dimension_names = [
    'height', 'z', 'alt', 'pres', 'pressure', 'air_pressure', 'altitude']
for name in tuple(vertical_dimension_names):
    vertical_dimension_names.append(name + '_on_full_levels')
    vertical_dimension_names.append(name + '_on_half_levels')
vertical_dimension_names.extend(['lev', 'levels', 'half_levels'])

horizontal_dimension_names = (
    'x', 'y', 'lon', 'lat', 'longitude', 'latitude',
)


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

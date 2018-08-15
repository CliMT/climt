import numpy as np
from sympl import DataArray, InvalidPropertyDictError
import xarray as xr


class IncompatibleDimensionsError(Exception):
    pass


def combine_dims(dims1, dims2):
    """
    Takes in two dims specifications and returns a single specification that
    satisfies both, if possible. Raises an InvalidDimensionsError if not.

    Parameters
    ----------
    dims1 : iterable of str
    dims2 : iterable of str

    Returns
    -------
    dims : iterable of str

    Raises
    ------
    IncompatibleDimensionsError
        If the two dims specifications cannot be combined
    """
    if '*' in dims2 and '*' not in dims1:
        dims1, dims2 = dims2, dims1
    if '*' in dims1:
        i_wildcard = dims1.index('*')
        dims_out = list(dims1)
        if '*' not in dims2:
            # dims2 is a hard definition, so we remove the wildcard
            dims_out.pop(i_wildcard)
            dims2_unmatched_dims = set(dims1).difference(dims2).difference(['*'])
            if len(dims2_unmatched_dims) != 0:
                raise IncompatibleDimensionsError(
                    'Dimensions {} and {} are incompatible'.format(dims1, dims2))
        dims1_unmatched_dims = set(dims2).difference(dims1)
        for dim in reversed(dims2):  # reversed cancels reversal from insert
            if dim in dims1_unmatched_dims:
                dims_out.insert(i_wildcard, dim)
        return tuple(dims_out)
    elif set(dims1) != set(dims2):
        raise IncompatibleDimensionsError(
            'Dimensions {} and {} are incompatible'.format(dims1, dims2))
    else:  # dims1 and dims2 contain the same dimensions
        dims_out = tuple(dims1)
    return dims_out


def combine_dims_dict_list(dims_dict_list):
    """

    Args:
        dims_dict_list : list of dict
            Dictionaries of name: dims pairs to be combined.

    Returns:
        dims_dict : dict
            name: dims pairs that are compatible with all pairs in the
            combined lists, if possible.

    Raises:
        IncompatibleDimensionsError
            If the dimensions dictionaries cannot be combined due to
            incompatible dimensions between the dictionaries.
    """
    return_dict = {}
    for dims_dict in dims_dict_list:
        for name, dims in dims_dict.items():
            if name not in return_dict:
                return_dict[name] = dims
            else:
                return_dict[name] = combine_dims(return_dict[name], dims)
    return return_dict


def get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard):
    for name, dims in dims_dict_without_wildcard.items():
        if '*' in dims:
            raise ValueError(
                'dims_dict_without_wildcard should have no wildcards, but '
                'wildcard is present for {}'.format(name)
            )
    cannot_be_wildcard_dims = set()
    wildcard_dims = []
    for name, dims_with_wildcard in dims_dict_with_wildcard.items():
        cannot_be_wildcard_dims.update(dims_with_wildcard)
        if name in dims_dict_without_wildcard:
            for dim in dims_dict_without_wildcard[name]:
                if dim not in dims_with_wildcard and dim not in wildcard_dims:
                    wildcard_dims.append(dim)
    bad_wildcard_dims = cannot_be_wildcard_dims.intersection(wildcard_dims)
    if len(bad_wildcard_dims) > 0:
        raise IncompatibleDimensionsError('')
    return tuple(wildcard_dims)


def get_dims_dict_from_state(state):
    return {name: array.dims for name, array in state.items() if name != 'time'}


def get_dims_dict_from_properties(properties, fallback_properties=None):
    return_dict = {}
    for name, value in properties.items():
        if 'dims' in value:
            return_dict[name] = tuple(value['dims'])
        elif fallback_properties is not None:
            try:
                return_dict[name] = fallback_properties[name]['dims']
            except KeyError:
                raise InvalidPropertyDictError(
                    'No dimensions defined for {}'.format(name)
                )
        else:
            raise InvalidPropertyDictError(
                'No dimensions defined for {}'.format(name)
            )
    return return_dict


def expand_wildcard_dim(dims_dict_with_wildcard, wildcard_dim_list):
    return_dict = {}
    for name, dims in dims_dict_with_wildcard.items():
        if '*' in dims:
            out_dims = list(dims)
            i_wildcard = dims.index('*')
            out_dims.pop(i_wildcard)
            for wildcard_dim in reversed(wildcard_dim_list):
                out_dims.insert(i_wildcard, wildcard_dim)
            return_dict[name] = tuple(out_dims)
        else:
            return_dict[name] = tuple(dims)
    return return_dict


def get_dim_lengths(state):
    return_dict = {}
    for name, array in state.items():
        if name != 'time':
            for dim, length in zip(array.dims, array.shape):
                if dim not in return_dict:
                    return_dict[dim] = length
                elif return_dict[dim] != length:
                    raise IncompatibleDimensionsError(
                        'dimension {} has conflicting lengths {} and {}'.format(
                            dim, return_dict[dim], length
                        )
                    )
    return return_dict


def get_initial_dims(component_list, grid_state):
    last_dims_dict = {}
    dims_dict = get_dims_dict_from_state(grid_state)
    while dims_dict != last_dims_dict:
        last_dims_dict = dims_dict
        dims_dict_with_wildcard = combine_dims_dict_list(*[get_dims_dict_from_properties(c.input_properties) for c in component_list])
        wildcard_dims = get_wildcard_dims(
            dims_dict_with_wildcard,
            last_dims_dict,
        )
        dims_dict = expand_wildcard_dim(dims_dict_with_wildcard, wildcard_dims)
    return dims_dict


def broadcast_state(state, dims_dict, dim_lengths):
    return_state = {}
    for name, array in state.items():
        if name in dims_dict:
            return_state[name] = transform_array_to_dims(
                array, dims_dict[name], dim_lengths)
        else:
            return_state[name] = state[name]
    return return_state


def transform_array_to_dims(array, out_dims, dim_lengths):
    for dim, length in zip(array.dims, array.shape):
        if dim_lengths[dim] != length:
            raise IncompatibleDimensionsError(
                'array has length {} for dimension {} but dim_lengths has {}'.format(
                    length, dim, dim_lengths[dim]
                )
            )
    missing_dims = set(array.dims).difference(out_dims)
    if len(missing_dims) != 0:
        raise IncompatibleDimensionsError(
            'array has dimensions {} not included in out_dims'.format(missing_dims)
        )
    if tuple(out_dims) == array.dims:
        return_array = array
    else:
        new_shape = tuple(dim_lengths[dim_name] for dim_name in out_dims)
        return_array = DataArray(
            np.empty(new_shape, dtype=array.dtype),
            dims=out_dims,
            attrs=array.attrs
        )
        _, broadcast_array = xr.broadcast(return_array, array)
        return_array[:] = broadcast_array
    return return_array


# def convert_state_units(state, units_dict):
#     pass
#
#
# def get_units_dict_from_properties(properties):
#     pass
#
#
# def get_arrays_from_state(state):
#     pass
#
#
# def insert_wildcard_dims(dims, wildcard_dims):
#     pass

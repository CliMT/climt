import unittest
from climt._core.properties import (
    combine_dims_dict_list, get_wildcard_dims, get_dims_dict_from_state,
    get_dims_dict_from_properties, expand_wildcard_dim, get_dim_lengths,
    # get_initial_dims, convert_state_units, get_units_dict_from_properties,
    broadcast_state, transform_array_to_dims,  # get_arrays_from_state,
    IncompatibleDimensionsError
)
from sympl import DataArray
import numpy as np
import pytest


class CombineDimsDictListTests(unittest.TestCase):

    def test_combine_one_empty_dict(self):
        return_value = combine_dims_dict_list([{}])
        assert return_value == {}

    def test_combine_many_empty_dict(self):
        return_value = combine_dims_dict_list([{}, {}, {}])
        assert return_value == {}

    def test_combine_one_nonempty_dict(self):
        dims_dict = {'name1': ('dim1', 'dim2'), 'name2': ('*',)}
        return_value = combine_dims_dict_list([dims_dict])
        assert return_value == dims_dict

    def test_combine_two_dicts_one_nonempty(self):
        dims_dict = {'name1': ('dim1', 'dim2'), 'name2': ('*',)}
        return_value = combine_dims_dict_list([dims_dict, {}])
        assert return_value == dims_dict
        return_value = combine_dims_dict_list([{}, dims_dict])
        assert return_value == dims_dict

    def test_combine_two_compatible_dicts(self):
        dims_dict_1 = {
            'name1': ('mid_levels', '*'),
        }
        dims_dict_2 = {
            'name1': ('mid_levels', 'latitude', 'longitude'),
            'name2': ('mid_levels', 'latitude', 'longitude'),
        }
        return_value = combine_dims_dict_list([dims_dict_1, dims_dict_2])
        assert return_value == dims_dict_2
        return_value = combine_dims_dict_list([dims_dict_2, dims_dict_1])
        assert return_value == dims_dict_2

    def test_combine_two_compatible_dicts_both_wildcard(self):
        dims_dict_1 = {
            'name1': ('mid_levels', '*'),
        }
        dims_dict_2 = {
            'name1': ('mid_levels', 'latitude', 'longitude', '*'),
            'name2': ('mid_levels', 'latitude', 'longitude'),
        }
        return_value = combine_dims_dict_list([dims_dict_1, dims_dict_2])
        assert return_value == dims_dict_2
        return_value = combine_dims_dict_list([dims_dict_2, dims_dict_1])
        assert return_value == dims_dict_2

    def test_combine_three_compatible_dicts(self):
        dims_dict_1 = {
            'name1': ('mid_levels', '*'),
        }
        dims_dict_2 = {
            'name1': ('mid_levels', 'latitude', '*'),
        }
        dims_dict_3 = {
            'name1': ('*', 'longitude', 'mid_levels'),
        }
        return_value = combine_dims_dict_list(
            [dims_dict_1, dims_dict_2, dims_dict_3])
        assert len(return_value) == 1
        assert set(return_value['name1']) == {'mid_levels', 'latitude', '*', 'longitude'}

    def test_combine_two_incompatible_dicts(self):
        dims_dict_1 = {
            'name1': ('mid_levels', 'latitude'),
        }
        dims_dict_2 = {
            'name1': ('mid_levels', 'longitude'),
        }
        with pytest.raises(IncompatibleDimensionsError):
            combine_dims_dict_list([dims_dict_1, dims_dict_2])

    def test_combine_incompatible_dicts_with_wildcard(self):
        dims_dict_1 = {
            'name1': ('mid_levels', 'latitude', '*'),
        }
        dims_dict_2 = {
            'name1': ('mid_levels', 'longitude'),
        }
        with pytest.raises(IncompatibleDimensionsError):
            combine_dims_dict_list([dims_dict_1, dims_dict_2])


class GetWildcardDimsTests(unittest.TestCase):

    def test_empty_dicts(self):
        result = get_wildcard_dims({}, {})
        assert len(result) == 0

    def test_empty_dims_dict_with_wildcard(self):
        dims_dict_with_wildcard = {}
        dims_dict_without_wildcard = {
            'name1': ['dim1', 'dim2']
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert len(result) == 0

    def test_empty_dims_dict_without_wildcard(self):
        dims_dict_without_wildcard = {}
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*']
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert len(result) == 0

    def test_no_overlapping_names(self):
        dims_dict_without_wildcard = {
            'name2': ['dim1', 'dim2', 'dim3', 'dim4'],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*']
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert len(result) == 0

    def test_invalid_wildcard(self):
        dims_dict_without_wildcard = {
            'name2': ['dim1', 'dim2', 'dim3', '*'],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*']
        }
        with pytest.raises(ValueError):
            get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)

    def test_incompatible_wildcards(self):
        dims_dict_without_wildcard = {
            'name1': ['dim1', 'dim2'],
            'name2': ['dim1', 'dim2']
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*'],
            'name2': ['dim1', '*'],
        }
        with pytest.raises(IncompatibleDimensionsError):
            get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)

    def test_one_overlapping_name(self):
        dims_dict_without_wildcard = {
            'name1': ['dim1', 'dim2', 'dim3', 'dim4'],
            'name2': ['dim5'],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*']
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert len(result) == 2
        assert result == ('dim3', 'dim4')

    def test_two_overlapping_names_same_wildcard_dims(self):
        dims_dict_without_wildcard = {
            'name1': ['dim1', 'dim2', 'dim3', 'dim4'],
            'name2': ['dim4', 'dim3', 'dim5', 'dim6'],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*'],
            'name2': ['*', 'dim5', 'dim6'],
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert len(result) == 2
        assert set(result) == {'dim3', 'dim4'}

    def test_two_overlapping_names_different_wildcard_dims(self):
        dims_dict_without_wildcard = {
            'name1': ['dim1', 'dim2', 'dim3'],
            'name2': ['dim4', 'dim5', 'dim6'],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', 'dim2', '*'],
            'name2': ['*', 'dim5', 'dim6'],
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert len(result) == 2
        assert set(result) == {'dim3', 'dim4'}

    def test_extra_dims_in_wildcard_def_no_wildcards(self):
        dims_dict_without_wildcard = {
            'name1': [],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', '*']
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert result == ()

    def test_extra_dims_in_wildcard_def_with_wildcard(self):
        dims_dict_without_wildcard = {
            'name1': ['dim2'],
        }
        dims_dict_with_wildcard = {
            'name1': ['dim1', '*']
        }
        result = get_wildcard_dims(dims_dict_with_wildcard, dims_dict_without_wildcard)
        assert result == ('dim2',)


class TestGetDimsDictFromState(unittest.TestCase):

    def test_empty(self):
        result = get_dims_dict_from_state({})
        assert len(result) == 0
        assert isinstance(result, dict)

    def test_single(self):
        state = {
            'name1': DataArray(
                np.empty([3, 4]),
                dims=['dim1', 'dim2'],
                attrs={'units': 'm'},
            ),
        }
        result = get_dims_dict_from_state(state)
        assert result == {'name1': ('dim1', 'dim2')}

    def test_multiple(self):
        state = {
            'name1': DataArray(
                np.empty([3, 4]),
                dims=['dim1', 'dim2'],
                attrs={'units': 'm'},
            ),
            'name2': DataArray(
                np.empty([3]),
                dims=['dim1'],
                attrs={'units': 'degK'},
            ),
            'name3': DataArray(
                np.empty([5]),
                dims=['dim3'],
                attrs={'units': 'km'}
            )
        }
        result = get_dims_dict_from_state(state)
        assert result == {
            'name1': ('dim1', 'dim2'),
            'name2': ('dim1',),
            'name3': ('dim3',)
        }


class TestGetDimsDictFromProperties(unittest.TestCase):
    def test_empty(self):
        result = get_dims_dict_from_properties({})
        assert len(result) == 0
        assert isinstance(result, dict)

    def test_single(self):
        properties = {
            'name1': {
                'dims': ['dim1', 'dim2'],
                'units': 'm',
            }
        }
        result = get_dims_dict_from_properties(properties)
        assert result == {'name1': ('dim1', 'dim2')}

    def test_multiple(self):
        properties = {
            'name1': {
                'dims': ['dim1', 'dim2'],
                'units': 'm',
            },
            'name2': {
                'dims': ('dim1',),
                'units': 'degK',
            },
            'name3': {
                'dims': ['dim3'],
                'units': 'km',
            },
        }
        result = get_dims_dict_from_properties(properties)
        assert result == {
            'name1': ('dim1', 'dim2'),
            'name2': ('dim1',),
            'name3': ('dim3',)
        }


class ExpandWildcardDimTests(unittest.TestCase):

    def test_empty(self):
        result = expand_wildcard_dim({}, ['dim1'])
        assert len(result) == 0
        assert isinstance(result, dict)

    def test_single_name_single_wildcard_dim(self):
        result = expand_wildcard_dim(
            {'name1': ['*', 'dim2', 'dim3']},
            ['dim1'],
        )
        assert result == {'name1': ('dim1', 'dim2', 'dim3')}

    def test_single_name_multiple_wildcard_dims(self):
        result = expand_wildcard_dim(
            {'name1': ['*', 'dim2', 'dim3']},
            ['dim0', 'dim1'],
        )
        assert result == {'name1': ('dim0', 'dim1', 'dim2', 'dim3')}
        result = expand_wildcard_dim(
            {'name1': ['dim2', '*', 'dim3']},
            ['dim0', 'dim1'],
        )
        assert result == {'name1': ('dim2', 'dim0', 'dim1', 'dim3')}

    def test_multiple_name_single_wildcard_dim(self):
        result = expand_wildcard_dim(
            {
                'name1': ['*', 'dim2', 'dim3'],
                'name2': ['dim4', '*'],
            },
            ['dim1'],
        )
        assert result == {
            'name1': ('dim1', 'dim2', 'dim3'),
            'name2': ('dim4', 'dim1'),
        }


class GetDimLengthsTests(unittest.TestCase):

    def test_empty(self):
        result = get_dim_lengths({})
        assert len(result) == 0
        assert isinstance(result, dict)

    def test_single(self):
        result = get_dim_lengths({
            'name1': DataArray(
                np.zeros([3, 4]),
                dims=['dim1', 'dim2'],
                attrs={'units': 'm'},
            ),
        })
        assert result == {'dim1': 3, 'dim2': 4}

    def test_multiple_compatible(self):
        result = get_dim_lengths({
            'name1': DataArray(
                np.zeros([3, 4]),
                dims=['dim1', 'dim2'],
                attrs={'units': 'm'},
            ),
            'name2': DataArray(
                np.zeros([4]),
                dims=['dim2'],
                attrs={'units': 'degK'},
            ),
            'name3': DataArray(
                np.zeros([2, 5]),
                dims=['dim3', 'dim4'],
                attrs={'units': 'km'}
            ),
        })
        assert result == {'dim1': 3, 'dim2': 4, 'dim3': 2, 'dim4': 5}

    def test_multiple_incompatible(self):
        with pytest.raises(IncompatibleDimensionsError):
            get_dim_lengths({
                'name1': DataArray(
                    np.zeros([3, 4]),
                    dims=['dim1', 'dim2'],
                    attrs={'units': 'm'},
                ),
                'name2': DataArray(
                    np.zeros([5]),
                    dims=['dim2'],
                    attrs={'units': 'degK'},
                ),
                'name3': DataArray(
                    np.zeros([2, 5]),
                    dims=['dim3', 'dim4'],
                    attrs={'units': 'km'}
                ),
            })


class TransformArrayToDimsTests(unittest.TestCase):

    def test_0d_no_dim_change(self):
        array = DataArray(
            1.,
            dims=[],
            attrs={'units': 'm'}
        )
        out_dims = []
        dim_lengths = {}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result is array

    def test_1d_no_dim_change(self):
        array = DataArray(
            np.zeros([5]),
            dims=['mid_levels'],
            attrs={'units': 'm'}
        )
        out_dims = ['mid_levels']
        dim_lengths = {'mid_levels': 5}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result is array

    def test_1d_invalid_dimension_length(self):
        array = DataArray(
            np.zeros([5]),
            dims=['mid_levels'],
            attrs={'units': 'm'}
        )
        out_dims = ['mid_levels']
        dim_lengths = {'mid_levels': 4}
        with pytest.raises(IncompatibleDimensionsError):
            transform_array_to_dims(array, out_dims, dim_lengths)

    def test_2d_no_dim_change(self):
        array = DataArray(
            np.zeros([5, 3]),
            dims=['mid_levels', 'longitude'],
            attrs={'units': 'm'}
        )
        out_dims = ['mid_levels', 'longitude']
        dim_lengths = {'mid_levels': 5, 'longitude': 3}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result is array

    def test_2d_invalid_dimension_removal(self):
        array = DataArray(
            np.zeros([5, 3]),
            dims=['mid_levels', 'longitude'],
            attrs={'units': 'm'}
        )
        out_dims = ['mid_levels']
        dim_lengths = {'mid_levels': 5, 'longitude': 3}
        with pytest.raises(IncompatibleDimensionsError):
            transform_array_to_dims(array, out_dims, dim_lengths)

    def test_0d_to_2d(self):
        array = DataArray(
            5.,
            dims=[],
            attrs={'units': 'm'}
        )
        out_dims = ['mid_levels', 'longitude']
        dim_lengths = {'mid_levels': 5, 'longitude': 3}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result.shape == (5, 3)
        assert result.dims == ('mid_levels', 'longitude')
        assert np.all(result.values == 5.)

    def test_1d_to_2d(self):
        array = DataArray(
            np.zeros([5]) + 4.,
            dims=['mid_levels'],
            attrs={'units': 'm'}
        )
        out_dims = ['mid_levels', 'longitude']
        dim_lengths = {'mid_levels': 5, 'longitude': 3}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result.shape == (5, 3)
        assert result.dims == ('mid_levels', 'longitude')
        assert np.all(result.values == 4.)

    def test_1d_to_3d(self):
        array = DataArray(
            np.zeros([5]) + 4.,
            dims=['mid_levels'],
            attrs={'units': 'm'}
        )
        out_dims = ['latitude', 'mid_levels', 'longitude']
        dim_lengths = {'latitude': 2, 'mid_levels': 5, 'longitude': 3}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result.shape == (2, 5, 3)
        assert result.dims == ('latitude', 'mid_levels', 'longitude')
        assert np.all(result.values == 4.)

    def test_3d_transpose(self):
        array = DataArray(
            np.zeros([5, 2, 3]) + 4.,
            dims=['mid_levels', 'latitude', 'longitude'],
            attrs={'units': 'm'}
        )
        out_dims = ['latitude', 'longitude', 'mid_levels']
        dim_lengths = {'latitude': 2, 'mid_levels': 5, 'longitude': 3}
        result = transform_array_to_dims(array, out_dims, dim_lengths)
        assert result.shape == (2, 3, 5)
        assert result.dims == ('latitude', 'longitude', 'mid_levels')
        assert np.all(result.values == 4.)


class BroadcastStateTests(unittest.TestCase):

    def test_empty_state(self):
        state = {}
        dims_dict = {}
        dim_lengths = {}
        result = broadcast_state(state, dims_dict, dim_lengths)
        assert len(result) == 0
        assert isinstance(result, dict)

    def test_single_item_state(self):
        state = {
            'name1': DataArray(
                np.zeros([5, 2, 3]) + 4.,
                dims=['mid_levels', 'latitude', 'longitude'],
                attrs={'units': 'm'}
            )
        }
        dims_dict = {'name1': ['latitude', 'longitude', 'mid_levels']}
        dim_lengths = {'latitude': 2, 'mid_levels': 5, 'longitude': 3}
        result = broadcast_state(state, dims_dict, dim_lengths)
        array_result = result['name1']
        assert array_result.shape == (2, 3, 5)
        assert array_result.dims == ('latitude', 'longitude', 'mid_levels')
        assert np.all(array_result.values == 4.)

    def test_single_item_state_from_flat(self):
        state = {
            'name1': DataArray(
                4.,
                dims=[],
                attrs={'units': 'm'}
            )
        }
        dims_dict = {'name1': ['mid_levels']}
        dim_lengths = {'mid_levels': 5}
        result = broadcast_state(state, dims_dict, dim_lengths)
        array_result = result['name1']
        assert array_result.shape == (5,)
        assert array_result.dims == ('mid_levels',)
        assert np.all(array_result.values == 4.)

    def test_multiple_item_state(self):
        state = {
            'name1': DataArray(
                np.zeros([5, 2, 3]) + 4.,
                dims=['mid_levels', 'latitude', 'longitude'],
                attrs={'units': 'm'}
            ),
            'name2': DataArray(
                5.,
                dims=[],
                attrs={'units': 'degK'}
            ),
        }
        dims_dict = {
            'name1': ['latitude', 'longitude', 'mid_levels'],
            'name2': ['mid_levels']
        }
        dim_lengths = {'latitude': 2, 'mid_levels': 5, 'longitude': 3}
        result = broadcast_state(state, dims_dict, dim_lengths)
        array_result = result['name1']
        assert array_result.shape == (2, 3, 5)
        assert array_result.dims == ('latitude', 'longitude', 'mid_levels')
        assert np.all(array_result.values == 4.)
        array_result = result['name2']
        assert array_result.shape == (5,)
        assert array_result.dims == ('mid_levels',)
        assert np.all(array_result.values == 5.)

import pytest
import abc
import os
from glob import glob
import xarray as xr
import numpy as np
from climt import (
    DataArray, HeldSuarez, GrayLongwaveRadiation,
    Frierson06LongwaveOpticalDepth)

cache_folder = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'cached_component_output')


def cache_dictionary(dictionary, filename):
    dataset = xr.Dataset(dictionary)
    dataset.to_netcdf(filename, engine='scipy')


def load_dictionary(filename):
    dataset = xr.open_dataset(filename, engine='scipy')
    return dict(dataset.data_vars)


class ComponentBase(object):

    @abc.abstractmethod
    def get_input_state(self):
        pass

    @abc.abstractmethod
    def get_component_instance(self):
        pass

    def get_cached_output(self):
        cache_filename_list = sorted(glob(
            os.path.join(
                cache_folder,
                '{}-*.cache'.format(
                    self.__class__.__name__))))
        if len(cache_filename_list) > 0:
            return_list = []
            for filename in cache_filename_list:
                return_list.append(load_dictionary(filename))
            if len(return_list) > 1:
                return tuple(return_list)
            elif len(return_list) == 1:
                return return_list[0]
        else:
            return None

    def cache_output(self, output):
        if not isinstance(output, tuple):
            output = (output,)
        for i in range(len(output)):
            cache_filename = os.path.join(
                cache_folder, '{}-{}.cache'.format(self.__class__.__name__, i))
            cache_dictionary(output[i], cache_filename)

    def test_output_matches_cached_output(self):
        state = self.get_input_state()
        component = self.get_component_instance()
        output = component(state)
        cached_output = self.get_cached_output()
        if cached_output is None:
            self.cache_output(output)
            raise AssertionError(
                'Failed due to no cached output, cached current output')
        else:
            compare_outputs(output, cached_output)

    def test_component_listed_inputs_are_accurate(self):
        state = self.get_input_state()
        component = self.get_component_instance()
        input_state = {}
        for key in component.inputs:
            input_state[key] = state[key]
        output = component(state)
        cached_output = self.get_cached_output()
        if cached_output is not None:
            compare_outputs(output, cached_output)

    def test_consistent_dim_length(self):
        input_state = self.get_input_state()
        assert_dimension_lengths_are_consistent(input_state)
        component = self.get_component_instance()
        output = component(input_state)
        if isinstance(output, tuple):
            # Check diagnostics/tendencies/outputs are consistent with one
            # another
            test_state = {}
            for state in output:
                test_state.update(state)
            assert_dimension_lengths_are_consistent(test_state)
        else:
            test_state = output  # if not a tuple assume it's a dict
            assert_dimension_lengths_are_consistent(test_state)


def assert_dimension_lengths_are_consistent(state):
    dimension_lengths = {}
    for name, value in state.items():
        for i, dim_name in enumerate(value.dims):
            if dim_name in dimension_lengths:
                assert dimension_lengths[dim_name] == value.shape[i]
            else:
                dimension_lengths[dim_name] = value.shape[i]


def compare_outputs(current, cached):
    if isinstance(current, tuple) and isinstance(cached, tuple):
        for i in range(len(current)):
            compare_one_state_pair(current[i], cached[i])
    elif (not isinstance(current, tuple)) and (not isinstance(cached, tuple)):
        compare_one_state_pair(current, cached)
    else:
        raise AssertionError('Different number of dicts returned than cached.')


def compare_one_state_pair(current, cached):
    for key in current.keys():
        assert np.all(current[key].values == cached[key].values)
        for attr in current[key].attrs:
            assert current[key].attrs[attr] == cached[key].attrs[attr]
        for attr in cached[key].attrs:
            assert attr in current[key].attrs
        assert current[key].dims == cached[key].dims
    for key in cached.keys():
        assert key in current.keys()


class TestHeldSuarez(ComponentBase):

    def get_input_state(self):
        random = np.random.RandomState(0)
        return {
            'latitude': DataArray(
                random.randn(3), dims=['lat'], attrs={'units': 'degrees_N'}),
            'air_pressure': DataArray(
                random.rand(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'hPa'},),
            'sigma': DataArray(
                np.linspace(0., 1., num=6), dims=['lev'], attrs={'units': ''}),
            'air_temperature': DataArray(
                270. + random.randn(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'degK'}),
            'eastward_wind': DataArray(
                random.randn(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'm/s'}),
            'northward_wind': DataArray(
                random.randn(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'm/s'}),
        }

    def get_component_instance(self):
        return HeldSuarez()


class TestHeldSuarezCachedCoordinates(ComponentBase):

    def get_input_state(self):
        random = np.random.RandomState(0)
        return {
            'latitude': DataArray(
                np.linspace(-90, 90, num=3),
                dims=['lat'], attrs={'units': 'degrees_N'}),
            'air_pressure': DataArray(
                random.rand(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'hPa'},),
            'sigma': DataArray(
                np.linspace(0., 1., num=6), dims=['lev'], attrs={'units': ''}),
            'air_temperature': DataArray(
                270. + random.randn(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'degK'}),
            'eastward_wind': DataArray(
                random.randn(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'm/s'}),
            'northward_wind': DataArray(
                random.randn(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'm/s'}),
        }

    def get_component_instance(self):
        random = np.random.RandomState(0)
        return HeldSuarez(
            latitude=DataArray(
                np.linspace(-90, 90, num=3),
                dims=['lat'], attrs={'units': 'degrees_N'}),
            air_pressure=DataArray(
                random.rand(2, 3, 6), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'hPa'},),
            sigma=DataArray(
                np.linspace(0., 1., num=6), dims=['lev'], attrs={'units': ''}))


class TestFrierson06LongwaveOpticalDepth(ComponentBase):

    def get_input_state(self):
        state = {
            'latitude': DataArray(
                np.linspace(-90, 90, num=10),
                dims=['lat'], attrs={'units': 'degrees_N'}),
            'sigma_on_interface_levels': DataArray(
                np.linspace(0., 1., num=6),
                dims=['interface_levels'], attrs={'units': ''}),
        }
        return state

    def get_component_instance(self):
        return Frierson06LongwaveOpticalDepth()


class TestGrayLongwaveRadiation(ComponentBase):

    def get_component_instance(self):
        return GrayLongwaveRadiation()

    def get_input_state(self):
        random = np.random.RandomState(1)
        nx, ny, nz = 4, 4, 10
        state = {
            'longwave_optical_depth_on_interface_levels': DataArray(
                (np.linspace(0, 6, nz+1)[None, None, :] *
                 (1 + 0.1*random.randn(nx, ny))[:, :, None]),
                dims=['x', 'y', 'interface_levels'], attrs={'units': ''},
            ),
            'air_temperature': DataArray(
                5*random.randn(nx, ny, nz) + 270.,
                dims=['x', 'y', 'mid_levels'], attrs={'units': 'degK'},
            ),
            'air_pressure_on_interface_levels': DataArray(
                np.linspace(1e5, 0, nz+1),
                dims=['interface_levels'], attrs={'units': 'Pa'},
            ),
            'surface_temperature': DataArray(
                5 * random.randn(nx, ny) + 270.,
                dims=['x', 'y'], attrs={'units': 'degK'},
            )}
        return state


class TestGrayLongwaveRadiation1D(ComponentBase):

    def get_component_instance(self):
        return GrayLongwaveRadiation()

    def get_input_state(self):
        random = np.random.RandomState(1)
        nz = 10
        state = {
            'longwave_optical_depth_on_interface_levels': DataArray(
                np.linspace(0, 6, nz+1),
                dims=['interface_levels'], attrs={'units': ''},
            ),
            'air_temperature': DataArray(
                5*random.randn(nz) + 270.,
                dims=['mid_levels'], attrs={'units': 'degK'},
            ),
            'air_pressure_on_interface_levels': DataArray(
                np.linspace(1e5, 0, nz+1),
                dims=['interface_levels'], attrs={'units': 'Pa'},
            ),
            'surface_temperature': DataArray(
                270.,
                dims=[], attrs={'units': 'degK'},
            )}
        return state


if __name__ == '__main__':
    pytest.main([__file__])
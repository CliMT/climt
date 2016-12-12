import pytest
import abc
import os
from glob import glob
import xarray as xr
import numpy as np
from climt import HeldSuarez, DataArray

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
            return return_list
        else:
            return None

    def cache_output(self, output):
        for i in range(len(output)):
            cache_filename = os.path.join(
                cache_folder, '{}-{}.cache'.format(self.__class__.__name__, i))
            cache_dictionary(output[i], cache_filename)

    def compare_outputs(self, current, cached):
        for i in range(len(current)):
            for key in current[i].keys():
                assert np.all(current[i][key].values == cached[i][key].values)
                for attr in current[i][key].attrs:
                    assert current[i][key].attrs[attr] == cached[i][key].attrs[attr]
                for attr in cached[i][key].attrs:
                    assert attr in current[i][key].attrs
                assert current[i][key].dims == cached[i][key].dims
            for key in cached[i].keys():
                assert key in current[i].keys()

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
            self.compare_outputs(output, cached_output)


class TestHeldSuarez(ComponentBase):

    def get_input_state(self):
        random = np.random.RandomState(0)
        return {
            'latitude': DataArray(
                random.randn(10), dims=['lat'], attrs={'units': 'degrees_N'}),
            'air_pressure': DataArray(
                random.rand(5, 10, 20), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'hPa'},),
            'sigma': DataArray(
                np.linspace(0., 1., num=20), dims=['lev'], attrs={'units': ''}),
            'air_temperature': DataArray(
                270. + random.randn(5, 10, 20), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'degK'}),
            'eastward_wind': DataArray(
                random.randn(5, 10, 20), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'm/s'}),
            'northward_wind': DataArray(
                random.randn(5, 10, 20), dims=['lon', 'lat', 'lev'],
                attrs={'units': 'm/s'}),
        }

    def get_component_instance(self):
        return HeldSuarez()


if __name__ == '__main__':
    pytest.main([__file__])

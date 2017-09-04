import pytest
import abc
import os
from glob import glob
import xarray as xr
import numpy as np
from climt import (
    HeldSuarez, GrayLongwaveRadiation,
    Frierson06LongwaveOpticalDepth, GridScaleCondensation,
    BergerSolarInsolation, SimplePhysics, RRTMGLongwave,
    RRTMGShortwave, SlabSurface, EmanuelConvection,
    DcmipInitialConditions, GfsDynamicalCore, ClimtSpectralDynamicalCore,
    IceSheet)
import climt
from sympl import (
    DataArray, Implicit, TimeStepper, set_dimension_names
)
from datetime import datetime, timedelta
os.environ['NUMBA_DISABLE_JIT'] = '1'

vertical_dimension_names = [
    'interface_levels', 'mid_levels', 'full_levels']
set_dimension_names(x='longitude', y='latitude', z=vertical_dimension_names)

cache_folder = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'cached_component_output')


def cache_dictionary(dictionary, filename):
    dataset = xr.Dataset(dictionary)
    dataset.to_netcdf(filename, engine='scipy')


def load_dictionary(filename):
    dataset = xr.open_dataset(filename, engine='scipy')
    return dict(dataset.data_vars)


def state_3d_to_1d(state):
    return_state = {}
    for name, value in state.items():
        if name is 'time':
            return_state[name] = value
        else:
            dim_list = []
            for i, dim in enumerate(value.dims):
                if dim in vertical_dimension_names:
                    dim_list.append(slice(0, value.shape[i]))
                else:
                    dim_list.append(0)
            return_state[name] = value[tuple(dim_list)]
    return return_state


def transpose_state(state, dims=None):
    return_state = {}
    for name, value in state.items():
        if name is 'time':
            return_state[name] = state[name]
        else:
            if dims is None:
                return_state[name] = state[name].transpose()
            else:
                return_state[name] = state[name].transpose(*dims)
    return return_state


def call_with_timestep_if_needed(
        component, state, timestep=timedelta(seconds=10.)):

    if isinstance(component, IceSheet):
        output, diag = component(state, timestep=timestep)
        diag.pop('snow_ice_temperature_poly')
        return output, diag
    if isinstance(component, ClimtSpectralDynamicalCore):
        return component(state)
    elif isinstance(component, (Implicit, TimeStepper)):
        return component(state, timestep=timestep)
    else:
        return component(state)


class ComponentBase(object):

    @abc.abstractmethod
    def get_3d_input_state(self):
        pass

    @abc.abstractmethod
    def get_component_instance(self, state_modification_func=lambda x: x):
        pass

    def get_cached_output(self,):
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
        state = self.get_3d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output()
        if cached_output is None:
            self.cache_output(output)
            raise AssertionError(
                'Failed due to no cached output, cached current output')
        else:
            compare_outputs(output, cached_output)

    def test_1d_output_matches_cached_output(self):
        state = state_3d_to_1d(self.get_3d_input_state())
        component = self.get_component_instance(
            state_modification_func=state_3d_to_1d)
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output()
        if cached_output is None:
            raise AssertionError(
                'Failed due to no cached output.')
        else:
            if isinstance(cached_output, dict):
                compare_outputs(output, state_3d_to_1d(cached_output))
            else:
                cached_output_1d = []
                for cached_state in cached_output:
                    cached_output_1d.append(state_3d_to_1d(cached_state))
                compare_outputs(output, tuple(cached_output_1d))

    def test_reversed_state_gives_same_output(self):
        state = transpose_state(self.get_3d_input_state())
        component = self.get_component_instance(
            state_modification_func=transpose_state)
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output()
        if cached_output is None:
            raise AssertionError(
                'Failed due to no cached output.')
        else:
            compare_outputs(output, cached_output)

    def test_component_listed_inputs_are_accurate(self):
        state = self.get_3d_input_state()
        component = self.get_component_instance()
        input_state = {}
        for key in component.inputs:
            input_state[key] = state[key]
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output()
        if cached_output is not None:
            compare_outputs(output, cached_output)

    def test_consistent_dim_length(self):
        input_state = self.get_3d_input_state()
        assert_dimension_lengths_are_consistent(input_state)
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, input_state)
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
        if name != 'time':
            for i, dim_name in enumerate(value.dims):
                try:
                    if dim_name in dimension_lengths:
                        assert dimension_lengths[dim_name] == value.shape[i]
                    else:
                        dimension_lengths[dim_name] = value.shape[i]
                except AssertionError as err:
                    raise AssertionError(
                        'Inconsistent length on dimension {} for value {}:'
                        '{}'.format(dim_name, name, err))


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
        try:
            assert np.all(np.isclose(current[key].values, cached[key].values))
            for attr in current[key].attrs:
                assert current[key].attrs[attr] == cached[key].attrs[attr]
            for attr in cached[key].attrs:
                assert attr in current[key].attrs
            assert current[key].dims == cached[key].dims
        except AssertionError as err:
            raise AssertionError('Error for {}: {}'.format(key, err))
    for key in cached.keys():
        assert key in current.keys()


class TestHeldSuarez(ComponentBase):

    def get_3d_input_state(self):
        random = np.random.RandomState(0)
        return {
            'latitude': DataArray(
                random.randn(3), dims=['latitude'], attrs={'units': 'degrees_N'}),
            'air_pressure': DataArray(
                random.rand(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'hPa'},),
            'surface_air_pressure': DataArray(
                random.rand(2, 3), dims=['longitude', 'latitude'],
                attrs={'units': 'hPa'},),
            'air_temperature': DataArray(
                270. + random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'degK'}),
            'eastward_wind': DataArray(
                random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'm/s'}),
            'northward_wind': DataArray(
                random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'm/s'}),
        }

    def get_component_instance(self, state_modification_func=lambda x: x):
        return HeldSuarez()


class TestHeldSuarezCachedCoordinates(ComponentBase):

    def get_3d_input_state(self):
        random = np.random.RandomState(0)
        return {
            'latitude': DataArray(
                np.linspace(-90, 90, num=3),
                dims=['latitude'], attrs={'units': 'degrees_N'}),
            'air_pressure': DataArray(
                random.rand(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'hPa'},),
            'surface_air_pressure': DataArray(
                random.rand(2, 3), dims=['longitude', 'latitude'],
                attrs={'units': 'hPa'},),
            'air_temperature': DataArray(
                270. + random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'degK'}),
            'eastward_wind': DataArray(
                random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'm/s'}),
            'northward_wind': DataArray(
                random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
                attrs={'units': 'm/s'}),
        }

    def get_component_instance(self, state_modification_func=lambda x: x):
        return HeldSuarez()

    def test_1d_output_matches_cached_output(self):
        state = state_3d_to_1d(self.get_3d_input_state())
        component = self.get_component_instance()
        output = component(state)
        cached_output = self.get_cached_output()
        if cached_output is None:
            raise AssertionError(
                'Failed due to no cached output.')
        else:
            if isinstance(cached_output, dict):
                compare_outputs(output, state_3d_to_1d(cached_output))
            else:
                cached_output_1d = []
                for cached_state in cached_output:
                    cached_output_1d.append(state_3d_to_1d(cached_state))
                compare_outputs(output, tuple(cached_output_1d))


def test_hs_without_latitude():

    hs = HeldSuarez()

    random = np.random.RandomState(0)
    input_state = {
        'air_pressure': DataArray(
            random.rand(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
            attrs={'units': 'hPa'},),
        'surface_air_pressure': DataArray(
            random.rand(2, 3), dims=['longitude', 'latitude'],
            attrs={'units': 'hPa'},),
        'air_temperature': DataArray(
            270. + random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
            attrs={'units': 'degK'}),
        'eastward_wind': DataArray(
            random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
            attrs={'units': 'm/s'}),
        'northward_wind': DataArray(
            random.randn(2, 3, 6), dims=['longitude', 'latitude', 'mid_levels'],
            attrs={'units': 'm/s'})}

    with pytest.raises(IndexError) as excinfo:
        hs(input_state)
    assert 'quantity labeled' in str(excinfo.value)


class TestFrierson06LongwaveOpticalDepth(ComponentBase):

    def get_3d_input_state(self):
        state = {
            'latitude': DataArray(
                np.linspace(-90, 90, num=10),
                dims=['latitude'], attrs={'units': 'degrees_N'}),
            'sigma_on_interface_levels': DataArray(
                np.linspace(0., 1., num=6),
                dims=['interface_levels'], attrs={'units': ''}),
        }
        return state

    def get_1d_input_state(self):
        state_3d = self.get_3d_input_state()
        return {
            'latitude': state_3d['latitude'][0],
            'sigma_on_interface_levels': state_3d['sigma_on_interface_levels']
        }

    def get_component_instance(self, state_modification_func=lambda x: x):
        return Frierson06LongwaveOpticalDepth()


class TestGrayLongwaveRadiation(ComponentBase):

    def get_component_instance(self, state_modification_func=lambda x: x):
        return GrayLongwaveRadiation()

    def get_3d_input_state(self):
        component = self.get_component_instance()
        state = climt.get_default_state([component])
        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestGridScaleCondensation(ComponentBase):

    def get_component_instance(self, state_modification_func=lambda x: x):
        return GridScaleCondensation()

    def get_3d_input_state(self):
        nx, ny, nz = 2, 3, 10
        p_interface = DataArray(
            np.linspace(1e5, 0, nz+1),
            dims=['interface_levels'], attrs={'units': 'Pa'},
        )
        p = DataArray(
            0.5*(p_interface.values[1:] + p_interface.values[:-1]),
            dims=['mid_levels'], attrs={'units': 'Pa'})
        random = np.random.RandomState(0)
        return {
            'air_pressure': p,
            'air_temperature': DataArray(
                270. + random.randn(nx, ny, nz),
                dims=['longitude', 'latitude', 'mid_levels'], attrs={'units': 'degK'}),
            'specific_humidity': DataArray(
                random.rand(nx, ny, nz)*15.,
                dims=['longitude', 'latitude', 'mid_levels'], attrs={'units': 'kg/kg'}),
            'air_pressure_on_interface_levels': p_interface,
        }


class TestBergerSolarInsolation(ComponentBase):

    def get_component_instance(self, state_modification_func=lambda x: x):
        return BergerSolarInsolation()

    def get_3d_input_state(self):
        nx = 5
        ny = 10
        return {
            'time': datetime(2016, 12, 20, 6),
            'longitude': DataArray(
                np.linspace(-90, 90, nx, endpoint=False),
                dims=['longitude'], attrs={'units': 'degree_E'}),
            'latitude': DataArray(
                np.linspace(-180., 180., num=ny),
                dims=['latitude'], attrs={'units': 'degrees_north'}),
        }


class TestBergerSolarInsolationDifferentTime(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return BergerSolarInsolation()

    def get_3d_input_state(self):
        nx = 5
        ny = 10
        return {
            'time': datetime(1916, 12, 20, 6),
            'longitude': DataArray(
                np.linspace(-90, 90, nx, endpoint=False),
                dims=['longitude'], attrs={'units': 'degree_E'}),
            'latitude': DataArray(
                np.linspace(-180., 180., num=ny),
                dims=['latitude'], attrs={'units': 'degrees_north'}),
        }


class TestBergerSolarInsolationWithSolarConstant(ComponentBase):

    def get_component_instance(self, state_modification_func=lambda x: x):
        return BergerSolarInsolation()

    def get_3d_input_state(self):
        nx = 5
        ny = 10
        return {
            'time': datetime(2016, 12, 20, 6),
            'longitude': DataArray(
                np.linspace(-90, 90, nx, endpoint=False),
                dims=['longitude'], attrs={'units': 'degree_E'}),
            'latitude': DataArray(
                np.linspace(-180., 180., num=ny),
                dims=['latitude'], attrs={'units': 'degrees_north'}),
            'solar_constant': DataArray(
                1364.*np.ones(1),
                dims=[''], attrs={'units': 'W m^-2'}),

        }


def test_berger_insolation_with_bad_solar_constant():

    berger = BergerSolarInsolation()
    nx = 5
    ny = 10
    input_state = {
        'time': datetime(2016, 12, 20, 6),
        'longitude': DataArray(
            np.linspace(-90, 90, nx, endpoint=False),
            dims=['longitude'], attrs={'units': 'degree_E'}),
        'latitude': DataArray(
            np.linspace(-180., 180., num=ny),
            dims=['latitude'], attrs={'units': 'degrees_north'}),
        'solar_constant': DataArray(
            1364.*np.ones((2)),
            dims=['latitude'], attrs={'units': 'W m^-2'}),

    }

    with pytest.raises(ValueError) as excinfo:
        berger(input_state)
    assert 'Solar constant should' in str(excinfo.value)


class TestSimplePhysics(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return SimplePhysics()

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestSimplePhysicsPrognostic(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        component = SimplePhysics()
        return component.get_prognostic_version(timedelta(minutes=10))

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestRRTMGLongwave(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return RRTMGLongwave()

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestRRTMGLongwaveWithClouds(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return RRTMGLongwave(cloud_optical_properties=1)

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestRRTMGShortwave(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return RRTMGShortwave()

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestSlabSurface(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return SlabSurface()

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestEmanuel(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        emanuel = EmanuelConvection()
        emanuel.current_time_step = timedelta(seconds=300)
        return emanuel

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestDcmip(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return DcmipInitialConditions()

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


def testDcmipOptions():

    dcmip = DcmipInitialConditions()

    state = climt.get_default_state([dcmip],
                                    y=dict(label='latitude',
                                           values=np.linspace(0, 60, 20),
                                           units='degrees_north'))

    dry_state = dcmip(state)
    moist_state = dcmip(state, moist_simulation=True)
    not_perturbed_state = dcmip(state, add_perturbation=False)

    tropical_cyclone_state = dcmip(state, type_of_output='tropical_cyclone',
                                   moist_simulation=True)

    assert not np.all(np.isclose(dry_state['specific_humidity'].values,
                                 moist_state['specific_humidity'].values))

    assert not np.all(np.isclose(dry_state['eastward_wind'].values,
                                 not_perturbed_state['eastward_wind'].values))

    assert np.all(np.isclose(tropical_cyclone_state['surface_air_pressure'].values - 1.015e5,
                             np.zeros(not_perturbed_state['surface_air_pressure'].values.shape)))


class TestGfsDycore(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        return GfsDynamicalCore(number_of_longitudes=68,
                                number_of_latitudes=32)

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component], x=component.grid_definition['x'],
            y=component.grid_definition['y'],
            mid_levels=component.grid_definition['mid_levels'],
            interface_levels=component.grid_definition['interface_levels'])

        dcmip = climt.DcmipInitialConditions()
        out = dcmip(state, add_perturbation=True)
        state.update(out)

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestIceSheet(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        ice = IceSheet()
        return ice

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            x=dict(label='longtiude', values=np.linspace(0, 2, 4), units='degrees_east'),
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestIceSheetSeaIce(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        ice = IceSheet()
        return ice

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            x=dict(label='longtiude', values=np.linspace(0, 2, 4), units='degrees_east'),
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        state['area_type'][:] = 'sea_ice'
        state['sea_ice_thickness'][:] = 4
        state['surface_snow_thickness'][:] = 3

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestIceSheetLandIce(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        ice = IceSheet()
        return ice

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            x=dict(label='longtiude', values=np.linspace(0, 2, 4), units='degrees_east'),
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        state['area_type'][:] = 'land_ice'
        state['land_ice_thickness'][:] = 4
        state['surface_snow_thickness'][:] = 3

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


class TestIceSheetLand(ComponentBase):
    def get_component_instance(self, state_modification_func=lambda x: x):
        ice = IceSheet()
        return ice

    def get_3d_input_state(self):

        component = self.get_component_instance()
        state = climt.get_default_state(
            [component],
            x=dict(label='longtiude', values=np.linspace(0, 2, 4), units='degrees_east'),
            y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

        state['area_type'][:] = 'land'
        state['surface_snow_thickness'][:] = 3

        return state

    def test_1d_output_matches_cached_output(self):
        assert True


def test_ice_sheet_too_high():

    ice = IceSheet()

    state = climt.get_default_state(
        [ice],
        x=dict(label='longtiude', values=np.linspace(0, 2, 4), units='degrees_east'),
        y=dict(label='latitude', values=np.linspace(0, 2, 4), units='degrees_north'))

    state['area_type'][:] = 'land_ice'
    state['land_ice_thickness'][:] = 8
    state['surface_snow_thickness'][:] = 3

    with pytest.raises(ValueError) as excinfo:
        ice(state, timedelta(seconds=100))

    assert 'exceeds maximum value' in str(excinfo.value)


if __name__ == '__main__':
    pytest.main([__file__])

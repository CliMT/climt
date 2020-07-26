import pytest
import abc
import os
import sys
from glob import glob
import xarray as xr
import numpy as np
import logging
from climt import (
    HeldSuarez, GrayLongwaveRadiation,
    Frierson06LongwaveOpticalDepth, GridScaleCondensation,
    BergerSolarInsolation, SimplePhysics, RRTMGLongwave,
    RRTMGShortwave, SlabSurface, EmanuelConvection,
    DcmipInitialConditions, GFSDynamicalCore, BucketHydrology,
    IceSheet, Instellation, DryConvectiveAdjustment,
    get_grid)
import climt
from sympl import (
    Stepper, TendencyStepper, TimeDifferencingWrapper,
    ImplicitTendencyComponent, UpdateFrequencyWrapper, DataArray,
    TendencyComponent, AdamsBashforth
)
from sympl._core.tracers import reset_tracers, reset_packers
from datetime import datetime, timedelta
os.environ['NUMBA_DISABLE_JIT'] = '1'

vertical_dimension_names = [
    'interface_levels', 'mid_levels', 'full_levels']

cache_folder = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'cached_component_output')


def cache_dictionary(dictionary, filename):
    dataset = xr.Dataset(dictionary)
    dataset.to_netcdf(filename, engine='scipy')


def load_dictionary(filename):
    dataset = xr.open_dataset(filename, engine='scipy')
    return_dict = dict(dataset.data_vars)
    return_dict.update(dataset.coords)
    return return_dict


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
    np.random.seed(0)
    if isinstance(component, (Stepper, TendencyStepper, ImplicitTendencyComponent)):
        return component(state, timestep=timestep)
    else:
        return component(state)


class ComponentBase(object):

    def setUp(self):
        reset_tracers()
        reset_packers()
        super(ComponentBase, self).setUp()

    @abc.abstractmethod
    def get_component_instance(self):
        pass

    def get_cache_filename(self, descriptor, i):
        return '{}-{}-{}.cache'.format(self.__class__.__name__, descriptor, i)

    def get_cached_output(self, descriptor):
        cache_filename_list = sorted(glob(
            os.path.join(
                cache_folder,
                self.get_cache_filename(descriptor, '*'))))
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

    def cache_output(self, output, descriptor):
        if not isinstance(output, tuple):
            output = (output,)
        for i in range(len(output)):
            cache_filename = os.path.join(
                cache_folder, self.get_cache_filename(descriptor, i))
            cache_dictionary(output[i], cache_filename)

    def assert_valid_output(self, output):
        if isinstance(output, dict):
            output = [output]
        for i, out_dict in enumerate(output):
            for name, value in out_dict.items():
                try:
                    if name != 'time' and np.any(np.isnan(value)):
                        raise AssertionError(
                            'NaN produced in output {} from dict {}'.format(name, i))
                except TypeError:  # raised if cannot run isnan on dtype of value
                    pass


class ComponentBaseColumn(ComponentBase):

    def get_1d_input_state(self, component=None):
        if component is None:
            component = self.get_component_instance()
        return climt.get_default_state(
            [component], grid_state=get_grid(nx=None, ny=None, nz=30))

    def test_column_output_matches_cached_output(self):
        state = self.get_1d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output('column')
        if cached_output is None:
            self.cache_output(output, 'column')
            raise AssertionError(
                'Failed due to no cached output, cached current output.')
        else:
            compare_outputs(output, cached_output)

    def test_no_nans_in_column_output(self):
        state = self.get_1d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        self.assert_valid_output(output)

    def test_column_stepping_output_matches_cached_output(self):
        component = self.get_component_instance()
        if isinstance(component, (TendencyComponent, ImplicitTendencyComponent)):
            component = AdamsBashforth(self.get_component_instance())
            state = self.get_1d_input_state(component)
            output = call_with_timestep_if_needed(component, state)
            cached_output = self.get_cached_output('column_stepping')
            if cached_output is None:
                self.cache_output(output, 'column_stepping')
                raise AssertionError(
                    'Failed due to no cached output, cached current output.')
            else:
                compare_outputs(output, cached_output)


class ComponentBase3D(ComponentBase):

    def get_3d_input_state(self, component=None):
        if component is None:
            component = self.get_component_instance()
        return climt.get_default_state(
            [component], grid_state=get_grid(nx=32, ny=16, nz=28))

    def test_3d_output_matches_cached_output(self):
        state = self.get_3d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output('3d')
        if cached_output is None:
            self.cache_output(output, '3d')
            raise AssertionError(
                'Failed due to no cached output, cached current output.')
        else:
            compare_outputs(output, cached_output)

    def test_3d_stepping_output_matches_cached_output(self):
        component = self.get_component_instance()
        if isinstance(component, (TendencyComponent, ImplicitTendencyComponent)):
            component = AdamsBashforth(component)
            state = self.get_3d_input_state(component)
            output = call_with_timestep_if_needed(component, state)
            cached_output = self.get_cached_output('3d_stepping')
            if cached_output is None:
                self.cache_output(output, '3d_stepping')
                raise AssertionError(
                    'Failed due to no cached output, cached current output.')
            else:
                compare_outputs(output, cached_output)

    def test_no_nans_in_3D_output(self):
        state = self.get_3d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        self.assert_valid_output(output)

    def test_reversed_state_gives_same_output(self):
        state = self.get_3d_input_state()
        for name, value in state.items():
            if isinstance(value, (timedelta, datetime)):
                pass
            elif len(value.dims) == 3:
                state[name] = value.transpose(value.dims[2], value.dims[1], value.dims[0])
            elif len(value.dims) == 2:
                state[name] = value.transpose(value.dims[1], value.dims[0])
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output('3d')
        if cached_output is None:
            raise AssertionError(
                'Failed due to no cached output.')
        else:
            compare_outputs(output, cached_output)

    def test_transposed_state_gives_same_output(self):
        state = self.get_3d_input_state()
        for name, value in state.items():
            if isinstance(value, (timedelta, datetime)):
                pass
            elif len(value.dims) == 3:
                state[name] = value.transpose(value.dims[2], value.dims[0], value.dims[1])
            elif len(value.dims) == 2:
                state[name] = value.transpose(value.dims[1], value.dims[0])
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output('3d')
        if cached_output is None:
            raise AssertionError(
                'Failed due to no cached output.')
        else:
            compare_outputs(output, cached_output)


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
        if key == 'time':
            assert key in cached.keys()
        else:
            try:
                if not np.all(current[key] == cached[key]):
                    assert np.all(np.isclose(current[key] - cached[key], 0.))
                for attr in current[key].attrs:
                    assert current[key].attrs[attr] == cached[key].attrs[attr]
                for attr in cached[key].attrs:
                    assert attr in current[key].attrs
                assert set(current[key].dims) == set(cached[key].dims)
            except AssertionError as err:
                raise AssertionError('Error for {}: {}'.format(key, err))
    for key in cached.keys():
        assert key in current.keys()


class TestHeldSuarez(ComponentBase3D, ComponentBaseColumn):

    def get_component_instance(self):
        return HeldSuarez()


class TestFrierson06LongwaveOpticalDepth(ComponentBaseColumn, ComponentBase3D):

    def get_component_instance(self):
        return Frierson06LongwaveOpticalDepth()


class TestGrayLongwaveRadiation(ComponentBaseColumn, ComponentBase3D):

    def get_component_instance(self):
        return GrayLongwaveRadiation()


class TestGridScaleCondensation(ComponentBaseColumn, ComponentBase3D):

    def get_component_instance(self):
        return GridScaleCondensation()


class TestBergerSolarInsolation(ComponentBaseColumn, ComponentBase3D):

    def get_component_instance(self):
        return BergerSolarInsolation()

    def test_no_nans_in_2d_output(self):
        state = {
            'time': datetime(1998, 7, 13),
            'latitude': DataArray(
                np.linspace(-90, 90, 30),
                dims=['latitude'],
                attrs={'units': 'degrees_N'}
            ),
            'longitude': DataArray(
                np.linspace(0, 360, 60),
                dims=['longitude'],
                attrs={'units': 'degrees_E'}
            ),
        }
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        self.assert_valid_output(output)


class TestSimplePhysics(ComponentBaseColumn, ComponentBase3D):

    def get_component_instance(self):
        return SimplePhysics()


class TestSimplePhysicsImplicitPrognostic(ComponentBaseColumn, ComponentBase3D):

    def get_component_instance(self):
        component = TimeDifferencingWrapper(SimplePhysics())
        return component


class TestRRTMGLongwave(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGLongwave()


class TestRRTMGLongwaveMCICA(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGLongwave(mcica=True)

    def get_3d_input_state(self, component=None):
        if component is None:
            component = self.get_component_instance()
        state = climt.get_default_state([component],
                                        grid_state=climt.get_grid(nx=10, ny=5))
        state['cloud_area_fraction_in_atmosphere_layer'][16:19] = 0.5
        state['mass_content_of_cloud_ice_in_atmosphere_layer'][16:19] = 0.3
        return state

    def test_rrtmg_logging(self, caplog):
        caplog.set_level(logging.INFO)
        RRTMGLongwave(mcica=True, cloud_overlap_method='clear_only')
        assert 'no clouds' in caplog.text
        caplog.clear()

        RRTMGLongwave(mcica=True, cloud_optical_properties='single_cloud_type')
        assert "must be 'direct_input' or " \
               "'liquid_and_ice_clouds'" in caplog.text

    def test_transposed_state_gives_same_output(self):
        return

    def test_reversed_state_gives_same_output(self):
        return


class TestRRTMGLongwaveWithClouds(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGLongwave(cloud_optical_properties='single_cloud_type')


class TestRRTMGLongwaveWithExternalInterfaceTemperature(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGLongwave(calculate_interface_temperature=False)


class TestRRTMGShortwave(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGShortwave()


class TestRRTMGShortwaveMCICA(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGShortwave(mcica=True)

    def get_3d_input_state(self, component=None):
        if component is None:
            component = self.get_component_instance()
        state = climt.get_default_state([component],
                                        grid_state=climt.get_grid(nx=3, ny=2, nz=15))
        state['cloud_area_fraction_in_atmosphere_layer'][10:12] = 0.5
        state['mass_content_of_cloud_ice_in_atmosphere_layer'][10:12] = 0.3
        return state

    def test_transposed_state_gives_same_output(self):
        return

    def test_reversed_state_gives_same_output(self):
        return

    def test_rrtmg_logging(self, caplog):
        caplog.set_level(logging.INFO)
        RRTMGShortwave(mcica=True, cloud_overlap_method='clear_only')
        assert 'no clouds' in caplog.text
        caplog.clear()

        RRTMGShortwave(mcica=True, cloud_optical_properties='single_cloud_type')
        assert "must be 'direct_input' or " \
               "'liquid_and_ice_clouds'" in caplog.text
        caplog.clear()

        RRTMGShortwave(mcica=True,
                       cloud_optical_properties='liquid_and_ice_clouds',
                       cloud_ice_properties='ebert_curry_one')
        assert "not be set to 'ebert_curry_one'" in caplog.text
        caplog.clear()

        RRTMGShortwave(mcica=True,
                       cloud_optical_properties='liquid_and_ice_clouds',
                       cloud_liquid_water_properties='radius_independent_absorption')
        assert "must be set to 'radius_dependent_absorption'" in caplog.text


class TestSlabSurface(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return SlabSurface()


class TestBucketHydrology(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return BucketHydrology()


class TestEmanuel(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        emanuel = EmanuelConvection()
        return emanuel


class TestDcmip(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return DcmipInitialConditions()


def test_dcmip_options():

    state = climt.get_default_state([DcmipInitialConditions()],
                                    grid_state=get_grid(nx=64, ny=64, nz=10))

    dry_state = DcmipInitialConditions(moist=False)(state)
    moist_state = DcmipInitialConditions(moist=True)(state)
    not_perturbed_state = DcmipInitialConditions(moist=False, add_perturbation=False)(state)
    tropical_cyclone_state = DcmipInitialConditions(moist=True, condition_type='tropical_cyclone')(state)

    assert not np.all(np.isclose(dry_state['specific_humidity'].values,
                                 moist_state['specific_humidity'].values))

    assert not np.all(np.isclose(dry_state['eastward_wind'].values,
                                 not_perturbed_state['eastward_wind'].values))

    assert not np.all(np.isclose(tropical_cyclone_state['surface_air_pressure'].values - 1.015e5,
                                 np.zeros(not_perturbed_state['surface_air_pressure'].values.shape)))


class TestIceSheet(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return IceSheet()


class TestIceSheetLand(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        ice = IceSheet()
        return ice

    def get_3d_input_state(self):
        state = super(TestIceSheetLand, self).get_3d_input_state()

        state['area_type'].values[:] = 'land'
        state['surface_snow_thickness'].values[:] = 3

        return state
#
#
# def test_ice_sheet_too_high():
#
#     ice = IceSheet()
#
#     state_array = climt.get_default_state([ice])
#
#     state_array['area_type'].values = 'land_ice'
#     state_array['land_ice_thickness'].values = 8
#     state_array['surface_snow_thickness'].values = 3
#
#     with pytest.raises(ValueError) as excinfo:
#         ice(state_array, timedelta(seconds=100))
#
#     assert 'exceeds maximum value' in str(excinfo.value)


class TestInstellation(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return Instellation()


class TestDryConvection(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return DryConvectiveAdjustment()


@pytest.mark.skip("fails on CI, no idea why")
class TestFullMoistGFSDycoreWithPhysics(ComponentBase3D):

    def get_component_instance(self):
        # Create Radiation Prognostic
        radiation = climt.RRTMGLongwave()
        # Create Convection Prognostic
        convection = climt.EmanuelConvection()
        # Create a SimplePhysics Prognostic
        boundary_layer = TimeDifferencingWrapper(
            climt.SimplePhysics()
        )
        return GFSDynamicalCore(
            [radiation, convection, boundary_layer]
        )


@pytest.mark.skip("fails on CI, no idea why")
class TestGFSDycore(ComponentBase3D):

    def get_component_instance(self):
        return GFSDynamicalCore()


@pytest.mark.skip("fails on CI, no idea why")
class TestGFSDycoreWithDcmipInitialConditions(ComponentBase3D):

    def get_component_instance(self):
        return GFSDynamicalCore()

    def get_3d_input_state(self):
        state = climt.get_default_state(
            [self.get_component_instance()], grid_state=get_grid(nx=32, ny=32, nz=28))
        # state = super(TestGFSDycoreWithDcmipInitialConditions, self).get_3d_input_state()
        state.update(climt.DcmipInitialConditions(add_perturbation=True)(state))
        return state


@pytest.mark.skipif(sys.platform == 'win32',
                    reason="fails on appveyor, no idea why")
class TestGFSDycoreWithImplicitTendency(ComponentBase3D):

    def get_component_instance(self):
        return GFSDynamicalCore([EmanuelConvection()])

    def get_3d_input_state(self):
        state = climt.get_default_state(
            [self.get_component_instance()], grid_state=get_grid(nx=16, ny=16, nz=28))
        return state


@pytest.mark.skipif(sys.platform == 'win32',
                    reason="fails on appveyor, no idea why")
class TestGFSDycoreWithHeldSuarez(ComponentBase3D):
    def test_inputs_are_dry(self):
        component = self.get_component_instance()
        assert 'specific_humidity' not in component.input_properties.keys()
        assert 'specific_humidity_on_interface_levels' not in component.input_properties.keys()

    def get_component_instance(self):
        return GFSDynamicalCore([HeldSuarez()])

    def get_3d_input_state(self):
        state = climt.get_default_state(
            [self.get_component_instance()],
            grid_state=get_grid(nx=16, ny=16, nz=28))
        return state


@pytest.mark.skipif(sys.platform == 'win32',
                    reason="fails on appveyor, no idea why")
class TestGFSDycoreWithGrayLongwaveRadiation(ComponentBase3D):

    def get_component_instance(self):
        return GFSDynamicalCore([GrayLongwaveRadiation()])

    def get_3d_input_state(self):
        state = climt.get_default_state(
            [self.get_component_instance()],
            grid_state=get_grid(nx=16, ny=16, nz=28))
        return state


@pytest.mark.skipif(sys.platform == 'win32',
                    reason="fails on appveyor, no idea why")
class TestGFSDycoreWithRRTMGLongwave(ComponentBase3D):

    def get_component_instance(self):
        radiation = RRTMGLongwave()
        return GFSDynamicalCore([radiation], moist=True)

    def get_3d_input_state(self):
        state = climt.get_default_state(
            [self.get_component_instance()], grid_state=get_grid(nx=16, ny=16, nz=28))
        return state


def test_piecewise_constant_component():

    radiation = UpdateFrequencyWrapper(RRTMGLongwave(), timedelta(seconds=1000))

    state = climt.get_default_state([radiation])

    current_tendency, current_diagnostic = radiation(state)

    # Perturb state_array
    state['air_temperature'] += 3

    new_tendency, new_diagnostic = radiation(state)

    assert np.all(current_tendency['air_temperature'].values ==
                  new_tendency['air_temperature'].values)

    state['time'] += timedelta(seconds=1500)

    new_tendency, new_diagnostic = radiation(state)

    assert np.any(current_tendency['air_temperature'].values !=
                  new_tendency['air_temperature'].values)


if __name__ == '__main__':
    pytest.main([__file__])

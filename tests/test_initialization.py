from climt import (
    get_default_state, Frierson06LongwaveOpticalDepth, GrayLongwaveRadiation, HeldSuarez,
    GridScaleCondensation, BergerSolarInsolation, SimplePhysics, RRTMGLongwave, RRTMGShortwave,
    EmanuelConvection, SlabSurface, GFSDynamicalCore, DcmipInitialConditions, IceSheet,
    Instellation, get_grid
)
import random
from sympl import (
    TendencyComponent, DiagnosticComponent, Stepper,
    ImplicitTendencyComponent, TendencyStepper, TimeDifferencingWrapper
)
import numpy as np
import pytest
import unittest
from datetime import timedelta


def call_component(component, state):
    if isinstance(component, (DiagnosticComponent, TendencyComponent)):
        return component(state)
    elif isinstance(component, (Stepper, ImplicitTendencyComponent, TendencyStepper)):
        return component(state, timedelta(hours=1))
    else:
        raise AssertionError('Component is of unknown type')


class GetGridTests(unittest.TestCase):

    def assert_grid_quantities_present(self, state, latitude=False, longitude=False):
        grid_names = ['time', 'air_pressure', 'air_pressure_on_interface_levels',
                      'surface_air_pressure', 'height_on_ice_interface_levels']
        if latitude:
            grid_names.append('latitude')
        if longitude:
            grid_names.append('longitude')
        for name in grid_names:
            if name not in state:
                raise AssertionError(
                    'Grid quantity {} is not present in state_array.'.format(name))

    def assert_grid_quantities_have_dimensions(self, state, dim_names):
        state_dim_names = []
        for name, value in state.items():
            if name != 'time':
                state_dim_names.extend(value.dims)
        extra_dims = set(state_dim_names).difference(dim_names)
        missing_dims = set(dim_names).difference(state_dim_names)
        assert len(extra_dims) == 0
        assert len(missing_dims) == 0

    def assert_grid_dimensions_have_lengths(self, state, dim_lengths):
        state_dim_lengths = {}
        for name, value in state.items():
            if name != 'time':
                for dim_name, length in zip(value.dims, value.shape):
                    if dim_name not in state_dim_lengths:
                        state_dim_lengths[dim_name] = length
                    elif state_dim_lengths[dim_name] != length:
                        print(
                            {name: zip(value.dims, value.shape)
                             for name, value in state.items()
                             if name != 'time'})
                        raise AssertionError(
                            'Inconsistent lengths {} and {} for dimension {}'.format(
                                length, state_dim_lengths[dim_name], dim_name))
        for dim_name, length in dim_lengths.items():
            if state_dim_lengths[dim_name] != length:
                raise AssertionError(
                    'Want length {} for dimension {} but instead have '
                    'length {}'.format(
                        length, dim_name, state_dim_lengths[dim_name])
                )

    def test_get_default_grid(self):
        grid = get_grid()
        self.assert_grid_quantities_present(grid)
        self.assert_grid_quantities_have_dimensions(
            grid, ['lat', 'lon', 'mid_levels', 'interface_levels', 'ice_interface_levels'])

    def test_get_1d_vertical_grid(self):
        grid = get_grid(nz=20)
        self.assert_grid_quantities_present(grid)
        self.assert_grid_quantities_have_dimensions(
            grid, ['lat', 'lon', 'mid_levels', 'interface_levels', 'ice_interface_levels'])
        self.assert_grid_dimensions_have_lengths(
            grid, {'mid_levels': 20, 'interface_levels': 21}
        )

    def test_get_3d_grid(self):
        grid = get_grid(nx=4, ny=6, nz=20)
        self.assert_grid_quantities_present(grid, latitude=True, longitude=True)
        self.assert_grid_quantities_have_dimensions(
            grid, ['mid_levels', 'interface_levels', 'lat', 'lon', 'ice_interface_levels'])
        self.assert_grid_dimensions_have_lengths(
            grid, {'mid_levels': 20, 'interface_levels': 21, 'lat': 6, 'lon': 4}
        )

    def test_get_3d_grid_custom_dim_names(self):
        grid = get_grid(nx=3, ny=8, nz=20, x_name='name1', y_name='name2')
        self.assert_grid_quantities_present(grid, latitude=True, longitude=True)
        self.assert_grid_quantities_have_dimensions(
            grid, ['mid_levels', 'interface_levels', 'name1', 'name2', 'ice_interface_levels'])
        self.assert_grid_dimensions_have_lengths(
            grid, {'mid_levels': 20, 'interface_levels': 21, 'name1': 3, 'name2': 8}
        )

    def test_get_1d_grid_custom_surface_pressure(self):
        grid = get_grid(nz=20, p_surf_in_Pa=0.9e5)
        self.assert_grid_quantities_present(grid)
        self.assert_grid_quantities_have_dimensions(
            grid, ['lat', 'lon', 'mid_levels', 'interface_levels', 'ice_interface_levels'])
        self.assert_grid_dimensions_have_lengths(
            grid, {'mid_levels': 20, 'interface_levels': 21}
        )
        p = grid['air_pressure'].to_units('Pa')
        p_interface = grid['air_pressure_on_interface_levels'].to_units('Pa')
        assert grid['surface_air_pressure'].to_units('Pa') == 0.9e5
        assert np.isclose(p_interface[0], 0.9e5)
        assert np.all(p_interface[1:].values < p_interface[:-1].values)
        assert np.all(p[1:].values < p[:-1].values)
        assert np.all(p[:].values < p_interface[:-1].values)
        assert np.all(p[:].values > p_interface[1:].values)


def assert_state_is_full(state, component):
    for dict_name in ('output_properties', 'tendency_properties', 'diagnostic_proprerties'):
        if hasattr(component, dict_name):
            for quantity_name, properties in getattr(component, dict_name).items():
                if quantity_name not in component.input_properties.keys():
                    continue
                elif 'dims' in properties.keys():
                    dims = properties['dims']
                else:
                    dims = component.input_properties[quantity_name]['dims']
                missing_dims = set(dims).difference(
                    ['*'] + list(state[quantity_name].dims))
                assert len(missing_dims) == 0, '{} is missing {} dims {}'.format(
                    quantity_name, dict_name, missing_dims)


def create_default_test_for(cls):
    def test_component(self):
        component = cls()
        state = get_default_state([component])
        call_component(component, state)
    test_component.__name__ = 'test_{}'.format(cls.__name__)
    return test_component


def create_1d_grid_test_for(cls):
    def test_component_1d_grid(self):
        grid = get_grid(nz=10)
        component = cls()
        state = get_default_state([component], grid_state=grid)
        assert_state_is_full(state, component)
        call_component(component, state)
    test_component_1d_grid.__name__ = 'test_{}_1d_grid'.format(cls.__name__)
    return test_component_1d_grid


def create_2d_grid_test_for(cls):
    def test_component_2d_grid(self):
        grid = get_grid(nx=3, nz=10)
        component = cls()
        state = get_default_state([component], grid_state=grid)
        assert_state_is_full(state, component)
        call_component(component, state)
    test_component_2d_grid.__name__ = 'test_{}_2d_grid'.format(cls.__name__)
    return test_component_2d_grid


def create_3d_grid_test_for(cls):
    def test_component_3d_grid(self):
        grid = get_grid(nx=3, ny=4, nz=10)
        component = cls()
        state = get_default_state([component], grid_state=grid)
        assert_state_is_full(state, component)
        call_component(component, state)
    test_component_3d_grid.__name__ = 'test_{}_3d_grid'.format(cls.__name__)
    return test_component_3d_grid


class ComponentQuantityInitializationTests(unittest.TestCase):

    component_classes = (
        Frierson06LongwaveOpticalDepth, GrayLongwaveRadiation, HeldSuarez,
        GridScaleCondensation, BergerSolarInsolation, SimplePhysics,
        RRTMGLongwave,
        RRTMGShortwave,
        EmanuelConvection, SlabSurface,
        DcmipInitialConditions, IceSheet,
        Instellation
    )

    pair_tests = 20
    triplet_tests = 20

    for cls in component_classes:
        func = create_default_test_for(cls)
        locals()[func.__name__] = func
        func = create_1d_grid_test_for(cls)
        locals()[func.__name__] = func
        func = create_2d_grid_test_for(cls)
        locals()[func.__name__] = func
        func = create_3d_grid_test_for(cls)
        locals()[func.__name__] = func

    def test_GFSDynamicalCore(self):
        grid = get_grid(nx=12, ny=16, nz=28)
        component = GFSDynamicalCore()
        state = get_default_state([component], grid_state=grid)
        call_component(component, state)

    def test_component_pairs(self):
        random.seed(0)
        for _ in range(self.pair_tests):
            i, j = random.sample(range(len(self.component_classes)), 2)
            component1 = self.component_classes[i]()
            component2 = self.component_classes[j]()
            print(component1.__class__.__name__, component2.__class__.__name__)
            state = get_default_state([component1, component2])
            call_component(component1, state)
            call_component(component2, state)

    def test_component_triplets(self):
        random.seed(0)
        for _ in range(self.triplet_tests):
            i, j, k = random.sample(range(len(self.component_classes)), 3)
            component1 = self.component_classes[i]()
            component2 = self.component_classes[j]()
            component3 = self.component_classes[k]()
            print(component1.__class__.__name__, component2.__class__.__name__, component3.__class__.__name__)
            state = get_default_state([component1, component2, component3])
            call_component(component1, state)
            call_component(component2, state)
            call_component(component3, state)


class TestFullMoistGFSDycoreWithPhysics(unittest.TestCase):

    def get_component_instance(self):
        # Create Radiation Prognostic
        radiation = RRTMGLongwave()
        # Create Convection Prognostic
        convection = EmanuelConvection()
        # Create a SimplePhysics Prognostic
        boundary_layer = TimeDifferencingWrapper(
            SimplePhysics()
        )
        return GFSDynamicalCore(
            [radiation, convection, boundary_layer]
        )

    def test_component_3d_grid(self):
        grid = get_grid(nx=16, ny=16, nz=16)
        component = self.get_component_instance()
        state = get_default_state([component], grid_state=grid)
        call_component(component, state)


class TestGFSDycoreWith32VerticalLevels(unittest.TestCase):

    def get_component_instance(self):
        return GFSDynamicalCore()

    def test_component_3d_grid(self):
        grid = get_grid(nx=16, ny=16, nz=32)
        component = self.get_component_instance()
        state = get_default_state([component], grid_state=grid)
        call_component(component, state)


def test_3d_initialization_is_full_based_on_wildcard():
    grid_state = get_grid(nx=10, ny=10, nz=20)
    rrtmg_shortwave = RRTMGShortwave()
    instellation = Instellation()
    state = get_default_state([rrtmg_shortwave, instellation], grid_state=grid_state)
    for quantity_name, properties in rrtmg_shortwave.input_properties.items():
        if '*' in properties['dims']:
            assert len(state[quantity_name].dims) == len(properties['dims']) + 1
        if tuple(properties['dims']) == ('*',):
            assert set(state[quantity_name].dims) == {'lat', 'lon'}
        elif tuple(properties['dims']) == ('mid_levels', '*'):
            assert state[quantity_name].dims[0] == 'mid_levels'
            assert set(state[quantity_name].dims[1:]) == {'lat', 'lon'}


if __name__ == '__main__':
    pytest.main([__file__])

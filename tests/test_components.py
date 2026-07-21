import abc
import logging
import os
import sys
from datetime import datetime, timedelta
from glob import glob

import numpy as np
import pytest
import unyt
import xarray as xr
from sympl import (
    AdamsBashforth,
    DataArray,
    ImplicitTendencyComponent,
    Stepper,
    TendencyComponent,
    TendencyStepper,
    TimeDifferencingWrapper,
    UpdateFrequencyWrapper,
    set_backend,
)
from sympl._core.tracers import reset_packers, reset_tracers

import climt
from climt import (
    BergerSolarInsolation,
    BucketHydrology,
    DcmipInitialConditions,
    DryConvectiveAdjustment,
    EmanuelConvection,
    Frierson06LongwaveOpticalDepth,
    GrayLongwaveRadiation,
    GridScaleCondensation,
    HeldSuarez,
    IceSheet,
    Instellation,
    LandIce,
    LandMask,
    RRTMGLongwave,
    RRTMGShortwave,
    SeaIce,
    SimplePhysics,
    SlabSurface,
    UnytTimeDelta,
    get_grid,
)
from climt._core.unyt_backend import UnytBackend, UnytStateContainer

set_backend(UnytBackend())

# Note: NUMBA_DISABLE_JIT must NOT be set here — the components use @njit
# kernels that are decorated at import time but compiled lazily. Setting
# NUMBA_DISABLE_JIT after import corrupts Numba's internal state.

vertical_dimension_names = ["interface_levels", "mid_levels", "full_levels"]

cache_folder = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "cached_component_output"
)


def cache_dictionary(dictionary, filename):
    xr_dict = {}
    for key, value in dictionary.items():
        if key == "time":
            xr_dict[key] = value
        elif hasattr(value, "dims") and hasattr(value, "values"):
            xr_dict[key] = xr.DataArray(
                value.values, dims=value.dims, attrs=value.attrs
            )
        else:
            xr_dict[key] = value
    dataset = xr.Dataset(xr_dict)
    dataset.to_netcdf(filename, engine="scipy")


def load_dictionary(filename):
    import pandas as pd
    import sympl

    # Open with decode_times=False: on some numpy/xarray builds (observed on
    # the py3.12 CI runners) xarray's C-level CF-datetime decoder segfaults
    # while decoding even a trivial "days since 2000-01-01" scalar. We decode
    # the scalar time coordinate ourselves below, which avoids that native
    # crash and reproduces the exact same datetime the decoder would return.
    dataset = xr.open_dataset(filename, engine="scipy", decode_times=False)
    return_dict = {}
    backend = sympl.get_backend()

    for name, var in dataset.variables.items():
        if name == "time":
            val = np.asarray(var.values)
            units = var.attrs.get("units", "")
            if getattr(val, "ndim", 1) == 0 and "since" in units:
                # e.g. units="days since 2000-01-01 00:00:00"
                interval, ref = units.split("since")
                ref_date = pd.Timestamp(ref.strip())
                return_dict[name] = (
                    ref_date + pd.Timedelta(f"{float(val)} {interval.strip()}")
                ).to_pydatetime()
            else:
                return_dict[name] = var.values
        else:
            units = var.attrs.get("units", "")
            if hasattr(backend, "create_quantity"):
                return_dict[name] = backend.create_quantity(
                    var.values, name, units, var.dims
                )
            else:
                return_dict[name] = var

    return return_dict


def state_3d_to_1d(state):
    return_state = {}
    for name, value in state.items():
        if name == "time":
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
        if name == "time":
            return_state[name] = state[name]
        else:
            if dims is None:
                return_state[name] = state[name].transpose()
            else:
                return_state[name] = state[name].transpose(*dims)
    return return_state


def call_with_timestep_if_needed(
    component, state, timestep=UnytTimeDelta(seconds=10.0)
):
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
        return "{}-{}-{}.cache".format(self.__class__.__name__, descriptor, i)

    def get_cached_output(self, descriptor):
        cache_filename_list = sorted(
            glob(os.path.join(cache_folder, self.get_cache_filename(descriptor, "*")))
        )
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
                cache_folder, self.get_cache_filename(descriptor, i)
            )
            cache_dictionary(output[i], cache_filename)

    def assert_valid_output(self, output):
        if isinstance(output, dict):
            output = [output]
        for i, out_dict in enumerate(output):
            for name, value in out_dict.items():
                try:
                    if name != "time" and np.any(np.isnan(value)):
                        raise AssertionError(
                            "NaN produced in output {} from dict {}".format(name, i)
                        )
                except TypeError:  # raised if cannot run isnan on dtype of value
                    pass


class ComponentBaseColumn(ComponentBase):
    def get_1d_input_state(self, component=None):
        if component is None:
            component = self.get_component_instance()
        return climt.get_default_state(
            [component], grid_state=get_grid(nx=None, ny=None, nz=30)
        )

    def test_column_output_matches_cached_output(self):
        state = self.get_1d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output("column")
        if cached_output is None:
            self.cache_output(output, "column")
            raise AssertionError(
                "Failed due to no cached output, cached current output."
            )
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
            cached_output = self.get_cached_output("column_stepping")
            if cached_output is None:
                self.cache_output(output, "column_stepping")
                raise AssertionError(
                    "Failed due to no cached output, cached current output."
                )
            else:
                compare_outputs(output, cached_output)


class ComponentBase3D(ComponentBase):
    def get_3d_input_state(self, component=None):
        if component is None:
            component = self.get_component_instance()
        return climt.get_default_state(
            [component], grid_state=get_grid(nx=32, ny=16, nz=28)
        )

    def test_3d_output_matches_cached_output(self):
        state = self.get_3d_input_state()
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output("3d")
        if cached_output is None:
            self.cache_output(output, "3d")
            raise AssertionError(
                "Failed due to no cached output, cached current output."
            )
        else:
            compare_outputs(output, cached_output)

    def test_3d_stepping_output_matches_cached_output(self):
        component = self.get_component_instance()
        if isinstance(component, (TendencyComponent, ImplicitTendencyComponent)):
            component = AdamsBashforth(component)
            state = self.get_3d_input_state(component)
            output = call_with_timestep_if_needed(component, state)
            cached_output = self.get_cached_output("3d_stepping")
            if cached_output is None:
                self.cache_output(output, "3d_stepping")
                raise AssertionError(
                    "Failed due to no cached output, cached current output."
                )
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
                state[name] = value.transpose(
                    value.dims[2], value.dims[1], value.dims[0]
                )
            elif len(value.dims) == 2:
                state[name] = value.transpose(value.dims[1], value.dims[0])
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output("3d")
        if cached_output is None:
            raise AssertionError("Failed due to no cached output.")
        else:
            compare_outputs(output, cached_output)

    def test_transposed_state_gives_same_output(self):
        state = self.get_3d_input_state()
        for name, value in state.items():
            if isinstance(value, (timedelta, datetime)):
                pass
            elif len(value.dims) == 3:
                state[name] = value.transpose(
                    value.dims[2], value.dims[0], value.dims[1]
                )
            elif len(value.dims) == 2:
                state[name] = value.transpose(value.dims[1], value.dims[0])
        component = self.get_component_instance()
        output = call_with_timestep_if_needed(component, state)
        cached_output = self.get_cached_output("3d")
        if cached_output is None:
            raise AssertionError("Failed due to no cached output.")
        else:
            compare_outputs(output, cached_output)


def compare_outputs(current, cached):
    if isinstance(current, tuple) and isinstance(cached, tuple):
        for i in range(len(current)):
            compare_one_state_pair(current[i], cached[i])
    elif (not isinstance(current, tuple)) and (not isinstance(cached, tuple)):
        compare_one_state_pair(current, cached)
    else:
        raise AssertionError("Different number of dicts returned than cached.")


def compare_one_state_pair(current, cached):
    for key in current.keys():
        if key == "time":
            assert key in cached.keys()
        else:
            try:
                curr_arr = np.asarray(current[key])
                cach_arr = np.asarray(cached[key])
                if (
                    set(current[key].dims) == set(cached[key].dims)
                    and current[key].dims != cached[key].dims
                ):
                    perm = [current[key].dims.index(d) for d in cached[key].dims]
                    curr_arr = curr_arr.transpose(perm)

                if not np.all(curr_arr == cach_arr):
                    assert np.all(np.isclose(curr_arr - cach_arr, 0.0))
                for attr in current[key].attrs:
                    if attr == "units":
                        sanitized_curr = (
                            str(current[key].attrs[attr])
                            .replace("^", "**")
                            .replace(" ", "*")
                        )
                        sanitized_cach = (
                            str(cached[key].attrs[attr])
                            .replace("^", "**")
                            .replace(" ", "*")
                        )
                        assert unyt.Unit(sanitized_curr) == unyt.Unit(sanitized_cach)
                    else:
                        assert current[key].attrs[attr] == cached[key].attrs[attr]
                for attr in cached[key].attrs:
                    assert attr in current[key].attrs
                assert set(current[key].dims) == set(cached[key].dims)
            except AssertionError as err:
                raise AssertionError("Error for {}: {}".format(key, err))
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
            "time": datetime(1998, 7, 13),
            "latitude": UnytStateContainer(
                np.linspace(-90, 90, 30),
                dims=["latitude"],
                attrs={"units": "degrees_N"},
            ),
            "longitude": UnytStateContainer(
                np.linspace(0, 360, 60),
                dims=["longitude"],
                attrs={"units": "degrees_E"},
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
        state = climt.get_default_state(
            [component], grid_state=climt.get_grid(nx=10, ny=5)
        )
        state["cloud_area_fraction_in_atmosphere_layer"][16:19] = 0.5
        state["mass_content_of_cloud_ice_in_atmosphere_layer"][16:19] = 0.3
        return state

    def test_rrtmg_logging(self, caplog):
        caplog.set_level(logging.INFO)
        RRTMGLongwave(mcica=True, cloud_overlap_method="clear_only")
        assert "no clouds" in caplog.text
        caplog.clear()

        RRTMGLongwave(mcica=True, cloud_optical_properties="single_cloud_type")
        assert "must be 'direct_input' or 'liquid_and_ice_clouds'" in caplog.text

    def test_transposed_state_gives_same_output(self):
        return

    def test_reversed_state_gives_same_output(self):
        return


class TestRRTMGLongwaveWithClouds(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return RRTMGLongwave(cloud_optical_properties="single_cloud_type")


class TestRRTMGLongwaveWithExternalInterfaceTemperature(
    ComponentBaseColumn, ComponentBase3D
):
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
        state = climt.get_default_state(
            [component], grid_state=climt.get_grid(nx=3, ny=2, nz=15)
        )
        state["cloud_area_fraction_in_atmosphere_layer"][10:12] = 0.5
        state["mass_content_of_cloud_ice_in_atmosphere_layer"][10:12] = 0.3
        return state

    def test_transposed_state_gives_same_output(self):
        return

    def test_reversed_state_gives_same_output(self):
        return

    def test_rrtmg_logging(self, caplog):
        caplog.set_level(logging.INFO)
        RRTMGShortwave(mcica=True, cloud_overlap_method="clear_only")
        assert "no clouds" in caplog.text
        caplog.clear()

        RRTMGShortwave(mcica=True, cloud_optical_properties="single_cloud_type")
        assert "must be 'direct_input' or 'liquid_and_ice_clouds'" in caplog.text
        caplog.clear()

        RRTMGShortwave(
            mcica=True,
            cloud_optical_properties="liquid_and_ice_clouds",
            cloud_ice_properties="ebert_curry_one",
        )
        assert "not be set to 'ebert_curry_one'" in caplog.text
        caplog.clear()

        RRTMGShortwave(
            mcica=True,
            cloud_optical_properties="liquid_and_ice_clouds",
            cloud_liquid_water_properties="radius_independent_absorption",
        )
        assert "must be set to 'radius_dependent_absorption'" in caplog.text


class TestSlabSurface(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return SlabSurface()

    def test_column_stepping_output_matches_cached_output(self):
        component = self.get_component_instance()
        if isinstance(component, (TendencyComponent, ImplicitTendencyComponent)):
            component = AdamsBashforth(self.get_component_instance())
            state = self.get_1d_input_state(component)
            state["surface_material_density"].values[:] = state["sea_water_density"]
            output = call_with_timestep_if_needed(component, state)
            cached_output = self.get_cached_output("column_stepping")
            if cached_output is None:
                self.cache_output(output, "column_stepping")
                raise AssertionError(
                    "Failed due to no cached output, cached current output."
                )
            else:
                compare_outputs(output, cached_output)

    def test_3d_stepping_output_matches_cached_output(self):
        component = self.get_component_instance()
        if isinstance(component, (TendencyComponent, ImplicitTendencyComponent)):
            component = AdamsBashforth(component)
            state = self.get_3d_input_state(component)
            state["surface_material_density"].values[:] = state["sea_water_density"]
            output = call_with_timestep_if_needed(component, state)
            cached_output = self.get_cached_output("3d_stepping")
            if cached_output is None:
                self.cache_output(output, "3d_stepping")
                raise AssertionError(
                    "Failed due to no cached output, cached current output."
                )
            else:
                compare_outputs(output, cached_output)


class TestBucketHydrology(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return BucketHydrology()


class TestBucketHydrologyTwoLayer(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.BucketHydrology(num_layers=2,
                                     moisture_diffusion_timescale=86400.0)


class TestEmanuel(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        emanuel = EmanuelConvection()
        return emanuel


class TestDcmip(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return DcmipInitialConditions()


def test_dcmip_options():

    state = climt.get_default_state(
        [DcmipInitialConditions()], grid_state=get_grid(nx=64, ny=64, nz=10)
    )

    dry_state = DcmipInitialConditions(moist=False)(state)
    moist_state = DcmipInitialConditions(moist=True)(state)
    not_perturbed_state = DcmipInitialConditions(moist=False, add_perturbation=False)(
        state
    )
    tropical_cyclone_state = DcmipInitialConditions(
        moist=True, condition_type="tropical_cyclone"
    )(state)

    assert not np.all(
        np.isclose(
            dry_state["specific_humidity"].values,
            moist_state["specific_humidity"].values,
        )
    )

    assert not np.all(
        np.isclose(
            dry_state["eastward_wind"].values,
            not_perturbed_state["eastward_wind"].values,
        )
    )

    assert not np.all(
        np.isclose(
            tropical_cyclone_state["surface_air_pressure"].values - 1.015e5,
            np.zeros(not_perturbed_state["surface_air_pressure"].values.shape),
        )
    )


class TestIceSheet(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return IceSheet()


class TestIceSheetLand(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        ice = IceSheet()
        return ice

    def get_3d_input_state(self):
        state = super(TestIceSheetLand, self).get_3d_input_state()

        state["area_type"].values[:] = "land"
        state["surface_snow_thickness"].values[:] = 3

        return state


class TestSeaIce(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.SeaIce()

    def get_3d_input_state(self, component=None):
        state = super(TestSeaIce, self).get_3d_input_state(component)
        state["area_type"].values[:] = "sea_ice"
        state["sea_ice_thickness"].values[:] = 1.0
        return state


class TestLandIce(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.LandIce()

    def get_3d_input_state(self, component=None):
        state = super(TestLandIce, self).get_3d_input_state(component)
        state["area_type"].values[:] = "land_ice"
        state["land_ice_thickness"].values[:] = 3.0
        return state


class TestDataOcean(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        import os
        import tempfile

        # write a tiny fixed dataset once to a temp path stored on the class
        if not hasattr(self.__class__, "_sst_path"):
            lat = np.arange(-88.0, 90.0, 8.0)
            lon = np.arange(4.0, 360.0, 8.0)
            data = np.repeat(
                (290.0 + 0 * lat[:, None] + 0 * lon[None, :])[None], 12, 0
            )
            ds = xr.Dataset(
                {"tos": (("time", "lat", "lon"), data)},
                coords={"time": np.arange(12), "lat": lat, "lon": lon},
            )
            ds["tos"].attrs["units"] = "K"
            p = os.path.join(tempfile.gettempdir(), "climt_test_sst.nc")
            ds.to_netcdf(p, engine="scipy")
            self.__class__._sst_path = p
        return climt.DataOcean(self.__class__._sst_path, sst_variable="tos")

    def get_1d_input_state(self, component=None):
        state = super(TestDataOcean, self).get_1d_input_state(component)
        state["time"] = datetime(2000, 1, 15, 12)
        return state

    def get_3d_input_state(self, component=None):
        state = super(TestDataOcean, self).get_3d_input_state(component)
        state["time"] = datetime(2000, 1, 15, 12)
        return state


class TestSecondBEST(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.SecondBEST()

    def get_1d_input_state(self, component=None):
        state = super(TestSecondBEST, self).get_1d_input_state(component)
        state["area_type"].values[:] = "land"
        return state

    def get_3d_input_state(self, component=None):
        state = super(TestSecondBEST, self).get_3d_input_state(component)
        state["area_type"].values[:] = "land"
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

    def test_piecewise_constant_component(self):
        radiation = UpdateFrequencyWrapper(RRTMGLongwave(), UnytTimeDelta(seconds=1000))
        state = climt.get_default_state([radiation])
        current_tendency, current_diagnostic = radiation(state)
        # Perturb state_array
        state["air_temperature"] += unyt.unyt_array(
            3, state["air_temperature"].attrs["units"]
        )

        new_tendency, new_diagnostic = radiation(state)

        assert np.all(
            current_tendency["air_temperature"].values
            == new_tendency["air_temperature"].values
        )

        state["time"] += timedelta(seconds=1500)

        new_tendency, new_diagnostic = radiation(state)

        assert np.any(
            current_tendency["air_temperature"].values
            != new_tendency["air_temperature"].values
        )


class TestLandMask(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.LandMask()


if __name__ == "__main__":
    pytest.main([__file__])

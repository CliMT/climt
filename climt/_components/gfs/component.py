from __future__ import division
from ..._core import ensure_contiguous_state
from sympl import (
    get_constant, TendencyStepper, initialize_numpy_arrays_with_properties,
    get_numpy_arrays_with_properties, ImplicitTendencyComponentComposite,
    get_tracer_names, restore_data_arrays_with_properties,
)
import numpy as np
import sys
import logging
try:
    from . import _gfs_dynamics
except ImportError as error:
    logging.warning(
        'Import failed. GFS dynamical core is likely not compiled and will not '
        'be available.'
    )
    print(error)


class GFSError(Exception):
    pass


def get_valid_properties(gfs_properties, prognostic_properties, property_type):
    return_dict = {}
    for name, properties in prognostic_properties.items():
        if name not in gfs_properties:
            return_dict[name] = properties
        elif 'dims' in prognostic_properties.keys():
            extra_dims = set(prognostic_properties['dims']).difference(['*'] + list(gfs_properties[name]['dims']))
            if len(extra_dims) != 0:
                raise GFSError(
                    'Cannot handle TendencyComponent with {} {} '
                    'that has extra dimensions {} not used by GFS'.format(
                        property_type, name, extra_dims)
                )
    return return_dict


class GFSDynamicalCore(TendencyStepper):
    """
    Climt interface to the GFS dynamical core. The GFS
    code is available on `github`_.

    .. _github:
       https://github.com/jswhit/gfs-dycore
    """

    _gfs_input_properties = {
        'latitude': {
            'units': 'degrees_N',
            'dims': ['lat', 'lon'],
        },
        'longitude': {
            'units': 'degrees_E',
            'dims': ['lat', 'lon'],
        },
        'air_temperature': {
            'units': 'degK',
            'dims': ['mid_levels', 'lat', 'lon'],
        },
        'atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels': {
            'units': 'dimensionless',
            'dims': ['interface_levels'],
            'alias': 'a_coord',
        },
        'atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels': {
            'units': 'dimensionless',
            'dims': ['interface_levels'],
            'alias': 'b_coord',
        },
        'air_pressure': {
            'dims': ['mid_levels', 'lat', 'lon'],
            'units': 'Pa'
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', 'lat', 'lon'],
            'units': 'Pa'
        },
        'surface_air_pressure': {
            'units': 'Pa',
            'dims': ['lat', 'lon'],
        },
        'eastward_wind': {
            'units': 'm s^-1',
            'dims': ['mid_levels', 'lat', 'lon'],
        },
        'northward_wind': {
            'units': 'm s^-1',
            'dims': ['mid_levels', 'lat', 'lon'],
        },
        'divergence_of_wind': {
            'units': 's^-1',
            'dims': ['mid_levels', 'lat', 'lon'],
        },
        'atmosphere_relative_vorticity': {
            'units': 's^-1',
            'dims': ['mid_levels', 'lat', 'lon'],
        },
        'surface_geopotential': {
            'units': 'm^2 s^-2',
            'dims': ['lat', 'lon'],
        },
    }

    _gfs_output_properties = {
        'air_temperature': {'units': 'degK'},
        'air_pressure': {
            'dims': ['mid_levels', 'lat', 'lon'],
            'units': 'Pa'
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', 'lat', 'lon'],
            'units': 'Pa'
        },
        'surface_air_pressure': {'units': 'Pa'},
        'eastward_wind': {'units': 'm s^-1'},
        'northward_wind': {'units': 'm s^-1'},
        'divergence_of_wind': {'units': 's^-1'},
        'atmosphere_relative_vorticity': {'units': 's^-1'},
    }

    _gfs_diagnostic_properties = {}

    input_properties = None
    output_properties = None
    diagnostic_properties = None

    uses_tracers = True
    tracer_dims = ('tracer', 'mid_levels', 'lat', 'lon')
    prepend_tracers = (('specific_humidity', 'kg/kg'),)

    @property
    def spectral_names(self):
        return (
            'eastward_wind', 'northward_wind', 'air_temperature',
            'surface_air_pressure') + get_tracer_names()

    def __init__(
            self,
            tendency_component_list=None,
            number_of_damped_levels=0,
            damping_timescale=2.*86400,
            moist=True,
            zero_negative_moisture=True,
    ):
        """
        Initialise the GFS dynamical core.

        Args:
            tendency_component_list (list of TendencyComponent, optional):
                The TendencyComponent objects to use for spectral time stepping.

            number_of_damped_levels (int, optional):
                The number of levels from the model top which are Rayleigh damped.

            damping_timescale (float, optional):
                The damping timescale in :math:`s` to use for top-of-model Rayleigh damping.

            zero_negative_moisture (bool, optional):
                If True, all negative values of moisture will be set to zero
                before tracers are returned.
        """
        tendency_component_list = tendency_component_list or []
        self._tendency_component = ImplicitTendencyComponentComposite(*tendency_component_list)
        bad_diagnostics = set(
            self._tendency_component.diagnostic_properties.keys()).intersection(
            self.spectral_names)
        if len(bad_diagnostics) > 0:
            raise GFSError(
                'Cannot use TendencyComponent components that produce {} as diagnostic '
                'output as these must only be stepped spectrally.'.format(
                    bad_diagnostics,
                )
            )
        bad_diagnostics = set(
            self._tendency_component.diagnostic_properties.keys()).intersection(
            self._gfs_diagnostic_properties.keys())
        if len(bad_diagnostics) > 0:
            raise GFSError(
                'Cannot use TendencyComponent components that produce {} as diagnostic'
                ' output as these are already diagnosed by GFS.'.format(
                    bad_diagnostics,
                )
            )

        self._update_constants()

        self._damping_levels = number_of_damped_levels
        self._tau_damping = damping_timescale
        self._zero_negative_moisture = zero_negative_moisture

        self.initialized = False

        self.input_properties = self._gfs_input_properties.copy()
        self.output_properties = self._gfs_output_properties.copy()
        self.diagnostic_properties = self._gfs_diagnostic_properties.copy()
        super(GFSDynamicalCore, self).__init__()
        self.input_properties.update(
            get_valid_properties(
                self._gfs_input_properties, self._tendency_component.input_properties, 'input')
        )
        self.output_properties.update(
            get_valid_properties(
                self._gfs_output_properties, self._tendency_component.tendency_properties, 'output')
        )
        self.diagnostic_properties.update(
            get_valid_properties(
                self._gfs_diagnostic_properties, self._tendency_component.diagnostic_properties, 'diagnostic')
        )

    def _update_constants(self):
        self._radius = get_constant('planetary_radius', 'm')
        self._omega = get_constant('planetary_rotation_rate', 's^-1')
        self._R = get_constant('universal_gas_constant', 'J/mole/K')
        self._Rd = get_constant('gas_constant_of_dry_air', 'J/kg/K')
        self._Rv = get_constant('gas_constant_of_vapor_phase', 'J/kg/K')
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._Cp = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K')
        self._Cvap = get_constant('heat_capacity_of_vapor_phase', 'J/kg/K')
        self._fvirt = (1 - self._Rd/self._Rv)/(self._Rd/self._Rv)
        self._dry_pressure = get_constant('reference_air_pressure', 'Pa')
        self._toa_pressure = get_constant('top_of_model_pressure', 'Pa')

    def _initialize_model(self, state, timestep):
        assert not self.initialized
        self._time_step = timestep
        self._num_tracers = state['tracers'].shape[0]
        self._num_levs, self._num_lats, self._num_lons = state['air_temperature'].shape

        self._truncation = int(self._num_lons/3 - 2)

        self._spectral_dim = int(
            (self._truncation + 1)*(self._truncation + 2)/2)

        _gfs_dynamics.set_time_step(
            self._time_step.total_seconds())

        _gfs_dynamics.set_constants(
            self._radius, self._omega, self._R, self._Rd, self._Rv, self._g,
            self._Cp, self._Cvap)

        _gfs_dynamics.set_model_grid(
            self._num_lats,
            self._num_lons,
            self._num_levs,
            self._truncation,
            self._spectral_dim,
            self._num_tracers,
            state['a_coord'][::-1],
            state['b_coord'][::-1],
        )

        logging.info('Initialising dynamical core, this could take some time...')

        gaussian_weights, area_weights, latitudes, longitudes = \
            _gfs_dynamics.init_model(
                self._dry_pressure,
                self._damping_levels,
                self._tau_damping,
                self._toa_pressure)

        np.testing.assert_allclose(latitudes*180./np.pi, state['latitude'])
        np.testing.assert_allclose(longitudes*180./np.pi, state['longitude'])

        logging.info('Done!')

        self._hash_u = 1000
        self._hash_v = 1000
        self._hash_temperature = 1000
        self._hash_surface_pressure = 1000
        self._hash_tracers = 1000

        self.initialized = True

    def __call__(self, state, timestep):
        """
        Gets diagnostics from the current model state and steps the state
        forward in time according to the timestep.

        Args
        ----
        state : dict
            A model state dictionary satisfying the input_properties of this
            object.
        timestep : timedelta
            The amount of time to step forward.

        Returns
        -------
        diagnostics : dict
            Diagnostics from the timestep of the input state.
        new_state : dict
            A dictionary whose keys are strings indicating
            state quantities and values are the value of those quantities
            at the timestep after input state.

        Raises
        ------
        KeyError
            If a required quantity is missing from the state.
        InvalidStateError
            If state is not a value input for the Stepper instance
            for other reasons.
        """
        self._check_self_is_initialized()
        self._input_checker.check_inputs(state)
        raw_state = get_numpy_arrays_with_properties(state, self.input_properties)
        if self.uses_tracers:
            raw_state['tracers'] = self._tracer_packer.pack(state)
        raw_state['time'] = state['time']

        tendencies, diagnostics = self._tendency_component(state, timestep)
        for name, value in tendencies.items():
            if name in self.input_properties.keys():
                tendencies[name] = value.to_units(self.input_properties[name]['units'] + ' s^-1')
        raw_diagnostics, raw_new_state = self.array_call(
            raw_state, timestep, prognostic_tendencies=tendencies)
        if self.uses_tracers:
            new_state = self._tracer_packer.unpack(raw_new_state.pop('tracers'), state)
        else:
            new_state = {}
        if self.tendencies_in_diagnostics:
            self._insert_tendencies_to_diagnostics(
                raw_state, raw_new_state, timestep, raw_diagnostics)

        diagnostics.update(restore_data_arrays_with_properties(
            raw_diagnostics, self.diagnostic_properties,
            state, self.input_properties, ignore_missing=True))
        new_state.update(restore_data_arrays_with_properties(
            raw_new_state, self.output_properties,
            state, self.input_properties, ignore_missing=True))

        gfs_output_quantities = list(self._gfs_output_properties.keys())
        for tracer in self.prepend_tracers:
            gfs_output_quantities.append(tracer[0])
        remaining = set(tendencies.keys()).difference(gfs_output_quantities)

        for name in remaining:
            new_state[name] = state[name] + tendencies[name]*timestep.total_seconds()
        for key in state.keys():
            if key not in new_state:
                new_state[key] = state[key]

        check_new_state = {
            name: quantity
            for (name, quantity) in new_state.items()
            if name in self.output_properties.keys()
        }
        self._diagnostic_checker.check_diagnostics(diagnostics)
        self._output_checker.check_outputs(check_new_state)
        return diagnostics, new_state

    @ensure_contiguous_state
    def array_call(self, state, timestep, prognostic_tendencies=None):
        """ Step the dynamical core by one step

        Args:
            state (dict): The state dictionary of numpy arrays.

        Returns:
            new_state, diagnostics (dict):
                The new state and associated diagnostics.
        """
        prognostic_tendencies = prognostic_tendencies or {}
        self._update_constants()
        nlev, nlat, nlon = state['air_temperature'].shape
        if nlat < 16:
            raise GFSError('GFS requires at least 16 latitudes.')
        if nlon < 12:
            raise GFSError('GFS requires at least 12 longitudes')
        if not self.initialized:
            self._initialize_model(state, timestep)
        if nlev != self._num_levs:
            raise GFSError(
                'Number of vertical levels may not change between successive '
                'calls to GFS. Last time was {}, this time is {}'.format(
                    self._num_levs, nlev)
            )
        if nlat != self._num_lats:
            raise GFSError(
                'Number of latitudes may not change between successive '
                'calls to GFS. Last time was {}, this time is {}'.format(
                    self._num_lats, nlat)
            )
        if nlon != self._num_lons:
            raise GFSError(
                'Number of longitudes may not change between successive '
                'calls to GFS. Last time was {}, this time is {}'.format(
                    self._num_lons, nlon)
            )
        if timestep.total_seconds() != self._time_step.total_seconds():
            raise GFSError(
                'GFSDynamicalCore can only be run with a constant timestep.'
            )

        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties,
            prepend_tracers=self.prepend_tracers,
            tracer_dims=self.tracer_dims,
        )

        lnsp = np.log(state['surface_air_pressure'])
        t_virt = state['air_temperature']*(
            1 + self._fvirt*state['tracers'][0, :, :, :])

        outputs['air_pressure_on_interface_levels'][:] = (
            state['air_pressure_on_interface_levels'][::-1, :, :])
        for name in (
                'air_pressure', 'tracers', 'eastward_wind', 'northward_wind',
                'divergence_of_wind', 'atmosphere_relative_vorticity',
                'surface_air_pressure',):
            if np.product(outputs[name].shape) > 0:
                outputs[name][:] = state[name]

        _gfs_dynamics.assign_grid_arrays(
            outputs['eastward_wind'],
            outputs['northward_wind'],
            t_virt,
            lnsp,
            outputs['tracers'],
            outputs['atmosphere_relative_vorticity'],
            outputs['divergence_of_wind'])

        _gfs_dynamics.assign_pressure_arrays(
            outputs['surface_air_pressure'],
            outputs['air_pressure'],
            outputs['air_pressure_on_interface_levels'])

        _gfs_dynamics.set_topography(state['surface_geopotential'])

        _gfs_dynamics.calculate_pressure()

        np.testing.assert_allclose(outputs['air_pressure'], state['air_pressure'])
        np.testing.assert_allclose(
            outputs['air_pressure_on_interface_levels'][::-1, :, :],
            state['air_pressure_on_interface_levels']
        )

        self._update_spectral_arrays(state)

        tendency_arrays = \
            self._get_tendency_arrays(
                prognostic_tendencies, state['air_temperature'].shape)

        # see Pg. 12 in gfsModelDoc.pdf
        virtual_temp_tend = tendency_arrays['air_temperature']*(
            1 + self._fvirt*state['tracers'][0, :, :, :]) + \
            self._fvirt*t_virt*tendency_arrays['tracers'][0, :, :, :]

        # dlnps/dt = (1/ps)*dps/dt
        lnps_tend = ((1. / state['surface_air_pressure']) *
                     tendency_arrays['surface_air_pressure'])

        _gfs_dynamics.assign_tendencies(
            tendency_arrays['eastward_wind'],
            tendency_arrays['northward_wind'],
            virtual_temp_tend,
            lnps_tend,
            tendency_arrays['tracers'])

        _gfs_dynamics.take_one_step()
        _gfs_dynamics.convert_to_grid()
        _gfs_dynamics.calculate_pressure()

        if self._zero_negative_moisture:
            set_negatives_to_zero(outputs['tracers'][0, :, :, :])
        outputs['air_temperature'][:] = t_virt/(
            1 + self._fvirt*outputs['tracers'][0, :, :, :])

        outputs['air_pressure_on_interface_levels'][:] = \
            outputs['air_pressure_on_interface_levels'][::-1, :, :]

        outputs['time'] = state['time']

        return {}, outputs

    def _get_tendency_arrays(self, tendencies, T_shape):
        out_arrays = {}
        out_arrays.update(self._get_3d_tendencies(tendencies, T_shape))
        out_arrays.update(self._get_2d_tendencies(tendencies, T_shape))
        out_arrays['tracers'] = self._get_tracer_tendencies(tendencies, T_shape)
        return out_arrays

    def _get_3d_tendencies(self, tendencies, T_shape):
        return self._get_tendencies(
            tendencies,
            ['air_temperature', 'eastward_wind', 'northward_wind'],
            ['mid_levels', 'lat', 'lon'],
            T_shape
        )

    def _get_2d_tendencies(self, tendencies, T_shape):
        return self._get_tendencies(
            tendencies,
            ['surface_air_pressure'],
            ['lat', 'lon'],
            T_shape[1:]
        )

    def _get_tendencies(self, tendencies, quantity_names, dims, shape):
        out_arrays = {}
        for name in quantity_names:
            if name in tendencies:
                property_dict = {
                    name: {
                        'dims': dims,
                        'units': tendencies[name].attrs['units'],
                    }
                }
                out_arrays.update(
                    get_numpy_arrays_with_properties(tendencies, property_dict))
            else:
                out_arrays[name] = np.zeros(shape)
        # must broadcast any singleton dimensions
        for name, current_array in out_arrays.items():
            if current_array.shape != shape:
                new_array = np.empty(shape)
                new_array[:] = current_array
                out_arrays[name] = new_array
        return out_arrays

    def _get_tracer_tendencies(self, tendencies, T_shape):
        return_array = np.empty([self._num_tracers] + list(T_shape))
        for i, name in enumerate(self._tracer_packer.tracer_names):
            if name in tendencies:
                property_dict = {
                    name: {
                        'dims': ['mid_levels', 'lat', 'lon'],
                        'units': tendencies[name].attrs['units'],
                    }
                }
                tend = get_numpy_arrays_with_properties(
                    tendencies, property_dict)[name]
                return_array[i, :, :, :] = tend
            else:
                return_array[i, :, :, :] = 0.
        return return_array

    def _update_spectral_arrays(self, state):
        """
        Checks the state to see if any arrays with spectral counterparts have
        been modified since they were last returned, and if they have will
        update the specctral counterpart to reflect the new state array.
        """
        u_hash = get_hash(state['eastward_wind'])
        v_hash = get_hash(state['northward_wind'])
        if (u_hash != self._hash_u or v_hash != self._hash_v):
            _gfs_dynamics.vrt_div_to_spectral()
            self._hash_u = u_hash
            self._hash_v = v_hash
        T_hash = get_hash(state['air_temperature'])
        if T_hash != self._hash_temperature:
            _gfs_dynamics.virtemp_to_spectral()
            self._hash_temperature = T_hash
        tracer_hash = get_hash(state['tracers'])
        if tracer_hash != self._hash_tracers:
            _gfs_dynamics.tracer_to_spectral()
            self._hash_tracers = tracer_hash
        p_surf_hash = get_hash(state['surface_air_pressure'])
        if p_surf_hash != self._hash_surface_pressure:
            _gfs_dynamics.lnps_to_spectral()
            self._hash_surface_pressure = p_surf_hash

    def __del__(self):
        """ call shutdown in fortran code """
        logging.info("Cleaning up dynamical core...")
        _gfs_dynamics.shut_down_model()
        logging.info("Done!")


def get_hash(array):
    if sys.version_info > (3, 0):
        return hash(array.data.tobytes())
    else:
        array.flags.writeable = False
        return hash(array.data)


def set_negatives_to_zero(array):
    array[array < 0] = 0
    return array


'''
def return_tendency_arrays_or_zeros(quantity_list, state, tendencies):

    tendency_list = []
    for quantity in quantity_list:
        if quantity in tendencies.keys():
            tendency_list.append(np.asfortranarray(tendencies[quantity].values))
        elif quantity in state.keys():
            tendency_list.append(np.zeros(state[quantity].shape, order='F'))
        else:
            raise IndexError("{} not found in input state or tendencies".format(quantity))

    return tendency_list
'''

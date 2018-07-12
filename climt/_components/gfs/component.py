from __future__ import division
from ..._core import numpy_version_of, ensure_contiguous_state
from sympl import (
    DataArray, get_constant, Implicit, initialize_numpy_arrays_with_properties,
    get_numpy_arrays_with_properties, AdamsBashforth, PrognosticComposite,
    Prognostic, get_tracer_names, restore_data_arrays_with_properties,

)
import numpy as np
import sys
from datetime import timedelta
import logging
try:
    from . import _gfs_dynamics
except ImportError:
    logging.warning(
        'Import failed. GFS dynamical core is likely not compiled and will not '
        'be available.'
    )


class GFSError(Exception):
    pass


class SelectivePrognostic(Prognostic):

    @property
    def input_properties(self):
        return self.prognostic.input_properties

    @property
    def tendency_properties(self):
        return_dict = {}
        return_dict.update(self.prognostic.tendency_properties)
        return self._filter_tendency_dict(return_dict)

    def _filter_tendency_dict(self, tendency_dict):
        return_dict = {}
        if self.ignore_names is not None:
            for name in set(tendency_dict.keys()).difference(self.ignore_names):
                return_dict[name] = tendency_dict[name]
        elif self.include_names is not None:
            for name in set(tendency_dict.keys()).intersection(self.include_names):
                return_dict[name] = tendency_dict[name]
        return return_dict

    @property
    def diagnostic_properties(self):
        return self.prognostic.diagnostic_properties

    def __init__(self, prognostic, include_names=None, ignore_names=None):
        if not (include_names is None or ignore_names is None):
            raise ValueError('Cannot give both include_names and ignore_names')
        self.prognostic = prognostic
        self.ignore_names = ignore_names
        self.include_names = include_names
        super(SelectivePrognostic, self).__init__()

    def __call__(self, *args, **kwargs):
        tendencies, diagnostics = self.prognostic(*args, **kwargs)
        return self._filter_tendency_dict(tendencies), diagnostics

    def array_call(self, state):
        tendencies, diagnostics = self.prognostic.array_call(state)
        return self._filter_tendency_dict(tendencies), diagnostics


class GFSDynamicalCore(Implicit):
    """
    Climt interface to the GFS dynamical core. The GFS
    code is available on `github`_.

    .. _github:
       https://github.com/jswhit/gfs-dycore
    """

    input_properties = {
        'latitude': {
            'units': 'degrees_N',
            'dims': ['latitude'],
        },
        'longitude': {
            'units': 'degrees_E',
            'dims': ['longitude'],
        },
        'air_temperature': {
            'units': 'degK',
            'dims': ['mid_levels', 'latitude', 'longitude'],
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
            'dims': ['mid_levels', 'latitude', 'longitude'],
            'units': 'Pa'
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', 'latitude', 'longitude'],
            'units': 'Pa'
        },
        'surface_air_pressure': {
            'units': 'Pa',
            'dims': ['latitude', 'longitude'],
        },
        'eastward_wind': {
            'units': 'm s^-1',
            'dims': ['mid_levels', 'latitude', 'longitude'],
        },
        'northward_wind': {
            'units': 'm s^-1',
            'dims': ['mid_levels', 'latitude', 'longitude'],
        },
        'divergence_of_wind': {
            'units': 's^-1',
            'dims': ['mid_levels', 'latitude', 'longitude'],
        },
        'atmosphere_relative_vorticity': {
            'units': 's^-1',
            'dims': ['mid_levels', 'latitude', 'longitude'],
        },
        'surface_geopotential': {
            'units': 'm^2 s^-2',
            'dims': ['latitude', 'longitude'],
        },
    }

    output_properties = {
        'air_temperature': {'units': 'degK'},
        'air_pressure': {
            'dims': ['mid_levels', 'latitude', 'longitude'],
            'units': 'Pa'
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', 'latitude', 'longitude'],
            'units': 'Pa'
        },
        'surface_air_pressure': {'units': 'Pa'},
        'eastward_wind': {'units': 'm s^-1'},
        'northward_wind': {'units': 'm s^-1'},
        'divergence_of_wind': {'units': 's^-1'},
        'atmosphere_relative_vorticity': {'units': 's^-1'},
    }

    diagnostic_properties = {}

    uses_tracers = True
    tracer_dims = ['tracer', 'mid_levels', 'latitude', 'longitude']

    @property
    def spectral_names(self):
        return (
            'eastward_wind', 'northward_wind', 'air_temperature',
            'surface_air_pressure') + get_tracer_names()

    def __init__(
            self,
            prognostic_list=None,
            timestepper=None,
            number_of_damped_levels=0,
            damping_timescale=2.*86400,
            moist=True
    ):
        """
        Initialise the GFS dynamical core.

        Args:
            number_of_damped_levels (int, optional):
                The number of levels from the model top which are Rayleigh damped.

            damping_timescale (float, optional):
                The damping timescale in :math:`s` to use for top-of-model Rayleigh damping.

            moist (bool, optional):
                Whether to account for and advect moisture. Default is True.

        """
        prognostic_list = prognostic_list or []
        self._prognostic = PrognosticComposite(*prognostic_list)
        self._prognostic_timestepper = timestepper or AdamsBashforth(
            SelectivePrognostic(
                self._prognostic,
                ignore_names=self.spectral_names,
            ),
            order=1,
        )
        bad_diagnostics = set(
            self._prognostic.diagnostic_properties.keys()).intersection(
            self.spectral_names)
        if len(bad_diagnostics) > 0:
            raise GFSError(
                'Cannot use Prognostic components that produce {} as diagnostic '
                'output as these must only be stepped spectrally.'.format(
                    bad_diagnostics,
                )
            )
        bad_diagnostics = set(
            self._prognostic.diagnostic_properties.keys()).intersection(
            self.diagnostic_properties.keys())
        if len(bad_diagnostics) > 0:
            raise GFSError(
                'Cannot use Prognostic components that produce {} as diagnostic'
                ' output as these are already diagnosed by GFS.'.format(
                    bad_diagnostics,
                )
            )

        self._update_constants()

        self._damping_levels = number_of_damped_levels
        self._tau_damping = damping_timescale

        if moist:
            self.prepend_tracers = [('specific_humidity', 'kg/kg')]
            self.moist = True
        else:
            self.prepend_tracers = []
            self.moist = False

        self.initialized = False
        super(GFSDynamicalCore, self).__init__()

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

    def _initialize_model(self, state, timestep):
        assert not self.initialized
        self._time_step = timestep
        self._num_tracers = state['tracers'].shape[0]
        if self._num_tracers > 0 and not self.moist:
            raise GFSError(
                'GFS does not support tracers when running as a dry model. '
                'It would assume the first tracer is specific humidity.')
        assert not self._num_tracers == 0 and self.moist
        self._num_levs, self._num_lats, self._num_lons = state['air_temperature'].shape

        if self._num_levs != 28:
            raise NotImplementedError(
                'GFS can currently only run with 28 vertical levels.')

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
            state['a_coord'],
            state['b_coord'],
        )

        logging.info('Initialising dynamical core, this could take some time...')

        gaussian_weights, area_weights, latitudes, longitudes = \
            _gfs_dynamics.init_model(
                self._dry_pressure,
                self._damping_levels,
                self._tau_damping)

        np.testing.assert_allclose(latitudes[:, 0]*180./np.pi, state['latitude'])
        np.testing.assert_allclose(longitudes[0, :]*180./np.pi, state['longitude'])

        logging.info('Done!')

        # Random array to slice variables
        self.initialise_state_signature()
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
            If state is not a valid input for the Implicit instance
            for other reasons.
        """
        self._check_self_is_initialized()
        self._input_checker.check_inputs(state)
        raw_state = get_numpy_arrays_with_properties(state, self.input_properties)
        if self.uses_tracers:
            raw_state['tracers'] = self._tracer_packer.pack(raw_state)
            for name in self._tracer_packer.tracer_names:
                raw_state.pop(name)
        raw_state['time'] = state['time']
        raw_diagnostics, raw_new_state = self.array_call(raw_state, timestep)
        if self.uses_tracers:
            raw_new_state.update(self._tracer_packer.unpack(raw_new_state['tracers']))
            raw_new_state.pop('tracers')
        if self.tendencies_in_diagnostics:
            self._insert_tendencies_to_diagnostics(
                raw_state, raw_new_state, timestep, raw_diagnostics)
        diagnostics = restore_data_arrays_with_properties(
            raw_diagnostics, self.diagnostic_properties,
            state, self.input_properties)
        new_state = restore_data_arrays_with_properties(
            raw_new_state, self.output_properties,
            state, self.input_properties)
        prog_diagnostics, prog_new_state = self._prognostic_timestepper(new_state, timestep)
        diagnostics.update(prog_diagnostics)
        self._diagnostic_checker.check_diagnostics(diagnostics)
        self._output_checker.check_outputs(prog_new_state)
        return diagnostics, prog_new_state

    @ensure_contiguous_state
    def array_call(self, state, timestep):
        """ Step the dynamical core by one step

        Args:
            state (dict): The state dictionary of numpy arrays.

        Returns:
            new_state, diagnostics (dict):
                The new state and associated diagnostics.
        """
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
            tracer_dims=self.tracer_dims
        )


        lnsp = np.log(state['surface_air_pressure'])
        if self.moist:
            t_virt = state['air_temperature']*(
                1 + self._fvirt*state['tracers'][0, :, :, :])
        else:
            t_virt = state['air_temperature'].copy()

        outputs['air_pressure_on_interface_levels'][:] = (
            state['air_pressure_on_interface_levels'][::-1, :, :])
        for name in (
                'air_pressure', 'tracers', 'eastward_wind', 'northward_wind',
                'divergence_of_wind', 'atmosphere_relative_vorticity',
                'surface_air_pressure',
            ):
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

        tendencies, _ = self._prognostic(state)

        temp_tend, u_tend, v_tend, ps_tend = \
            return_tendency_arrays_or_zeros(
                ['air_temperature',
                 'eastward_wind',
                 'northward_wind',
                 'surface_air_pressure'],
                state, tendencies)
        tracer_tend = self._tracer_packer.pack(tendencies)

        # see Pg. 12 in gfsModelDoc.pdf
        if self.moist:
            virtual_temp_tend = temp_tend*(
                1 + self._fvirt*state['tracers'][0, :, :, :]) + \
                self._fvirt*t_virt*tracer_tend[0, :, :, :]
        else:
            virtual_temp_tend = temp_tend

        # dlnps/dt = (1/ps)*dps/dt
        lnps_tend = (1. / state['surface_air_pressure'])*ps_tend

        _gfs_dynamics.assign_tendencies(u_tend, v_tend, virtual_temp_tend,
                                        lnps_tend, tracer_tend)

        _gfs_dynamics.take_one_step()
        _gfs_dynamics.convert_to_grid()
        _gfs_dynamics.calculate_pressure()

        if self.moist:
            # set_negatives_to_zero(outputs['tracers'][0, :, :, :])
            outputs['air_temperature'][:] = t_virt/(
                1 + self._fvirt*outputs['specific_humidity'])
        else:
            outputs['air_temperature'][:] = t_virt

        outputs['air_pressure_on_interface_levels'][:] = \
            outputs['air_pressure_on_interface_levels'][::-1, :, :]

        for quantity in tendencies.keys():
            if quantity not in self._climt_outputs.keys():
                # Step forward using Euler
                outputs[quantity] = state[quantity] + \
                    tendencies[quantity].values*self._time_step.total_seconds()

                for attrib in state[quantity].attrs:
                    outputs[quantity].attrs[attrib] = state[quantity].attrs[attrib]

        return {}, outputs

    def _update_spectral_arrays(self, state):
        """
        Checks the state to see if any arrays with spectral counterparts have
        been modified since they were last returned, and if they have will
        update the specctral counterpart to reflect the new state array.
        """
        u_hash = get_hash(state['eastward_wind'])
        v_hash = get_hash(state['northward_wind'])
        if (u_hash != self._hash_u or
            v_hash != self._hash_v):
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

    def initialise_state_signature(self):
        self._hash_u = 1000
        self._hash_v = 1000
        self._hash_temperature = 1000
        self._hash_surface_pressure = 1000
        self._hash_tracers = 1000

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

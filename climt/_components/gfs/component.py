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
        'air_pressure': {
            'units': 'Pa',
            'dims': ['mid_levels', 'latitude', 'longitude'],
        },
        'air_pressure_on_interface_levels': {
            'units': 'Pa',
            'dims': ['interface_levels', 'latitude', 'longitude'],
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
        'specific_humidity': {
            'units': 'g/g',
            'dims': ['mid_levels', 'latitude', 'longitude'],
        },
        'surface_geopotential': {
            'units': 'm^2 s^-2',
            'dims': ['latitude', 'longitude'],
        },
    }

    output_properties = {
        'air_temperature': {'units': 'degK'},
        'air_pressure': {'units': 'Pa'},
        'air_pressure_on_interface_levels': {'units': 'Pa'},
        'surface_air_pressure': {'units': 'Pa'},
        'specific_humidity': {'units': 'g/g'},
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
            self._num_tracers)

        logging.info('Initialising dynamical core, this could take some time...')

        gaussian_weights, area_weights, latitudes, longitudes, sigma, sigma_interface = \
            _gfs_dynamics.init_model(
                self._dry_pressure,
                self._damping_levels,
                self._tau_damping)

        np.testing.assert_allclose(latitudes[:, 0]*180./np.pi, state['latitude'])
        np.testing.assert_allclose(longitudes[0, :]*180./np.pi, state['longitude'])
        sigma_input = state['air_pressure'] / state['surface_air_pressure'][None, :, :]
        assert np.all(np.var(sigma_input, axis=(1, 2)) < 1e-7)  # constant sigma levels
        sigma_interface_input = state['air_pressure'] / state['surface_air_pressure'][None, :, :]
        assert np.all(np.var(sigma_interface_input, axis=(1, 2)) < 1e-7)  # constant sigma levels
        np.testing.assert_allclose(sigma_input[:, 0, 0], sigma)
        np.testing.assert_allclose(sigma_interface_input[:, 0, 0], sigma_interface)

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
        print(list(state.keys()))
        self._update_constants()
        if not self.initialized:
            self._initialize_model(state, timestep)
        if timestep.total_seconds() != self._time_step.total_seconds():
            raise GFSError(
                'GFSDynamicalCore can only be run with a constant timestep.'
            )
        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties,
            tracer_dims=self.tracer_dims
        )

        self._update_spectral_arrays(state)

        lnsp = np.log(state['surface_air_pressure'])
        t_virt = state['air_temperature']*(
            1 + self._fvirt*state['specific_humidity'])

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

        tendencies, _ = self._prognostic(state)

        temp_tend, q_tend, u_tend, v_tend, ps_tend = \
            return_tendency_arrays_or_zeros(
                ['air_temperature',
                 'specific_humidity',
                 'eastward_wind',
                 'northward_wind',
                 'surface_air_pressure'],
                state, tendencies)
        tracer_tend = self._tracer_packer.pack(tendencies)

        # see Pg. 12 in gfsModelDoc.pdf
        virtual_temp_tend = temp_tend*(
            1 + self._fvirt*state['specific_humidity']) + \
            self._fvirt*t_virt*q_tend

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

        self.store_current_state_signature(outputs)

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
        raise NotImplementedError()

    def initialise_state_signature(self):

        self._random_slice_x = np.random.randint(0, self._num_lons, size=(10, 10, 10))
        self._random_slice_y = np.random.randint(0, self._num_lats, size=(10, 10, 10))
        self._random_slice_z = np.random.randint(0, self._num_levs, size=(10, 10, 10))

        self._hash_u = 1000
        self._hash_v = 1000
        self._hash_temp = 1000
        self._hash_press = 1000
        self._hash_surf_press = 1000

    def calculate_state_signature(self, state):
        """ Calculates hash signatures from state """
        random_u = state['eastward_wind'][
            self._random_slice_z, self._random_slice_y, self._random_slice_x]
        hash_u = get_hash(random_u)
        random_v = state['northward_wind'][
            self._random_slice_z, self._random_slice_y, self._random_slice_x]
        hash_v = get_hash(random_v)
        random_temperature = state['air_temperature'][
            self._random_slice_z, self._random_slice_y, self._random_slice_x]
        hash_temperature = get_hash(random_temperature)
        random_pressure = state['air_pressure'][
            self._random_slice_z, self._random_slice_y, self._random_slice_x]
        hash_pressure = get_hash(random_pressure)
        random_ps = state['surface_air_pressure'][
            self._random_slice_y, self._random_slice_x]
        hash_ps = get_hash(random_ps)
        return hash_u, hash_v, hash_temperature, hash_pressure, hash_ps

    def state_is_modified_externally(self, state_arr):
        """ Function to check if grid space arrays have been modified outside the dynamical core """

        hash_u, hash_v, hash_temp, hash_press, hash_ps = self.calculate_state_signature(state_arr)

        if (
            (hash_u != self._hash_u) or
            (hash_v != self._hash_v) or
            (hash_press != self._hash_press) or
            (hash_ps != self._hash_surf_press) or
           (hash_temp != self._hash_temp)):
                logging.info('State modified, setting spectral arrays')
                self._hash_u = hash_u
                self._hash_v = hash_v
                self._hash_temp = hash_temp
                self._hash_surf_press = hash_ps
                self._hash_press = hash_press
                return True
        else:
            return False

    def store_current_state_signature(self, output_arr):
        """ Store state signature for comparison during next time step """

        hash_u, hash_v, hash_temp, hash_press, hash_ps = self.calculate_state_signature(output_arr)

        self._hash_u = hash_u
        self._hash_v = hash_v
        self._hash_temp = hash_temp
        self._hash_surf_press = hash_ps
        self._hash_press = hash_press

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

from __future__ import division
from ..._core import ClimtSpectralDynamicalCore, numpy_version_of, get_constant
from sympl import DataArray
import numpy as np
import sys
from datetime import timedelta
try:
    from . import _gfs_dynamics
except ImportError:
    print("Import failed. GFS dynamical core will not be available!")


class GFSDynamicalCore(ClimtSpectralDynamicalCore):
    """
    Climt interface to the GFS dynamical core. The GFS
    code is available on `github`_.

    :attribute area_weights: Weights used for calculating area averages.
    :attribute gauss_weights: Unnormalised weights used for calculating integrals along longitude.

    .. _github:
       https://github.com/jswhit/gfs-dycore
    """

    _climt_inputs = {
        'eastward_wind': 'm s^-1',
        'northward_wind': 'm s^-1',
        'air_temperature': 'degK',
        'surface_air_pressure': 'Pa',
        'air_pressure': 'Pa',
        'air_pressure_on_interface_levels': 'Pa',
        'specific_humidity': 'g/g',
        'surface_geopotential': 'm^2 s^-2',
        'atmosphere_relative_vorticity': 's^-1',
        'divergence_of_wind': 's^-1',
        'mole_fraction_of_ozone_in_air': 'dimensionless',
        'mass_content_of_cloud_ice_in_atmosphere_layer': 'g m^-2',
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': 'g m^-2',
        'gfs_tracers': 'dimensionless',
    }

    _climt_outputs = {
        'eastward_wind': 'm s^-1',
        'northward_wind': 'm s^-1',
        'air_temperature': 'degK',
        'air_pressure': 'Pa',
        'air_pressure_on_interface_levels': 'Pa',
        'surface_air_pressure': 'Pa',
        'specific_humidity': 'g/g',
        'atmosphere_relative_vorticity': 's^-1',
        'divergence_of_wind': 's^-1',
        'mole_fraction_of_ozone_in_air': 'dimensionless',
        'mass_content_of_cloud_ice_in_atmosphere_layer': 'g m^-2',
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': 'g m^-2',
        'gfs_tracers': 'dimensionless',
    }

    _climt_diagnostics = {
        # 'downward_air_velocity': 'm s^-1',
    }

    extra_dimensions = {'tracer_number': np.arange(4)}

    quantity_descriptions = {
        'gfs_tracers': {
            'dims': ['x', 'y', 'mid_levels', 'tracer_number'],
            'units': 'dimensionless',
            'default_value': 0.
        }
    }

    def __init__(
            self,
            number_of_latitudes=94,
            number_of_longitudes=198,
            number_of_levels=28,
            number_of_tracers=0,
            number_of_damped_levels=0,
            damping_timescale=2.*86400,
            time_step=1200.):
        """
        Initialise the GFS dynamical core.

        Args:

            number_of_latitudes (int, optional):
                The desired number of latitudes for the model. Note that
                not all combinations of latitudes and longitudes are
                acceptable. In particular, the number of latitudes must be
                :math:`\leq (longitudes)/2`.

            number_of_longitudes (int, optional):
                The desired number of longitudes. The resolution of the model in `Txx`
                notation is approximately :math:`xx = longitudes/3`. So, 192
                longitudes is T64, etc.,

            number_of_levels (int, optional):
                The desired number of levels. **Setting this option is not supported yet.**

            number_of_tracers (int, optional):
                The number of additional tracers to be used by the model. A minimum of
                four tracers are used for specific humidity, ozone and liquid and solid cloud condensate.
                This number indicates number of tracers beyond these four. These tracers
                will appear in the state dictionary in a :code:`DataArray` whose key is
                :code:`gfs_tracers` and dimensions are
                :code:`(number_of_longitudes, number_of_latitudes, number_of_levels,
                number_of_tracers)`.

            number_of_damped_levels (int, optional):
                The number of levels from the model top which are Rayleigh damped.

            damping_timescale (float, optional):
                The damping timescale in :math:`s` to use for top-of-model Rayleigh damping.

            time_step (float, optional):
                The time step to be used by the model in :math:`s`.

        """

        self._time_step = timedelta(seconds=time_step)

        self._radius = get_constant('planetary_radius', 'm')

        self._omega = get_constant('planetary_rotation_rate', 's^-1')

        self._R = get_constant('universal_gas_constant', 'J/mole/K')

        self._Rd = get_constant('gas_constant_of_dry_air', 'J/kg/K')

        self._Rv = get_constant('gas_constant_of_vapor_phase', 'J/kg/K')

        self._g = get_constant('gravitational_acceleration', 'm/s^2')

        self._Cp = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K')

        self._Cvap = get_constant('heat_capacity_of_vapor_phase', 'J/kg/K')

        self._fvirt = (1 - self._Rd/self._Rv)/(self._Rd/self._Rv)

        # Sanity Checks
        assert number_of_tracers >= 0
        assert number_of_levels > 0
        assert number_of_latitudes > 0
        assert number_of_longitudes > 0
        assert number_of_damped_levels >= 0

        self._num_lats = number_of_latitudes

        self._num_lons = number_of_longitudes

        self._num_levs = number_of_levels

        self._damping_levels = number_of_damped_levels
        self._tau_damping = damping_timescale

        # 4 tracers at least for water vapour, ozone and liquid and solid cloud condensate
        self._num_tracers = number_of_tracers + 4
        self.extra_dimensions['tracer_number'] = np.arange(self._num_tracers)

        self._dry_pressure = get_constant('reference_air_pressure', 'Pa')

        # Cannot set to new value currently.
        if self._num_levs != 28:
            raise NotImplementedError(
                'Setting levels is not supported yet!.')

        self._truncation = int(self._num_lons/3 - 2)

        self._spectral_dim = int(
            (self._truncation + 1)*(self._truncation + 2)/2)

        _gfs_dynamics.set_time_step(self._time_step.total_seconds())

        _gfs_dynamics.set_constants(self._radius, self._omega,
                                    self._R, self._Rd, self._Rv,
                                    self._g, self._Cp, self._Cvap)

        _gfs_dynamics.set_model_grid(self._num_lats,
                                     self._num_lons,
                                     self._num_levs,
                                     self._truncation,
                                     self._spectral_dim,
                                     self._num_tracers)

        print('Initialising dynamical core, this could take some time...')

        gaussian_weights, area_weights, latitudes, longitudes, sigma, sigma_interface = \
            _gfs_dynamics.init_model(self._dry_pressure,
                                     self._damping_levels,
                                     self._tau_damping)

        print('Done!')

        self.gauss_weights = DataArray(gaussian_weights,
                                       name='gauss_weights',
                                       dims=['latitude'],
                                       coords=[np.degrees(latitudes[0, :])],
                                       attrs=dict(
                                           units='dimensionless'))

        self.area_weights = DataArray(area_weights,
                                      name='area_weights',
                                      dims=['longitude', 'latitude'],
                                      coords=[np.degrees(longitudes[:, 0]),
                                              np.degrees(latitudes[0, :])],
                                      attrs=dict(
                                          units='dimensionless'))

        latitude = dict(label='latitude',
                        values=np.degrees(latitudes[0, :]),
                        units='degrees_north')

        longitude = dict(label='longitude',
                         values=np.degrees(longitudes[:, 0]),
                         units='degrees_east')

        sigma_levels = dict(label='sigma_levels',
                            values=sigma,
                            units='dimensionless')

        sigma_int_levels = dict(label='sigma_interface_levels',
                                values=sigma_interface,
                                units='dimensionless')

        self.grid_definition = dict(y=latitude, x=longitude,
                                    mid_levels=sigma_levels,
                                    interface_levels=sigma_int_levels)

        # Random array to slice variables
        self.initialise_state_signature()

    def __call__(self, state):
        """ Step the dynamical core by one step

        Args:
            state (dict): The state dictionary.

        Returns:
            new_state, diagnostics (dict):
                The new state and associated diagnostics.
        """

        raw_input_arrays = self.get_numpy_arrays_from_state('_climt_inputs', state)

        output_dict = self.create_state_dict_for('_climt_outputs', state)
        raw_output_arrays = numpy_version_of(output_dict)

        mylist = ['air_pressure']
        # mylist = self.outputs

        for quantity in mylist:
            raw_output_arrays[quantity][:] = raw_input_arrays[quantity][:]

        update_spectral_arrays = False
        if self.state_is_modified_externally(raw_input_arrays):
            for quantity in self._climt_outputs.keys():
                raw_output_arrays[quantity][:] = raw_input_arrays[quantity][:]
            update_spectral_arrays = True

        lnsp = np.log(raw_input_arrays['surface_air_pressure'])
        t_virt = raw_input_arrays['air_temperature']*(
            1 + self._fvirt*raw_input_arrays['specific_humidity'])

        raw_output_arrays['gfs_tracers'][:, :, :, 0] = raw_input_arrays['specific_humidity']
        raw_output_arrays['gfs_tracers'][:, :, :, 1] = \
            raw_input_arrays['mole_fraction_of_ozone_in_air']
        raw_output_arrays['gfs_tracers'][:, :, :, 2] = \
            raw_input_arrays['mass_content_of_cloud_liquid_water_in_atmosphere_layer']
        raw_output_arrays['gfs_tracers'][:, :, :, 3] = \
            raw_input_arrays['mass_content_of_cloud_ice_in_atmosphere_layer']
        raw_output_arrays['air_pressure_on_interface_levels'][:] = \
            np.asfortranarray(raw_input_arrays['air_pressure_on_interface_levels'][:, :, ::-1])

        _gfs_dynamics.assign_grid_arrays(
            raw_output_arrays['eastward_wind'],
            raw_output_arrays['northward_wind'],
            t_virt,
            lnsp,
            raw_output_arrays['gfs_tracers'],
            raw_output_arrays['atmosphere_relative_vorticity'],
            raw_output_arrays['divergence_of_wind'])

        _gfs_dynamics.assign_pressure_arrays(
            raw_output_arrays['surface_air_pressure'],
            raw_output_arrays['air_pressure'],
            raw_output_arrays['air_pressure_on_interface_levels'])

        _gfs_dynamics.set_topography(raw_input_arrays['surface_geopotential'])

        tendencies = {}
        diagnostics = {}
        if self.prognostics:
            tendencies, diagnostics = self.prognostics(state)

        (temp_tend, q_tend, u_tend, v_tend, ps_tend,
         ozone_tend, cloud_water_tend, cloud_ice_tend, tracer_tend) = \
            return_tendency_arrays_or_zeros(
                ['air_temperature',
                 'specific_humidity',
                 'eastward_wind',
                 'northward_wind',
                 'surface_air_pressure',
                 'mole_fraction_of_ozone_in_air',
                 'mass_content_of_cloud_liquid_water_in_atmosphere_layer',
                 'mass_content_of_cloud_ice_in_atmosphere_layer',
                 'gfs_tracers'],
                raw_input_arrays, tendencies)

        # see Pg. 12 in gfsModelDoc.pdf
        virtual_temp_tend = temp_tend*(
            1 + self._fvirt*raw_input_arrays['specific_humidity']) + \
            self._fvirt*t_virt*q_tend

        # dlnps/dt = (1/ps)*dps/dt
        lnps_tend = (1. / raw_input_arrays['surface_air_pressure'])*ps_tend

        tracer_tend[:, :, :, 0] = q_tend
        tracer_tend[:, :, :, 1] = ozone_tend
        tracer_tend[:, :, :, 2] = cloud_water_tend
        tracer_tend[:, :, :, 3] = cloud_ice_tend

        _gfs_dynamics.assign_tendencies(u_tend, v_tend, virtual_temp_tend, q_tend,
                                        lnps_tend, tracer_tend)

        if update_spectral_arrays:
            _gfs_dynamics.update_spectral_arrays()

        _gfs_dynamics.take_one_step()
        _gfs_dynamics.convert_to_grid()
        _gfs_dynamics.calculate_pressure()

        raw_output_arrays['specific_humidity'][:] = set_negatives_to_zero(
            raw_output_arrays['gfs_tracers'][:, :, :, 0])

        raw_output_arrays['air_temperature'][:] = t_virt/(
            1 + self._fvirt*raw_output_arrays['specific_humidity'])

        raw_output_arrays['air_pressure_on_interface_levels'][:] = \
            raw_output_arrays['air_pressure_on_interface_levels'][:, :, ::-1]

        raw_output_arrays['mole_fraction_of_ozone_in_air'][:] = \
            set_negatives_to_zero(raw_output_arrays['gfs_tracers'][:, :, :, 1])

        raw_output_arrays['mass_content_of_cloud_liquid_water_in_atmosphere_layer'][:] = \
            set_negatives_to_zero(raw_output_arrays['gfs_tracers'][:, :, :, 2])

        raw_output_arrays['mass_content_of_cloud_ice_in_atmosphere_layer'][:] = \
            set_negatives_to_zero(raw_output_arrays['gfs_tracers'][:, :, :, 3])

        self.store_current_state_signature(raw_output_arrays)

        for quantity in tendencies.keys():
            if quantity not in self._climt_outputs.keys():
                # Step forward using Euler
                output_dict[quantity] = state[quantity] + \
                    tendencies[quantity].values*self._time_step.total_seconds()

                for attrib in state[quantity].attrs:
                    output_dict[quantity].attrs[attrib] = state[quantity].attrs[attrib]

        for quantity in mylist:
            raw_input_arrays[quantity][:] = raw_output_arrays[quantity][:]

        return output_dict, diagnostics

    def initialise_state_signature(self):

        self._random_slice_x = np.random.randint(0, self._num_lons, size=(10, 10, 10))
        self._random_slice_y = np.random.randint(0, self._num_lats, size=(10, 10, 10))
        self._random_slice_z = np.random.randint(0, self._num_levs, size=(10, 10, 10))

        self._hash_u = 1000
        self._hash_v = 1000
        self._hash_temp = 1000
        self._hash_press = 1000
        self._hash_surf_press = 1000

    def calculate_state_signature(self, state_arr):
        """ Calculates hash signatures from state """
        random_u = state_arr['eastward_wind'][self._random_slice_x, self._random_slice_y,
                                              self._random_slice_z]

        if sys.version_info > (3, 0):
            hash_u = hash(random_u.data.tobytes())
        else:
            random_u.flags.writeable = False
            hash_u = hash(random_u.data)

        random_v = state_arr['northward_wind'][self._random_slice_x, self._random_slice_y,
                                               self._random_slice_z]
        if sys.version_info > (3, 0):
            hash_v = hash(random_v.data.tobytes())
        else:
            random_v.flags.writeable = False
            hash_v = hash(random_v.data)

        random_temp = state_arr['air_temperature'][self._random_slice_x, self._random_slice_y,
                                                   self._random_slice_z]
        if sys.version_info > (3, 0):
            hash_temp = hash(random_temp.data.tobytes())
        else:
            random_temp.flags.writeable = False
            hash_temp = hash(random_temp.data)

        random_pressure = state_arr['air_pressure'][self._random_slice_x, self._random_slice_y,
                                                    self._random_slice_z]
        if sys.version_info > (3, 0):
            hash_press = hash(random_pressure.data.tobytes())
        else:
            random_pressure.flags.writeable = False
            hash_press = hash(random_pressure.data)

        random_ps = state_arr['surface_air_pressure'][self._random_slice_x, self._random_slice_y]

        if sys.version_info > (3, 0):
            hash_ps = hash(random_ps.data.tobytes())
        else:
            random_ps.flags.writeable = False
            hash_ps = hash(random_ps.data)

        return hash_u, hash_v, hash_temp, hash_press, hash_ps

    def state_is_modified_externally(self, state_arr):
        """ Function to check if grid space arrays have been modified outside the dynamical core """

        hash_u, hash_v, hash_temp, hash_press, hash_ps = self.calculate_state_signature(state_arr)

        if (
            (hash_u != self._hash_u) or
            (hash_v != self._hash_v) or
            (hash_press != self._hash_press) or
            (hash_ps != self._hash_surf_press) or
           (hash_temp != self._hash_temp)):
                print('State modified, setting spectral arrays')
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
        print("Cleaning up dynamical core...")
        _gfs_dynamics.shut_down_model()
        print("Done!")


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

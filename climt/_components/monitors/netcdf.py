from ..._core.base_components import Monitor
from ..._core.exceptions import (
    DependencyException, InvalidStateException, IOException)
from ..._core.units import from_unit_to_another
from ..._core.array import DataArray
import xarray as xr
import os
import numpy as np
from datetime import timedelta
try:
    import netCDF4 as nc4
except ImportError:
    nc4 = None

if nc4 is None:
    # If dependency is not installed, use a dummy object that will alert the
    # user they need to install the dependency if they try to use it
    class NetCDFMonitor(Monitor):

        def __init__(self, filename):
            raise DependencyException(
                'netCDF4-python must be installed to use NetCDFMonitor')

        def store(self, state):
            pass

else:
    class NetCDFMonitor(Monitor):
        """A Monitor which caches stored states and then writes them to a
        NetCDF file when requested."""

        def __init__(
                self, filename, time_units='seconds', write_on_store=False):
            """
            Args:
                filename (str): The file to which the NetCDF file will be
                    written.
                time_units (str, optional): The units in which time will be
                    stored in the NetCDF file. Time is stored as an integer
                    number of these units.
                write_on_store (bool, optional): If True, stored changes are
                    immediately written to file. This can result in many file
                    open/close operations. Default is to write only when
                    the write() method is called directly.
            """
            self._cached_states = {}
            self._filename = filename
            self._time_units = time_units
            self._write_on_store = write_on_store

        def store(self, state):
            """
            Stores the given state in the Monitor and performs class-specific
            actions.

            Args:
                state (dict): A model state dictionary.

            Raises:
                InvalidStateException: if state is not a valid input for the
                    Diagnostic instance.
            """
            if state['time'] in self._cached_states:
                self._cached_states[state['time']].update(state)
            else:
                # This only copies the dictionary, not the arrays it contains
                state = state.copy()
                self._cached_states[state.pop('time')] = state
            if self._write_on_store:
                self.write()

        @property
        def _write_mode(self):
            if not os.path.isfile(self._filename):
                return 'w'
            else:
                return 'a'

        def _ensure_cached_states_have_same_keys(self):
            if len(self._cached_states) == 0:
                return  # trivially true
            reference_state = tuple(self._cached_states.values())[0]
            for state in self._cached_states.values():
                if state.keys() != reference_state.keys():
                    raise InvalidStateException(
                        'NetCDFMonitor was passed a different set of '
                        'quantities for different times')

        def _get_ordered_times_and_states(self):
            return zip(*sorted(self._cached_states.items(), key=lambda x: x[0]))

        def write(self):
            """Write all cached states to the NetCDF file. This will append
            to any existing NetCDF file."""
            self._ensure_cached_states_have_same_keys()
            with nc4.Dataset(self._filename, self._write_mode) as dataset:
                times, states = self._get_ordered_times_and_states()
                self._ensure_time_exists(dataset, times[0])
                it_start = dataset.dimensions['time'].size
                it_end = it_start + len(times)
                append_times_to_dataset(times, dataset, self._time_units)
                all_states = combine_states(states)
                for name, value in all_states.items():
                    ensure_variable_exists(dataset, name, value)
                    dataset.variables[name][
                        it_start:it_end, :] = value.values[:, :]
            self._cached_states = {}

        def _ensure_time_exists(self, dataset, possible_reference_time):
            ensure_dimension_exists(dataset, 'time', None)
            if 'time' not in dataset.variables:
                dataset.createVariable('time', np.int64, ('time',))
                if isinstance(possible_reference_time, timedelta):
                    dataset.variables['time'].setncattr(
                        'units', self._time_units)
                else:  # assume datetime
                    dataset.variables['time'].setncattr(
                        'units', '{} since {}'.format(
                            self._time_units, possible_reference_time))
                    dataset.variables['time'].setncattr(
                        'calendar', 'proleptic_gregorian')


def append_times_to_dataset(times, dataset, time_units):
    it_start = dataset.dimensions['time'].size
    it_end = it_start + len(times)
    if isinstance(times[0], timedelta):
        times_list = []
        for time in times:
            times_list.append(time.total_seconds())
        time_array = from_unit_to_another(
            np.array(times_list, 'seconds', time_units))
        dataset.variables['time'][it_start:it_end] = time_array[:]
    else:  # assume datetime
        dataset.variables['time'][it_start:it_end] = nc4.date2num(
            times, dataset.variables['time'].units,
            calendar='proleptic_gregorian'
        )


def combine_states(states):
    """Takes in an iterable of state dictionaries, and returns a single
    state dictionary whose values are arrays whose first dimension corresponds
    to which state."""
    return_dict = {}
    n_states = len(states)
    for name, value in states[0].items():
        return_dict[name] = DataArray(
            np.zeros((n_states,) + value.shape, dtype=value.values.dtype),
            dims=('time',) + value.dims, attrs=value.attrs)
    for i, state in enumerate(states):
        for name in state.keys():
            return_dict[name][i, :] = state[name][:]
    return return_dict


def ensure_variable_exists(dataset, name, data):
    """Dataset should be nc4.Dataset, name should be a string, and data should
    be a DataArray.

    Ensures there is a Variable in the dataset that corresponds to the given
    name and data, and creates it if not. Raises IOException if there is already
    a Variable but it is incompatible with the data."""
    if name not in dataset.variables:
        create_variable(dataset, name, data)
    else:
        ensure_variable_is_compatible(dataset.variables[name], name, data)


def create_variable(dataset, name, data):
    if isinstance(data, xr.DataArray):
        for i in range(len(data.dims)):
            try:
                if i == 0:  # time
                    ensure_dimension_exists(
                        dataset, data.dims[i], None)
                else:
                    ensure_dimension_exists(
                        dataset, data.dims[i], data.values.shape[i])
            except IOException as err:
                raise IOException(
                    'Error while creating {}: {}'.format(name, err))
        dataset.createVariable(
            name, data.values.dtype, data.dims)
        for key, value in data.attrs.items():
            dataset.variables[name].setncattr(key, value)
    else:
        raise TypeError('data must be of type DataArray')


def ensure_variable_is_compatible(variable, name, data):
    if variable.dimensions != data.dims:
        raise IOException(
            'Dimension in file is {} but on variable is {}'.format(
                variable.dimensions, data.dims))
    for key, value in data.attrs.items():
        if key not in variable.ncattrs():
            raise InvalidStateException(
                'State has attr {} for quantity {} but this is not '
                'present in the netCDF file'.format(key, name))
        elif value != variable.getncattr(key):
            raise InvalidStateException(
                'State has attr {} with value {} for quantity {} but '
                'the value in the netCDF file is {}'.format(
                    key, value, name,
                    variable.getncattr(key)))


def ensure_dimension_exists(dataset, dim_name, dim_length):
    if dim_name in dataset.dimensions:
        if dim_length is None:
            if not dataset.dimensions[dim_name].isunlimited():
                raise IOException(
                    'Dimension {} is unlimited in file but dim_length {} '
                    'is given'.format(dim_name, dim_length))
        elif dim_length != dataset.dimensions[dim_name].size:
            raise IOException(
                'Dimension {} is length {} in file but dim_length {} '
                'is given'.format(
                    dim_name, dataset.dimensions[dim_name].size, dim_length))
    else:
        dataset.createDimension(dim_name, dim_length)

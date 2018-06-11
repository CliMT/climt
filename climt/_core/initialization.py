from sympl import (DataArray, ensure_no_shared_keys, Diagnostic)
import numpy as np
import copy
from datetime import datetime
from scipy.interpolate import CubicSpline


class ConstantDefaultValue(Diagnostic):

    def __init__(
            self, input_dict, output_name, output_value, output_units):
        self.input_properties = {
            name: {
                'dims': ['*'],
                'units': units,
            } for name, units in input_dict.items()
        }
        self.diagnostic_properties = {
            output_name : {
                'dims': ['*'],
                'units': output_units,
            },
        }
        self.output_name = output_name
        self.output_value = output_value

    def array_call(self, state):
        if hasattr(self.value, 'dtype'):
            dtype = self.value.dtype
        else:
            dtype = type(self.value)
        sample_array = state.values()[0]
        out_array = np.empty(sample_array.shape, dtype=dtype)
        out_array[:] = self.value
        return {self.output_name: out_array}


class PressureFunctionDiagnostic(Diagnostic):
    """Defines a quantity as a function of pressure."""

    input_properties = {
        'air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
            'alias': 'p'
        }
    }

    diagnostic_properties = {}

    def __init__(
            self, output_name, output_function, output_units,
            mid_or_interface_levels='mid'):
        """

        Args:
            output_name : string
                Name of the output. If outputting on interface levels, do *not*
                include '_on_interface_levels' at the end of the name.
            output_function : func
                A function which takes a numpy array of pressure as its only
                argument and returns a numpy array of the output quantity with
                the same dimensions.
            output_units: string
                The units of the output quantity
            mid_or_interface_levels: 'mid' or 'interface'
                Whether this diagnostic should perform the calculation on
                mid levels or interface levels.
        """
        if mid_or_interface_levels == 'interface':
            output_name += '_on_interface_levels'
            self.input_properties = {
                'air_pressure_on_interface_levels': {
                    'dims': ['*'],
                    'units': 'Pa',
                    'alias': 'p'
                }
            }
        elif mid_or_interface_levels == 'mid':
            pass
        else:
            raise ValueError(
                "Argument mid_or_interface_levels must be 'mid' or 'interface'")
        self.diagnostic_properties = {
            output_name: {
                'dims': ['*'],
                'units': output_units
            }
        }
        self._output_function = output_function
        self._output_name = output_name
        super(PressureFunctionDiagnostic, self).__init__()

    def array_call(self, raw_state):
        return {
            self.output_name: self._output_function(raw_state['p'])
        }


def horizontal_broadcast_if_needed(output, nx=None, ny=None):
    """
    Take a tuple of numpy arrays output. If nx and ny are not None, broadcast
    each of those numpy arrays along new dimensions of the correct length
    before returning them.
    """
    output_list = list(output)
    if ny is not None:
        for i, var in enumerate(output_list):
            output_list[i] = np.repeat(var[None, :], ny, axis=0)
    if nx is not None:
        for i, var in enumerate(output_list):
            output_list[i] = np.repeat(var[None, :], nx, axis=0)
    return tuple(output_list)


def get_grid(
        nx=None, ny=None, nz=50, p_surf_in_Pa=1.0132e5, x_name='longitude',
        y_name='latitude'):
    p, p_interface, sigma, sigma_interface = horizontal_broadcast_if_needed(
        get_pressure_and_sigma_levels(nz, p_surf_in_Pa), nx, ny
    )
    horizontal_dims = []
    if nx is not None:
        horizontal_dims.append(x_name)
    if ny is not None:
        horizontal_dims.append(y_name)
    return_state = {
        'air_pressure':
            DataArray(
                p,
                dims=horizontal_dims + ['mid_levels'],
                attrs={'units': 'Pa'}),
        'air_pressure_on_interface_levels':
            DataArray(
                p_interface,
                dims=horizontal_dims + ['interface_levels'],
                attrs={'units': 'Pa'}),
        'sigma':
            DataArray(
                sigma,
                dims=horizontal_dims + ['mid_levels'],
                attrs={'units': 'Pa'}),
        'sigma_on_interface_levels':
            DataArray(
                sigma_interface,
                dims=horizontal_dims + ['interface_levels'],
                attrs={'units': 'Pa'}),
    }
    return return_state


def get_pressure_and_sigma_levels(nz, p_surf):
    sigma = np.linspace(0.998, 0.001, nz)
    sigma_interface = np.zeros([nz+1])
    sigma_interface[0] = 1.
    sigma_interface[-1] = 0.
    sigma_interface[1:-1] = 0.5 * (sigma[:-1] + sigma[1:])
    p = p_surf * sigma
    p_interface = p_surf * sigma_interface
    return p, p_interface, sigma, sigma_interface


default_values = {
    'air_temperature': {'value': 290., 'units': 'degK'},
    'surface_temperature': {'value': 300., 'units': 'degK'},
    'sea_surface_temperature': {'value': 300., 'units': 'degK'},
    'soil_surface_temperature': {'value': 300., 'units': 'degK'},
    'northward_wind': {'value': 0., 'units': 'm/s'},
    'eastward_wind': {'value': 0., 'units': 'm/s'},
    'divergence_of_wind': {'value': 0., 'units': 's^-1'},
    'atmosphere_relative_vorticity': {'value': 0., 'units': 's^-1'},
    'surface_geopotential': {'value': 0., 'units': 'm^2 s^-2'},
    'surface_longwave_emissivity': {'value': 0., 'units': 'dimensionless'},
    'specific_humidity': {'value': 0., 'units': 'dimensionless'},
    'mole_fraction_of_carbon_dioxide_in_air': {'value': 330e-6, 'units': 'dimensionless'},
    'mole_fraction_of_methane_in_air': {'value': 0., 'units': 'dimensionless'},
    'mole_fraction_of_nitrous_oxide_in_air': {'value': 0., 'units': 'dimensionless'},
    'mole_fraction_of_oxygen_in_air': {'value': 0.21, 'units': 'dimensionless'},
    'mole_fraction_of_nitrogen_in_air': {'value': 0.78, 'units': 'dimensionless'},
    'mole_fraction_of_hydrogen_in_air': {'value': 500e-9, 'units': 'dimensionless'},
    'mole_fraction_of_cfc11_in_air': {'value': 0., 'units': 'dimensionless'},
    'mole_fraction_of_cfc12_in_air': {'value': 0., 'units': 'dimensionless'},
    'mole_fraction_of_cfc22_in_air': {'value': 0., 'units': 'dimensionless'},
    'mole_fraction_of_carbon_tetrachloride_in_air': {'value': 0., 'units': 'dimensionless'},
    'cloud_area_fraction_in_atmosphere_layer': {'value': 0., 'units': 'dimensionless'},
    'shortwave_optical_thickness_due_to_aerosol': {'value': 0., 'units': 'dimensionless'},
    'longwave_optical_thickness_due_to_aerosol': {'value': 0., 'units': 'dimensionless'},
    'mass_content_of_cloud_ice_in_atmosphere_layer': {'value': 0., 'units': 'kg m^-2'},
    'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {'value': 0., 'units': 'kg m^-2'},
    'cloud_ice_particle_size': {'value': 20., 'units': 'micrometer'},
    'cloud_water_droplet_radius': {'value': 10., 'units': 'micrometer'},
    'longwave_optical_thickness_due_to_cloud': {'value': 0., 'units': 'dimensionless'},
    'shortwave_optical_thickness_due_to_cloud': {'value': 0., 'units': 'dimensionless'},
    'snow_and_ice_temperature_spline': {'value': CubicSpline(np.linspace(0, 50, 50), 260.*np.ones(50)), 'units': 'degK'},
    'cloud_base_mass_flux': {'value': 0, 'units': 'kg m^-2 s^-1'},
}


def aggregate_input_properties(component_list):
    """
    Takes in a list of objects with an input_properties dictionary as an
    attribute, and returns an input_properties dictionary that satisfies all
    of those components.
    """
    # This has already been implemented in Sympl, see how the composites do
    # it to make their input properties. May require exporting some currently
    # internal functions from Sympl.
    raise NotImplementedError()


def get_init_diagnostic(name, properties, grid_state):
    """
    Takes in a quantity name, properties dictionary, and grid state. Returns
    a Diagnostic object which calculates that quantity.
    """
    # First, check if the quantity is in the default_values dict, and return
    # a constructed ConstantDefaultValue Diagnostic if it is.
    on_interface = name[-20:] == '_on_interface_levels'
    if on_interface:
        base_name = name[:-20]
    else:
        base_name = name
    if base_name in default_values:
        p_name = 'air_pressure'
        if on_interface:
            p_name += '_on_interface_levels'
        return ConstantDefaultValue(
            {p_name: 'Pa'},
            name,
            default_values[base_name]['value'],
            default_values[base_name]['units']
        )
    # If it isn't, check if there is a diagnostic defined in some library of Diagnostic
    # classes (probably a list stored here) that can calculate the quantity,
    # and return that one.
    else:
        diagnostic = get_init_diagnostic(name)
        if diagnostic is not None:
            return diagnostic
        # If it still isn't, raise an exception because we can't calculate it.
        else:
            raise NotImplementedError('No initialization method for quantity name {}'.format(name))


def get_diagnostics_for(input_properties, grid_state):
    diagnostic_list = []
    for name, properties in input_properties.items():
        diagnostic_list.append(get_init_diagnostic(name, properties, grid_state))
    return diagnostic_list


def compute_all_diagnostics(state, diagnostic_list):
    return_dict = {}
    for diagnostic in diagnostic_list:
        return_dict.update(diagnostic(state))
    return return_dict


def get_default_state(component_list, grid_state=None):
    grid_state = grid_state or get_grid()
    input_properties = aggregate_input_properties(component_list)
    diagnostic_list = get_diagnostics_for(input_properties, grid_state)
    return_state = {}
    return_state.update(grid_state)
    return_state.update(compute_all_diagnostics(grid_state, diagnostic_list))
    return return_state

init_diagnostics = [
    PressureFunctionDiagnostic(
        'longwave_optical_depth',
        lambda p: 1. * (1. - p),
        'dimensionless',
        'interface',
    )
]

def get_init_diagnostic(name):
    for diag in init_diagnostics:
        if name in diag.diagnostic_properties:
            return diag
    return None

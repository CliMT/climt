from sympl import (
    DataArray, ensure_no_shared_keys, Diagnostic, combine_component_properties)
from .._components import RRTMGShortwave, RRTMGLongwave
import numpy as np
from datetime import datetime
from scipy.interpolate import CubicSpline
import pkg_resources
import numpy.linalg as la
from numpy.polynomial.legendre import legcompanion, legder, legval

class RRTMGLongwaveDefaultValues(Diagnostic):

    input_properties = {
        'air_pressure': {
            'dims': ['*', 'mid_levels'],
            'units': 'Pa',
        },
    }

    diagnostic_properties = {
        'surface_longwave_emissivity': {
            'dims': ['num_longwave_bands', '*'],
            'units': 'dimensionless',
        },
        'longwave_optical_thickness_due_to_cloud': {
            'dims': ['mid_levels', '*', 'num_longwave_bands'],
            'units': 'dimensionless',
        },
        'longwave_optical_thickness_due_to_aerosol': {
            'dims': ['mid_levels', '*', 'num_longwave_bands'],
            'units': 'dimensionless',
        },
    }

    def array_call(self, state):
        lw_band_shape = list(state['air_pressure'].shape) + [RRTMGLongwave.num_longwave_bands]
        diagnostics = {
            'surface_longwave_emissivity': np.ones(
                [state['air_pressure'].shape[0], RRTMGLongwave.num_longwave_bands]
            ),
            'longwave_optical_thickness_due_to_cloud': np.zeros(lw_band_shape),
            'longwave_optical_thickness_due_to_aerosol': np.zeros(lw_band_shape),
        }
        return diagnostics


class RRTMGShortwaveDefaultValues(Diagnostic):

    input_properties = {
        'air_pressure': {
            'dims': ['*', 'mid_levels'],
            'units': 'Pa',
        },
    }

    diagnostic_properties = {
        'shortwave_optical_thickness_due_to_cloud': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'cloud_asymmetry_parameter': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'cloud_forward_scattering_fraction': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'single_scattering_albedo_due_to_cloud': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'shortwave_optical_thickness_due_to_aerosol': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'aerosol_asymmetry_parameter': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'single_scattering_albedo_due_to_aerosol': {
            'dims': ['*', 'num_shortwave_bands', 'mid_levels'],
            'units': 'dimensionless',
        },
        'aerosol_optical_depth_at_55_micron': {
            'dims': ['*', 'num_ecmwf_aerosols', 'mid_levels'],
            'units': 'dimensionless',
        },
    }

    def array_call(self, state):
        p = state['air_pressure']
        sw_band_shape = [p.shape[0], RRTMGShortwave.num_shortwave_bands, p.shape[1]]
        ecmwf_aerosol_shape = [p.shape[0], RRTMGShortwave.num_ecmwf_aerosols, p.shape[1]]
        diagnostics = {
            'shortwave_optical_thickness_due_to_cloud': np.zeros(sw_band_shape),
            'cloud_asymmetry_parameter': 0.85 * np.ones(sw_band_shape),
            'cloud_forward_scattering_fraction': 0.8 * np.ones(sw_band_shape),
            'single_scattering_albedo_due_to_cloud': 0.9 * np.ones(sw_band_shape),
            'shortwave_optical_thickness_due_to_aerosol': np.zeros(sw_band_shape),
            'aerosol_asymmetry_parameter': np.zeros(sw_band_shape),
            'single_scattering_albedo_due_to_aerosol': 0.5 * np.ones(sw_band_shape),
            'aerosol_optical_depth_at_55_micron': np.zeros(ecmwf_aerosol_shape),
        }
        return diagnostics


class ConstantDefaultValue(Diagnostic):

    input_properties = {}
    diagnostic_properties = {}

    def __init__(
            self, output_name, output_value, output_units,
            dtype=None, **kwargs):
        if dtype is None:
            self._dtype = np.float64
        else:
            self._dtype = dtype
        self.diagnostic_properties = {
            output_name : {
                'dims': [],
                'units': output_units,
            },
        }
        self._output_name = output_name
        self._output_value = output_value
        super(ConstantDefaultValue, self).__init__(**kwargs)

    def array_call(self, state):
        return {self._output_name: np.array(self._output_value, dtype=self._dtype)}


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
            self._output_name: self._output_function(raw_state['p'])
        }


def horizontal_broadcast_if_needed(output, nx=None, ny=None):
    """
    Take a tuple of numpy arrays output. If nx and ny are not None, broadcast
    each of those numpy arrays along new dimensions of the correct length
    before returning them.
    """
    output_list = list(output)
    if ny is not None:
        for i in range(len(output_list)):
            output_list[i] = expand_new_last_dim(output_list[i], ny)
    if nx is not None:
        for i in range(len(output_list)):
            output_list[i] = expand_new_last_dim(output_list[i], nx)
    return tuple(output_list)


def expand_new_last_dim(var, length):
    indexer = tuple(slice(0, n) for n in var.shape) + (None,)
    return np.repeat(var[indexer], length, axis=-1)


def gaussian_latitudes(n):
    """Construct latitudes and latitude bounds for a Gaussian grid.

    Args:

    * n:
        The Gaussian grid number (half the number of latitudes in the
        grid.
    Returns:
        A 2-tuple where the first element is a length `n` array of
        latitudes (in degrees) and the second element is an `(n, 2)`
        array of bounds.
    """
    if abs(int(n)) != n:
        raise ValueError('n must be a non-negative integer')
    nlat = 2 * n
    # Create the coefficients of the Legendre polynomial and construct the
    # companion matrix:
    cs = np.array([0] * nlat + [1], dtype=np.int)
    cm = legcompanion(cs)
    # Compute the eigenvalues of the companion matrix (the roots of the
    # Legendre polynomial) taking advantage of the fact that the matrix is
    # symmetric:
    roots = la.eigvalsh(cm)
    roots.sort()
    # Improve the roots by one application of Newton's method, using the
    # solved root as the initial guess:
    fx = legval(roots, cs)
    fpx = legval(roots, legder(cs))
    roots -= fx / fpx
    # The roots should exhibit symmetry, but with a sign change, so make sure
    # this is the case:
    roots = (roots - roots[::-1]) / 2.
    # Compute the Gaussian weights for each interval:
    fm = legval(roots, cs[1:])
    fm /= np.abs(fm).max()
    fpx /= np.abs(fpx).max()
    weights = 1. / (fm * fpx)
    # Weights should be symmetric and sum to two (unit weighting over the
    # interval [-1, 1]):
    weights = (weights + weights[::-1]) / 2.
    weights *= 2. / weights.sum()
    # Calculate the bounds from the weights, still on the interval [-1, 1]:
    bounds1d = np.empty([nlat + 1])
    bounds1d[0] = -1
    bounds1d[1:-1] = -1 + weights[:-1].cumsum()
    bounds1d[-1] = 1
    # Convert the bounds to degrees of latitude on [-90, 90]:
    bounds1d = np.rad2deg(np.arcsin(bounds1d))
    bounds2d = np.empty([nlat, 2])
    bounds2d[:, 0] = bounds1d[:-1]
    bounds2d[:, 1] = bounds1d[1:]
    # Convert the roots from the interval [-1, 1] to latitude values on the
    # interval [-90, 90] degrees:
    latitudes = np.rad2deg(np.arcsin(roots))
    return latitudes, bounds2d


def get_grid(
        nx=None, ny=None, nz=50, p_surf_in_Pa=1.0132e5, x_name='longitude',
        y_name='latitude', latitude_grid='gaussian'):
    """
    Args:
        nx : int, optional
            Number of longitudinal points.
        ny : int, optional
            Number of latitudinal points.
        nz : int, optional
            Number of vertical mid levels.
        p_surf_in_Pa : float, optional
            Surface pressure in Pa.
        x_name : str, optional
            Name of latitudinal dimension
        y_name : str, optional
            Name of longitudinal dimension
        latitude_grid : 'gaussian' or 'regular'
            Type of spacing to use for the latitudinal grid.

    Returns:
        grid_state: dict
            A model state containing grid quantities.
    """
    p, p_interface, p_surface, sigma, sigma_interface = horizontal_broadcast_if_needed(
        get_pressure_and_sigma_levels(nz, p_surf_in_Pa), nx, ny,
    )
    return_state = {}
    horizontal_dims = []
    if nx is not None:
        return_state['longitude'] = DataArray(
            np.linspace(-180., 180., nx*2+1, endpoint=True)[1:-1:2],
            dims=[x_name],
            attrs={'units': 'degrees_east'},
        )
        horizontal_dims.insert(0, x_name)
    if ny is not None:
        if latitude_grid.lower() == 'regular':
            return_state['latitude'] = DataArray(
                np.linspace(-90., 90., ny*2+1, endpoint=True)[1:-1:2],
                dims=[y_name],
                attrs={'units': 'degrees_north'},
            )
        elif latitude_grid.lower() == 'gaussian':
            if ny % 2 != 0:
                raise ValueError('When latitude_grid is gaussian, ny must be even.')
            lat, lat_interface = gaussian_latitudes(ny / 2)
            return_state['latitude'] = DataArray(
                lat,
                dims=[y_name],
                attrs={'units': 'degrees_north'},
            )
        horizontal_dims.insert(0, y_name)
    return_state.update({
        'time': datetime(2000, 1, 1),
        'air_pressure': DataArray(
                p,
                dims=['mid_levels'] + horizontal_dims,
                attrs={'units': 'Pa'}),
        'air_pressure_on_interface_levels': DataArray(
                p_interface,
                dims=['interface_levels'] + horizontal_dims,
                attrs={'units': 'Pa'}),
        'surface_air_pressure': DataArray(
                p_surface,
                dims=horizontal_dims,
                attrs={'units': 'Pa'}),
        'sigma': DataArray(
                sigma,
                dims=['mid_levels'] + horizontal_dims,
                attrs={'units': 'dimensionless'}),
        'sigma_on_interface_levels': DataArray(
                sigma_interface,
                dims=['interface_levels'] + horizontal_dims,
                attrs={'units': 'dimensionless'}),
    })
    return return_state


def get_pressure_and_sigma_levels(nz, p_surface):
    sigma = np.linspace(0.998, 0.001, nz)
    sigma_interface = np.zeros([nz+1])
    sigma_interface[0] = 1.
    sigma_interface[-1] = 0.
    sigma_interface[1:-1] = 0.5 * (sigma[:-1] + sigma[1:])
    p = p_surface * sigma
    p_interface = p_surface * sigma_interface
    return p, p_interface, np.asarray(p_surface), sigma, sigma_interface


default_values = {
    'air_temperature': {'value': 290., 'units': 'degK'},
    'surface_temperature': {'value': 300., 'units': 'degK', 'grid': 'surface_air_pressure'},
    'sea_surface_temperature': {'value': 300., 'units': 'degK', 'grid': 'surface_air_pressure'},
    'soil_surface_temperature': {'value': 300., 'units': 'degK', 'grid': 'surface_air_pressure'},
    'northward_wind': {'value': 0., 'units': 'm/s'},
    'eastward_wind': {'value': 0., 'units': 'm/s'},
    'divergence_of_wind': {'value': 0., 'units': 's^-1'},
    'atmosphere_relative_vorticity': {'value': 0., 'units': 's^-1'},
    'surface_geopotential': {'value': 0., 'units': 'm^2 s^-2', 'grid': 'surface_air_pressure'},
    'surface_longwave_emissivity': {'value': 0., 'units': 'dimensionless', 'grid': 'surface_air_pressure'},
    'specific_humidity': {'value': 0., 'units': 'kg/kg'},
    'surface_specific_humidity': {'value': 0., 'units': 'kg/kg', 'grid': 'surface_air_pressure'},
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
    'snow_and_ice_temperature_spline': {'value': CubicSpline(np.linspace(0, 50, 50), 260.*np.ones(50)), 'units': 'degK', 'dtype': object},
    'cloud_base_mass_flux': {'value': 0., 'units': 'kg m^-2 s^-1'},
    'surface_thermal_capacity': {'value': 4.1813e3, 'units': 'J kg^-1 degK^-1'},
    'depth_of_slab_surface': {'value': 50., 'units': 'm'},
    'surface_material_density': {'value': 1000., 'units': 'kg m^-3'},
    'solar_cycle_fraction': {'value': 0., 'units': 'dimensionless'},
    'flux_adjustment_for_earth_sun_distance': {'value': 1.0, 'units': 'dimensionless'},
    'sea_water_density': {'value': 1.029e3, 'units': 'kg m^-3'},
    'surface_albedo_for_direct_shortwave': {'value': 0.06, 'units': 'dimensionless'},
    'surface_albedo_for_diffuse_shortwave': {'value': 0.06, 'units': 'dimensionless'},
    'surface_albedo_for_direct_near_infrared': {'value': 0.06, 'units': 'dimensionless'},
    'surface_albedo_for_diffuse_near_infrared': {'value': 0.06, 'units': 'dimensionless'},
    'soil_type': {'value': 'clay', 'units': 'dimensionless', 'dtype': 'a100'},
    'surface_roughness_length': {'value': 0.0002, 'units': 'dimensionless'},
    'surface_drag_coefficient_for_heat_in_air': {'value': 0.0012, 'units': 'dimensionless'},
    'surface_drag_coefficient_for_momentum_in_air': {'value': 0.0012, 'units': 'dimensionless'},
    'soil_temperature': {'value': 274., 'units': 'degK'},
    'zenith_angle': {'value': 0., 'units': 'radians'},
    'latitude': {'value': 0., 'units': 'degrees_north'},
    'longitude': {'value': 0., 'units': 'degrees_east'},
    'surface_downwelling_shortwave_flux': {'value': 0., 'units': 'W m^-2'},
    'surface_downwelling_longwave_flux': {'value': 0., 'units': 'W m^-2'},
    'surface_upwelling_shortwave_flux': {'value': 0., 'units': 'W m^-2'},
    'surface_upwelling_longwave_flux': {'value': 0., 'units': 'W m^-2'},
    'upward_heat_flux_at_ground_level_in_soil': {'value': 0., 'units': 'W m^-2'},
    'heat_flux_into_sea_water_due_to_sea_ice': {'value': 0., 'units': 'W m^-2'},
    'ocean_mixed_layer_thickness': {'value': 50., 'units': 'm'},
    'soil_layer_thickness': {'value': 50., 'units': 'm'},
    'heat_capacity_of_soil': {'value': 2000., 'units': 'J kg^-1 degK^-1'},
    'area_type': {'value': 'sea', 'units': 'dimensionless', 'dtype': 'a100'},
    'surface_upward_sensible_heat_flux': {'value': 0., 'units': 'W m^-2'},
    'surface_upward_latent_heat_flux': {'value': 0., 'units': 'W m^-2'},
    'land_ice_thickness': {'value': 0., 'units': 'm'},
    'sea_ice_thickness': {'value': 0., 'units': 'm'},
    'surface_snow_thickness': {'value': 0., 'units': 'm'},
}
for d in default_values.values():
    if 'grid' not in d.keys():
        d['grid'] = 'air_pressure'


def aggregate_input_properties(component_list):
    """
    Takes in a list of objects with an input_properties dictionary as an
    attribute, and returns an input_properties dictionary that satisfies all
    of those components.
    """
    return combine_component_properties(component_list, 'input_properties')


def get_init_diagnostic(name):
    """
    Takes in a quantity name. Returns a Diagnostic object which calculates that
    quantity from a grid state.
    """
    # First, check if the quantity is in the default_values dict, and return
    # a constructed ConstantDefaultValue Diagnostic if it is.
    if name in default_values:
        return ConstantDefaultValue(
            name,
            default_values[name]['value'],
            default_values[name]['units'],
            dtype=default_values[name].get('dtype', None)
        )
    # If it isn't, check if there is a diagnostic defined in some library of Diagnostic
    # classes (probably a list stored here) that can calculate the quantity,
    # and return that one.
    else:
        for diag in init_diagnostics:
            if name in diag.diagnostic_properties:
                return diag
    # If it still isn't, raise an exception because we can't calculate it.
    raise NotImplementedError('No initialization method for quantity name {}'.format(name))


def get_diagnostics_for(input_properties, grid_state):
    diagnostic_list = []
    for name in input_properties.keys():
        if name not in grid_state.keys():
            diagnostic_list.append(get_init_diagnostic(name))
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

def init_ozone(p):
    p_ref = np.linspace(0.998, 0.001, 30)
    ozone_ref = np.load(
        pkg_resources.resource_filename('climt._data', 'ozone_profile.npy')
    )
    spline = CubicSpline(p_ref[::-1], ozone_ref[::-1])  # x must be increasing
    return spline(p)

init_diagnostics = [
    PressureFunctionDiagnostic(
        'longwave_optical_depth',
        lambda p: 1. * (1. - p),
        'dimensionless',
        'interface',
    ),
    PressureFunctionDiagnostic(
        'mole_fraction_of_ozone_in_air',
        init_ozone,
        'mole/mole',
        'mid'
    ),
    RRTMGShortwaveDefaultValues(),
    RRTMGLongwaveDefaultValues(),
]

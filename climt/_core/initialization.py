from sympl import (
    DataArray, Diagnostic, combine_component_properties, get_constant,
)
from .._components import RRTMGShortwave, RRTMGLongwave
import numpy as np
from datetime import datetime
from scipy.interpolate import CubicSpline
import pkg_resources
import numpy.linalg as la
from numpy.polynomial.legendre import legcompanion, legder, legval

a_coord_spline = CubicSpline(
    np.linspace(0, 1, 29, endpoint=True),
    np.array([
        0.0, 0.00899999961, 11.6350002, 86.4570007, 292.576996, 701.453003,
        1381.526, 2389.20508, 3757.14209, 5481.14404, 7508.53223, 9731.80664,
        11991.4277, 14089.5352, 15812.9258, 16959.9941, 17364.6582, 16912.1309,
        15613.5645, 13671.4268, 11343.6543, 8913.7666, 6678.60791, 4844.03613,
        3376.51611, 2210.979, 1290.53296, 566.89801, 200.0
    ])
)

b_coord_spline = CubicSpline(
    np.linspace(0, 1, 29, endpoint=True),
    np.array([
        1.0, 0.988725841, 0.974401832, 0.955872416, 0.931749582, 0.900580883,
        0.860974848, 0.811784863, 0.752347112, 0.682746828, 0.604054928,
        0.518456697, 0.429195374, 0.340293199, 0.256084293, 0.180667043,
        0.117417611, 0.0686749965, 0.0349500999, 0.0143262697, 0.0039327601,
        0.000378680008, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ])
)

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
            'dims': ['num_longwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless',
        },
    }

    def array_call(self, state):
        ncol, nz = state['air_pressure'].shape
        diagnostics = {
            'surface_longwave_emissivity': np.ones(
                [RRTMGLongwave.num_longwave_bands, ncol]
            ),
            'longwave_optical_thickness_due_to_cloud': np.zeros(
                [nz, ncol, RRTMGLongwave.num_longwave_bands]),
            'longwave_optical_thickness_due_to_aerosol': np.zeros(
                [RRTMGLongwave.num_longwave_bands, nz, ncol]),
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
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless',
        },
        'cloud_asymmetry_parameter': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless',
        },
        'cloud_forward_scattering_fraction': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless',
        },
        'single_scattering_albedo_due_to_cloud': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless',
        },
        'shortwave_optical_thickness_due_to_aerosol': {
            'dims': ['num_shortwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless',
        },
        'aerosol_asymmetry_parameter': {
            'dims': ['num_shortwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless',
        },
        'single_scattering_albedo_due_to_aerosol': {
            'dims': ['num_shortwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless',
        },
        'aerosol_optical_depth_at_55_micron': {
            'dims': ['num_ecmwf_aerosols', 'mid_levels', '*'],
            'units': 'dimensionless',
        },
    }

    def array_call(self, state):
        ncol, nz = state['air_pressure'].shape
        diagnostics = {
            'shortwave_optical_thickness_due_to_cloud':
                np.zeros([nz, ncol, RRTMGShortwave.num_shortwave_bands]),
            'cloud_asymmetry_parameter':
                0.85 * np.ones([nz, ncol, RRTMGShortwave.num_shortwave_bands]),
            'cloud_forward_scattering_fraction':
                0.8 * np.ones([nz, ncol, RRTMGShortwave.num_shortwave_bands]),
            'single_scattering_albedo_due_to_cloud':
                0.9 * np.ones([nz, ncol, RRTMGShortwave.num_shortwave_bands]),
            'shortwave_optical_thickness_due_to_aerosol':
                np.zeros([nz, ncol, RRTMGShortwave.num_shortwave_bands]),
            'aerosol_asymmetry_parameter':
                np.zeros([RRTMGShortwave.num_shortwave_bands, nz, ncol]),
            'single_scattering_albedo_due_to_aerosol':
                0.5 * np.ones([RRTMGShortwave.num_shortwave_bands, nz, ncol]),
            'aerosol_optical_depth_at_55_micron':
                np.zeros([RRTMGShortwave.num_ecmwf_aerosols, nz, ncol]),
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


def leggauss(deg):
    """
    Gauss-Legendre quadrature.
    Computes the sample points and weights for Gauss-Legendre quadrature.
    These sample points and weights will correctly integrate polynomials of
    degree :math:`2*deg - 1` or less over the interval :math:`[-1, 1]` with
    the weight function :math:`f(x) = 1`.
    Parameters
    ----------
    deg : int
        Number of sample points and weights. It must be >= 1.
    Returns
    -------
    x : ndarray
        1-D ndarray containing the sample points.
    y : ndarray
        1-D ndarray containing the weights.
    Notes
    -----
    .. versionadded:: 1.7.0
    The results have only been tested up to degree 100, higher degrees may
    be problematic. The weights are determined by using the fact that
    .. math:: w_k = c / (L'_n(x_k) * L_{n-1}(x_k))
    where :math:`c` is a constant independent of :math:`k` and :math:`x_k`
    is the k'th root of :math:`L_n`, and then scaling the results to get
    the right value when integrating 1.
    """
    ideg = int(deg)
    if ideg != deg or ideg < 1:
        raise ValueError("deg must be a non-negative integer")

    # first approximation of roots. We use the fact that the companion
    # matrix is symmetric in this case in order to obtain better zeros.
    c = np.array([0]*deg + [1])
    m = np.polynomial.legendre.legcompanion(c)
    x = np.linalg.eigvalsh(m)

    # improve roots by one application of Newton
    dy = np.polynomial.legendre.legval(x, c)
    df = np.polynomial.legendre.legval(x, np.polynomial.legendre.legder(c))
    x -= dy/df

    # compute the weights. We scale the factor to avoid possible numerical
    # overflow.
    fm = np.polynomial.legendre.legval(x, c[1:])
    fm /= np.abs(fm).max()
    df /= np.abs(df).max()
    w = 1/(fm * df)

    # for Legendre we can also symmetrize
    w = (w + w[::-1])/2
    x = (x - x[::-1])/2

    # scale w to get the right value
    w *= 2. / w.sum()

    return x, w


def gaussian_latitudes(n):
    x, weights = leggauss(n)
    edges = np.empty([n+1])
    edges[0] = -1
    edges[1:-1] = -1 + weights[:-1].cumsum()
    edges[-1] = 1
    return -np.rad2deg(np.arcsin(x)), -np.rad2deg(np.arcsin(edges))


def get_grid(
        nx=None, ny=None, nz=28, p_surf_in_Pa=1.0132e5, x_name='longitude',
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
    return_state = get_hybrid_sigma_pressure_levels(nz)
    return_state['surface_air_pressure'] = DataArray(
        p_surf_in_Pa, dims=[], attrs={'units': 'Pa'}
    )
    return_state['time'] = datetime(2000, 1, 1)
    return_state.update(HybridSigmaPressureDiagnostic()(return_state))
    if nx is not None:
        return_state['longitude'] = DataArray(
            np.linspace(0., 360., nx*2, endpoint=False)[:-1:2],
            dims=[x_name],
            attrs={'units': 'degrees_east'},
        )
    if ny is not None:
        if latitude_grid.lower() == 'regular':
            return_state['latitude'] = DataArray(
                np.linspace(-90., 90., ny*2+1, endpoint=True)[1:-1:2],
                dims=[y_name],
                attrs={'units': 'degrees_north'},
            )
        elif latitude_grid.lower() == 'gaussian':
            lat, lat_interface = gaussian_latitudes(ny)
            return_state['latitude'] = DataArray(
                lat,
                dims=[y_name],
                attrs={'units': 'degrees_north'},
            )
    return return_state


class HybridSigmaPressureDiagnostic(Diagnostic):

    input_properties = {
        'atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels': {
            'units': 'dimensionless',
            'dims': ['interface_levels', '*'],
            'alias': 'a_coord',
        },
        'atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels': {
            'units': 'dimensionless',
            'dims': ['interface_levels', '*'],
            'alias': 'b_coord',
        },
        'surface_air_pressure': {
            'units': 'Pa',
            'dims': ['*'],
        },
    }

    diagnostic_properties = {
        'air_pressure': {
            'units': 'Pa',
            'dims': ['mid_levels', '*'],
        },
        'air_pressure_on_interface_levels': {
            'units': 'Pa',
            'dims': ['interface_levels', '*'],
        },
    }

    def array_call(self, state):
        p_interface = (
            state['a_coord'] +
            state['b_coord'] * state['surface_air_pressure'][None, :])
        delta_p = p_interface[1:, :] - p_interface[:-1, :]
        Rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1')
        rk = Rd/Cpd
        p = (
            (p_interface[1:, :]**(rk+1) - p_interface[:-1, :]**(rk+1)) / (
                (rk+1) * delta_p
            )
        ) ** (1./rk)
        assert not np.any(np.isnan(p))
        return {
            'air_pressure': p,
            'air_pressure_on_interface_levels': p_interface,
        }


def get_hybrid_sigma_pressure_levels(nz):
    a_interface = a_coord_spline(np.linspace(0., 1., nz+1, endpoint=True))
    b_interface = b_coord_spline(np.linspace(0., 1., nz+1, endpoint=True))
    return {
        'atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels':
            DataArray(
                a_interface,
                dims=['interface_levels'],
                attrs={'units': 'dimensionless'},
            ),
        'atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels':
            DataArray(
                b_interface,
                dims=['interface_levels'],
                attrs={'units': 'dimensionless'},
            ),
    }


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
    'mass_content_of_cloud_ice_in_atmosphere_layer': {'value': 0., 'units': 'kg m^-2'},
    'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {'value': 0., 'units': 'kg m^-2'},
    'cloud_ice_particle_size': {'value': 20., 'units': 'micrometer'},
    'cloud_water_droplet_radius': {'value': 10., 'units': 'micrometer'},
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

for diag in init_diagnostics:
    output_overlap = set(diag.diagnostic_properties.keys()).intersection(
        default_values.keys())
    if len(output_overlap) > 0:
        raise AssertionError(
            'Initialization diagnostic {} outputs quantity(ies) {} for which a '
            'default value is also present. The default value should probably '
            'be removed.'.format(diag.__class__.__name__, output_overlap))
    for diag2 in init_diagnostics:
        output_overlap = set(diag.diagnostic_properties.keys()).intersection(
            diag2.diagnostic_properties.keys())
        if diag2 is diag:
            pass
        elif len(output_overlap) > 0:
            raise AssertionError(
                'Initialization diagnostics {} and {} both outputs quantity(ies)'
                ' {}. One of these should probably '
                'be removed.'.format(
                    diag.__class__.__name__,
                    diag2.__class__.__name__,
                    output_overlap
                )
            )

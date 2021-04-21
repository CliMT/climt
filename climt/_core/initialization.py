from sympl import (
    DataArray, DiagnosticComponent, combine_component_properties, get_constant,
    set_constant
)
from .._components import RRTMGShortwave, RRTMGLongwave
import numpy as np
from datetime import datetime
from scipy.interpolate import CubicSpline
import pkg_resources


def get_atmosphere_grid(grid_state,
                        interface=False,
                        horizontal=False):
    y, x = grid_state['latitude'].shape
    y_dim, x_dim = grid_state['latitude'].dims
    z = grid_state['atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels'].shape[0]

    if horizontal:
        return (y, x), (y_dim, x_dim)
    if interface:
        return (z, y, x), ('interface_levels', y_dim, x_dim)
    else:
        return (z-1, y, x), ('mid_levels', y_dim, x_dim)


def get_surface_grid(grid_state,
                     interface=False,
                     horizontal=False):
    return tuple(grid_state['latitude'].shape), tuple(grid_state['latitude'].dims)


def get_land_grid(grid_state,
                  interface=False,
                  horizontal=False):
    y, x = grid_state['latitude'].shape
    y_dim, x_dim = grid_state['latitude'].dims

    if horizontal:
        return (y, x), (y_dim, x_dim)
    else:
        raise NotImplementedError(
            '3D land grids are not yet supported')


def get_ocean_grid(grid_state,
                   interface=False,
                   horizontal=False):
    y, x = grid_state['latitude'].shape
    y_dim, x_dim = grid_state['latitude'].dims

    if horizontal:
        return (y, x), (y_dim, x_dim)
    else:
        raise NotImplementedError(
            '3D ocean grids are not yet supported')


def get_ice_grid(grid_state,
                 interface=False,
                 horizontal=False):
    y, x = grid_state['latitude'].shape
    y_dim, x_dim = grid_state['latitude'].dims
    z = grid_state['height_on_ice_interface_levels'].shape[0]

    if horizontal:
        return (y, x), (y_dim, x_dim)
    if interface:
        return (z, y, x), ('ice_interface_levels', y_dim, x_dim)
    else:
        return (z-1, y, x), ('ice_mid_levels', y_dim, x_dim)


def get_scalar_grid(grid_state,
                    interface=False,
                    horizontal=False):
    return (), ()


domain_shape_descriptor = {
    'atmosphere': get_atmosphere_grid,
    'surface': get_surface_grid,
    'land': get_land_grid,
    'ocean': get_ocean_grid,
    'ice': get_ice_grid,
    'scalar': get_scalar_grid,
}


class RRTMGLongwaveDefaultValues(DiagnosticComponent):

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


class RRTMGShortwaveDefaultValues(DiagnosticComponent):

    input_properties = {
        'air_pressure': {
            'dims': ['mid_levels', '*'],
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
        nz, ncol = state['air_pressure'].shape
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
                np.zeros([RRTMGShortwave.num_shortwave_bands, nz, ncol]),
            'aerosol_asymmetry_parameter':
                np.zeros([RRTMGShortwave.num_shortwave_bands, nz, ncol]),
            'single_scattering_albedo_due_to_aerosol':
                0.5 * np.ones([RRTMGShortwave.num_shortwave_bands, nz, ncol]),
            'aerosol_optical_depth_at_55_micron':
                np.zeros([RRTMGShortwave.num_ecmwf_aerosols, nz, ncol]),
        }
        return diagnostics


class ConstantDefaultValue(DiagnosticComponent):

    input_properties = {}
    diagnostic_properties = {}

    def __init__(
            self, output_name, output_value, output_units,
            dtype=None, domain=None, **kwargs):
        if dtype is None:
            self._dtype = np.float64
        else:
            self._dtype = dtype
        self._output_name = output_name
        self._output_value = output_value
        self._output_units = output_units
        self.store_domain_properties(domain)

        self.diagnostic_properties = {
            output_name: {
                'dims': ['*'],
                'units': output_units,
            },
        }
        super(ConstantDefaultValue, self).__init__(**kwargs)

    def store_domain_properties(self, domain):
        self._interface = False
        self._horizontal = False
        if domain is None:
            self._domain = 'scalar'
            return

        domain_and_type = domain.split('_')
        self._domain = domain_and_type[0]

        if len(domain_and_type) > 1:
            if domain_and_type[1] == 'horizontal':
                self._horizontal = True
            elif domain_and_type[1] == 'interface':
                self._interface = True
            else:
                NotImplementedError(
                    'Unknown domain descriptor {}'.format(domain))

    def __call__(self, grid_state):
        shape, dims = domain_shape_descriptor[self._domain](
            grid_state, self._interface, self._horizontal)

        quantity_np_array = np.broadcast_to(
            np.array(self._output_value, dtype=self._dtype),
            shape).copy()

        quantity_array = DataArray(
            quantity_np_array,
            dims=dims,
            name=self._output_name,
            attrs=dict(units=self._output_units))

        return {self._output_name: quantity_array}

    def array_call(self, state):
        return


class PressureFunctionDiagnosticComponent(DiagnosticComponent):
    """Defines a quantity as a function of pressure."""

    input_properties = {
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa',
            'alias': 'p'
        },
        'surface_air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
            'alias': 'ps'
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
        vertical_dimension = 'mid_levels'
        if mid_or_interface_levels == 'interface':
            vertical_dimension = 'interface_levels'
            output_name += '_on_interface_levels'
            self.input_properties = {
                'air_pressure_on_interface_levels': {
                    'dims': ['interface_levels', '*'],
                    'units': 'Pa',
                    'alias': 'p'
                },
                'surface_air_pressure': {
                    'dims': ['*'],
                    'units': 'Pa',
                    'alias': 'ps'
                }
            }
        elif mid_or_interface_levels == 'mid':
            pass
        else:
            raise ValueError(
                "Argument mid_or_interface_levels must be 'mid' or 'interface'")
        self.diagnostic_properties = {
            output_name: {
                'dims': [vertical_dimension, '*'],
                'units': output_units
            }
        }
        self._output_function = output_function
        self._output_name = output_name
        super(PressureFunctionDiagnosticComponent, self).__init__()

    def array_call(self, raw_state):
        return {
            self._output_name: self._output_function(raw_state['p'], raw_state['ps'])
        }


'''
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
'''


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
        nx=None, ny=None, nz=28, n_ice_interface_levels=10,
        p_surf_in_Pa=None, p_toa_in_Pa=None,
        proportion_sigma_levels=0.1,
        proportion_isobaric_levels=0.25,
        x_name='lon', y_name='lat',
        latitude_grid='gaussian'):
    """
    Args:
        nx : int, optional
            Number of longitudinal points.
        ny : int, optional
            Number of latitudinal points.
        nz : int, optional
            Number of vertical mid levels.
        n_ice_interface_levels (int, optional): Number of vertical
            interface levels to use for ice. Use None to disable the ice
            vertical grid.
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
    if p_surf_in_Pa is None:
        p_surf_in_Pa = get_constant('reference_air_pressure', 'Pa')

    if p_toa_in_Pa is None:
        p_toa_in_Pa = get_constant('top_of_model_pressure', 'Pa')
    else:
        set_constant('top_of_model_pressure', p_toa_in_Pa, 'Pa')

    if nx is None:
        nx = 1
    if ny is None:
        ny = 1

    return_state = get_hybrid_sigma_pressure_levels(nz+1,
                                                    p_surf_in_Pa,
                                                    p_toa_in_Pa,
                                                    proportion_isobaric_levels,
                                                    proportion_sigma_levels)
    return_state['surface_air_pressure'] = DataArray(
        np.ones((ny, nx))*p_surf_in_Pa, dims=[y_name, x_name], attrs={'units': 'Pa'}
    )
    return_state['time'] = datetime(2000, 1, 1)
    return_state.update(HybridSigmaPressureDiagnosticComponent()(return_state))
    if nx is not None:
        two_dim_lons = np.ones((ny, nx))
        two_dim_lons[:] = np.linspace(0., 360., nx*2, endpoint=False)[:-1:2][np.newaxis, :]
        return_state['longitude'] = DataArray(
            two_dim_lons,
            dims=[y_name, x_name],
            attrs={'units': 'degrees_east'},
        )
    if ny is not None:
        two_dim_lats = np.ones((ny, nx))
        if latitude_grid.lower() == 'regular':
            two_dim_lats[:] = np.linspace(-90., 90., ny*2+1, endpoint=True)[1:-1:2][:, np.newaxis]
            return_state['latitude'] = DataArray(
                two_dim_lats,
                dims=[y_name, x_name],
                attrs={'units': 'degrees_north'},
            )
        elif latitude_grid.lower() == 'gaussian':
            lat, lat_interface = gaussian_latitudes(ny)
            two_dim_lats[:] = lat[:, np.newaxis]
            return_state['latitude'] = DataArray(
                two_dim_lats,
                dims=[y_name, x_name],
                attrs={'units': 'degrees_north'},
            )
        else:
            raise ValueError(
                'latitude_grid can be either regular or gaussian. ' +
                'Other grid types are currently not supported.')
    if n_ice_interface_levels is not None:
        return_state['height_on_ice_interface_levels'] = DataArray(
            np.zeros(n_ice_interface_levels),
            dims=['ice_interface_levels'],
            attrs={'units': 'm'},
        )
    return return_state


class HybridSigmaPressureDiagnosticComponent(DiagnosticComponent):

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
        model_top_pressure = get_constant('top_of_model_pressure', 'Pa')
        p_interface = (
            state['a_coord'] +
            state['b_coord'] * (state['surface_air_pressure'][None, :] - model_top_pressure))
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


def get_hybrid_sigma_pressure_levels(num_levels=28,
                                     reference_pressure=1e5,
                                     model_top_pressure=20,
                                     proportion_isobaric_levels=0.25,
                                     proportion_sigma_levels=0.1):

    '''
    Calculate the values of ak and bk using the algorithm from
    `[Eckermann]`_ (2008) for their NEWHYB2 coordinate system.

    The pressure thickness distribution is given by a sine shaped
    curve with a maximum at the middle of the pressure range. The
    nominal lower surface pressure of 1000 mb is used to calculate
    the coordinates, but this does not put any restriction on the
    actual surface pressure.

    Args:
        num_levels: int, optional
            The number of vertical levels
        reference_pressure: float, optional
            The reference surface pressure
        model_top_pressure: float, optional
            The pressure at the top of the model
        proportion_isobaric_levels: float, optional
            The proportion of levels where bk = 0
        proportion_sigma_levels: float, optional
            The proportion of levels where ak = 0

    Returns:
        hybrid_coords: dict
            dictionary containing ak, bk, and sigma levels

    Reference:
        Eckermann, S: Hybrid \sigma-p Coordinate Choices for a Global Model,
        Monthly Weather Review Jan 2009.

    .. _[Eckermann]:
        https://journals.ametsoc.org/doi/pdf/10.1175/2008MWR2537.1

    '''

    thickness_dist = np.sin(np.linspace(0.1, np.pi-0.1, num_levels-1))

    thickness_dist /= np.sum(thickness_dist)
    thickness_dist *= (reference_pressure - model_top_pressure)

    pressure_levels = np.zeros(num_levels)
    pressure_levels[0] = model_top_pressure
    pressure_levels[1:] = model_top_pressure + np.cumsum(thickness_dist)

    sigma_interface = (pressure_levels - model_top_pressure)/(reference_pressure -
                                                              model_top_pressure)

    ak = np.zeros(num_levels)
    bk = np.zeros(num_levels)

    num_isobaric_levels = int(proportion_isobaric_levels*num_levels)
    num_sigma_levels = int(proportion_sigma_levels*num_levels)

    ak[0:num_isobaric_levels] = pressure_levels[0:num_isobaric_levels]

    isobaric_sigma_level = sigma_interface[num_isobaric_levels-1]

    for level in range(num_isobaric_levels, num_levels - num_sigma_levels):

        sigma_value = sigma_interface[level]
        b_level = (sigma_value - isobaric_sigma_level)/(1 - isobaric_sigma_level)

        r_level = get_exponent_for_sigma(b_level, num_sigma_levels)

        bk[level] = b_level**r_level

        ak[level] = model_top_pressure + (sigma_value - bk[level])*(reference_pressure -
                                                                    model_top_pressure)

    for level in range(num_levels - num_sigma_levels, num_levels):

        sigma_value = sigma_interface[level]
        bk[level] = (sigma_interface[level] - isobaric_sigma_level)/(1 - isobaric_sigma_level)
        ak[level] = model_top_pressure + (sigma_value - bk[level])*(reference_pressure -
                                                                    model_top_pressure)

    ak = ak[::-1]
    bk = bk[::-1]

    return {
        'atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels':
            DataArray(
                ak,
                dims=['interface_levels'],
                attrs={'units': 'dimensionless'},
            ),
        'atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels':
            DataArray(
                bk,
                dims=['interface_levels'],
                attrs={'units': 'dimensionless'},
            ),
    }


def get_exponent_for_sigma(b_half, num_sigma_levels):

    # These values correspond to the NEWHYB2 coordinate
    r_p = 2.2
    r_sigma = 1.0
    S = 5

    if num_sigma_levels > 0:
        r_sigma = 1
    else:
        r_sigma = 1.35

    return r_p + (r_sigma - r_p)*np.arctan(S*b_half)/np.arctan(S)


default_values = {
    'air_temperature': {'value': 290., 'units': 'degK', 'domain': 'atmosphere'},
    'northward_wind': {'value': 0., 'units': 'm/s', 'domain': 'atmosphere'},
    'eastward_wind': {'value': 0., 'units': 'm/s', 'domain': 'atmosphere'},
    'divergence_of_wind': {'value': 0., 'units': 's^-1', 'domain': 'atmosphere'},
    'atmosphere_relative_vorticity': {'value': 0., 'units': 's^-1', 'domain': 'atmosphere'},
    'specific_humidity': {'value': 0., 'units': 'kg/kg', 'domain': 'atmosphere'},
    'mole_fraction_of_carbon_dioxide_in_air': {'value': 330e-6, 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_methane_in_air': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_nitrous_oxide_in_air': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_oxygen_in_air': {'value': 0.21, 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_nitrogen_in_air': {'value': 0.78, 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_hydrogen_in_air': {'value': 500e-9, 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_cfc11_in_air': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_cfc12_in_air': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_cfc22_in_air': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mole_fraction_of_carbon_tetrachloride_in_air': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'cloud_area_fraction_in_atmosphere_layer': {'value': 0., 'units': 'dimensionless', 'domain': 'atmosphere'},
    'mass_content_of_cloud_ice_in_atmosphere_layer': {'value': 0., 'units': 'kg m^-2', 'domain': 'atmosphere'},
    'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {'value': 0., 'units': 'kg m^-2', 'domain': 'atmosphere'},
    'cloud_ice_particle_size': {'value': 20., 'units': 'micrometer', 'domain': 'atmosphere'},
    'cloud_water_droplet_radius': {'value': 10., 'units': 'micrometer', 'domain': 'atmosphere'},
    'cloud_base_mass_flux': {'value': 0., 'units': 'kg m^-2 s^-1', 'domain': 'atmosphere_horizontal'},
    'zenith_angle': {'value': 0., 'units': 'radians', 'domain': 'atmosphere_horizontal'},
    'downwelling_shortwave_flux_in_air': {'value': 0., 'units': 'W m^-2', 'domain':
                                          'atmosphere_interface'},
    'downwelling_longwave_flux_in_air': {'value': 0., 'units': 'W m^-2', 'domain':
                                         'atmosphere_interface'},
    'upwelling_shortwave_flux_in_air': {'value': 0., 'units': 'W m^-2', 'domain':
                                        'atmosphere_interface'},
    'upwelling_longwave_flux_in_air': {'value': 0., 'units': 'W m^-2', 'domain':
                                       'atmosphere_interface'},


    'surface_specific_humidity': {'value': 0., 'units': 'kg/kg', 'domain': 'surface'},
    'surface_temperature': {'value': 300., 'units': 'degK', 'domain': 'surface'},
    'soil_surface_temperature': {'value': 300., 'units': 'degK', 'domain': 'surface'},
    'surface_geopotential': {'value': 0., 'units': 'm^2 s^-2', 'domain': 'surface'},
    'surface_thermal_capacity': {'value': 4.1813e3, 'units': 'J kg^-1 degK^-1', 'domain': 'surface'},
    'depth_of_slab_surface': {'value': 50., 'units': 'm', 'domain': 'surface'},
    'surface_material_density': {'value': 1000., 'units': 'kg m^-3', 'domain': 'surface'},
    'surface_albedo_for_direct_shortwave': {'value': 0.06, 'units': 'dimensionless', 'domain': 'surface'},
    'surface_albedo_for_diffuse_shortwave': {'value': 0.06, 'units': 'dimensionless', 'domain': 'surface'},
    'surface_albedo_for_direct_near_infrared': {'value': 0.06, 'units': 'dimensionless', 'domain': 'surface'},
    'surface_albedo_for_diffuse_near_infrared': {'value': 0.06, 'units': 'dimensionless', 'domain': 'surface'},
    'surface_roughness_length': {'value': 0.0002, 'units': 'dimensionless', 'domain': 'surface'},
    'surface_drag_coefficient_for_heat_in_air': {'value': 0.0012, 'units': 'dimensionless', 'domain': 'surface'},
    'surface_drag_coefficient_for_momentum_in_air': {'value': 0.0012, 'units': 'dimensionless', 'domain': 'surface'},
    'area_type': {'value': 'sea', 'units': 'dimensionless', 'dtype': 'a100', 'domain': 'surface'},
    'surface_upward_sensible_heat_flux': {'value': 0., 'units': 'W m^-2', 'domain': 'surface'},
    'surface_upward_latent_heat_flux': {'value': 0., 'units': 'W m^-2', 'domain': 'surface'},

    'soil_type': {'value': 'clay', 'units': 'dimensionless', 'dtype': 'a100', 'domain': 'land_horizontal'},
    'soil_temperature': {'value': 274., 'units': 'degK', 'domain': 'land_horizontal'},
    'soil_layer_thickness': {'value': 50., 'units': 'm', 'domain': 'land_horizontal'},
    'upward_heat_flux_at_ground_level_in_soil': {'value': 0., 'units': 'W m^-2', 'domain': 'land_horizontal'},
    'heat_capacity_of_soil': {'value': 2000., 'units': 'J kg^-1 degK^-1', 'domain': 'land_horizontal'},

    'sea_water_density': {'value': 1.029e3, 'units': 'kg m^-3', 'domain': 'ocean_horizontal'},
    'sea_surface_temperature': {'value': 300., 'units': 'degK', 'domain': 'ocean_horizontal'},
    'ocean_mixed_layer_thickness': {'value': 50., 'units': 'm', 'domain': 'ocean_horizontal'},

    'snow_and_ice_temperature': {'value': 270., 'units': 'degK', 'domain': 'ice_interface'},
    'heat_flux_into_sea_water_due_to_sea_ice': {'value': 0., 'units': 'W m^-2', 'domain': 'ice_horizontal'},
    'land_ice_thickness': {'value': 0., 'units': 'm', 'domain': 'ice_horizontal'},
    'sea_ice_thickness': {'value': 0., 'units': 'm', 'domain': 'ice_horizontal'},
    'surface_snow_thickness': {'value': 0., 'units': 'm', 'domain': 'ice_horizontal'},

    'solar_cycle_fraction': {'value': 0., 'units': 'dimensionless', 'domain': None},
    'flux_adjustment_for_earth_sun_distance': {'value': 1.0, 'units': 'dimensionless', 'domain': None},
    'lwe_thickness_of_soil_moisture_content': {'value': 0, 'units': 'm', 'domain': 'surface'},
    'convective_precipitation_rate': {'value': 0., 'units': 'mm day^-1', 'domain': 'surface'},
    'stratiform_precipitation_rate': {'value': 0., 'units': 'm s^-1', 'domain': 'surface'},
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


def get_init_diagnostic(name, grid_state):
    """
    Takes in a quantity name. Returns a DiagnosticComponent object which calculates that
    quantity from a grid state.
    """
    # First, check if the quantity is in the default_values dict, and return
    # a constructed ConstantDefaultValue DiagnosticComponent if it is.
    if name in default_values:
        return ConstantDefaultValue(
            name,
            default_values[name]['value'],
            default_values[name]['units'],
            dtype=default_values[name].get('dtype', None),
            domain=default_values[name]['domain'],
        )
    elif name[-20:] == '_on_interface_levels' and name[:-20] in default_values:
        return ConstantDefaultValue(
            name,
            default_values[name[:-20]]['value'],
            default_values[name[:-20]]['units'],
            dtype=default_values[name[:-20]].get('dtype', None),
            domain=default_values[name[:-20]]['domain'] + '_interface',
        )
    # If it isn't, check if there is a diagnostic defined in some library of DiagnosticComponent
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
            diagnostic_list.append(get_init_diagnostic(name, grid_state))
    return diagnostic_list


def compute_all_diagnostics(state, diagnostic_list):
    return_dict = {}
    for diagnostic in diagnostic_list:
        return_dict.update(diagnostic(state))
    return return_dict


def get_default_state(
        component_list, grid_state=None, n_ice_interface_levels=30):
    """
    Retrieves a reasonable initial state for the set of components given.

    Args:
        component_list (list): Components for which to retrieve an
            initial state.
        grid_state (dict, optional): An initial state containing grid
            quantities. If none is given, a default will be created.
        n_ice_interface_levels (int, optional): Number of vertical
            interface levels to use for ice. Use None to disable the ice
            vertical grid.

    Returns:
        default_state (dict): A reasonable initial state.
    """
    grid_state = grid_state or get_grid(
        n_ice_interface_levels=n_ice_interface_levels)
    input_properties = aggregate_input_properties(component_list)
    diagnostic_list = get_diagnostics_for(input_properties, grid_state)
    return_state = {}
    return_state.update(grid_state)
    return_state.update(compute_all_diagnostics(grid_state, diagnostic_list))
    # return_state = broadcast_dims_to_match_component_properties(return_state, component_list)
    return return_state


def init_ozone(p, ps):
    p_ref = 1e5*np.linspace(0.998, 0.001, 30)
    ozone_ref = np.load(
        pkg_resources.resource_filename('climt._data', 'ozone_profile.npy')
    )
    spline = CubicSpline(p_ref[::-1], ozone_ref[::-1])  # x must be increasing
    return spline(p)


init_diagnostics = [
    PressureFunctionDiagnosticComponent(
        'longwave_optical_depth',
        lambda p, ps: 1. * (1. - p/ps),
        'dimensionless',
        'interface',
    ),
    PressureFunctionDiagnosticComponent(
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

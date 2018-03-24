from sympl import Diagnostic, DataArray
from .._core import get_constant
import xarray as xr
try:
    from ._berger_solar_insolation import get_solar_parameters, get_orbital_parameters
except ImportError:
    print('Import failed. Insolation will not be available!')


class BergerSolarInsolation(Diagnostic):
    """Determines solar insolation using spectral solutions for orbital
    constants from Berger 1978. This is the same approach taken by CAM 3.
    """
    inputs = ('time', 'longitude', 'latitude')
    diagnostics = (
        'solar_insolation',
        'solar_zenith_angle',
        'obliquity',
        'eccentricity',
        'normalized_earth_sun_distance',
    )

    def __init__(self):
        """
        Initialise insolation object
        """
        self._orbital_parameters = {}
        self._solar_constant = get_constant('stellar_irradiance', 'W/m^2')

    def __call__(self, state):
        """
        Gets diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            diagnostics (dict): A dictionary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Prognostic instance.
        """
        lat, lon = xr.broadcast(
            state['latitude'].to_units('degrees_N'),
            state['longitude'].to_units('degrees_E'))
        if 'solar_constant' in state:
            if state['solar_constant'].size > 1:
                raise ValueError(
                    'Solar constant should be a scalar (inside a DataArray)')
            else:
                solar_constant = state['solar_constant'].values.item()
        else:
            solar_constant = self._solar_constant
        solar_insolation, solar_zenith_angle, obliquity, eccentricity, rho = \
            self._driver(state, lat, lon, solar_constant)
        diagnostics = {
            'solar_insolation': DataArray(
                solar_insolation.reshape(lat.shape),
                dims=lat.dims, attrs={'units': 'W m^-2'}
            ),
            'solar_zenith_angle': DataArray(
                solar_zenith_angle.reshape(lat.shape),
                dims=lat.dims, attrs={'units': 'radians'}
            ),
            'obliquity': DataArray(
                obliquity, dims=[], attrs={'units': ''}
            ),
            'eccentricity': DataArray(
                eccentricity, dims=[], attrs={'units': ''}
            ),
            'normalized_earth_sun_distance': DataArray(
                rho, dims=[], attrs={'units': ''}
            ),
        }
        return diagnostics

    def _driver(self, state, lat, lon, solar_constant):
        # We can compute orbital parameters yearly because they change very
        # slowly
        year = state['time'].year
        if year not in self._orbital_parameters:
            self._orbital_parameters[year] = get_orbital_parameters(
                float(year - 1950))
        lambda_m0, eccentricity, omega_tilde, obliquity = (
            self._orbital_parameters[year])
        return get_solar_parameters(
            lambda_m0,
            eccentricity,
            omega_tilde,
            obliquity,
            years_since_vernal_equinox(state['time']),
            fractional_day(state['time']),
            lat.values.flatten(),
            lon.values.flatten(),
            solar_constant,)


def years_since_vernal_equinox(dt):
    """Fractional years since last March 20, noon UTC (assumed time of
    vernal equinox)."""
    year_start = type(dt)(dt.year, 3, 20, 12)
    year_end = type(dt)(dt.year+1, 3, 20, 12)
    return (dt - year_start).total_seconds() / (year_end - year_start).total_seconds()


def fractional_day(dt):
    day_start = type(dt)(dt.year, dt.month, dt.day)
    return (dt - day_start).total_seconds() / (24.*60.*60.)

from sympl import DiagnosticComponent, get_constant
import logging
try:
    from ._berger_solar_insolation import get_solar_parameters, get_orbital_parameters
except ImportError:
    logging.warning(
        'Import failed. Insolation is likely not compiled and will not be '
        'available.'
    )


class BergerSolarInsolation(DiagnosticComponent):
    """Determines solar insolation using spectral solutions for orbital
    constants from Berger 1978. This is the same approach taken by CAM 3.
    """
    input_properties = {
        'longitude': {
            'dims': ['*'],
            'units': 'degrees_east',
        },
        'latitude': {
            'dims': ['*'],
            'units': 'degrees_north',
        },
    }

    diagnostic_properties = {
        'solar_insolation': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'solar_zenith_angle': {
            'dims': ['*'],
            'units': 'radians',
        },
        'obliquity': {
            'dims': [],
            'units': 'radians',
        },
        'eccentricity': {
            'dims': [],
            'units': 'radians',
        },
        'normalized_earth_sun_distance': {
            'dims': [],
            'units': 'dimensionless',
        },
    }

    def __init__(self, **kwargs):
        self._orbital_parameters = {}
        super(BergerSolarInsolation, self).__init__(**kwargs)

    def array_call(self, state):
        solar_constant = get_constant('stellar_irradiance', 'W/m^2')
        solar_insolation, solar_zenith_angle, obliquity, eccentricity, rho = \
            self._driver(
                state['time'], state['latitude'], state['longitude'], solar_constant)
        return {
            'solar_insolation': solar_insolation,
            'solar_zenith_angle': solar_zenith_angle,
            'obliquity': obliquity,
            'eccentricity': eccentricity,
            'normalized_earth_sun_distance': rho,
        }

    def _driver(self, time, lat, lon, solar_constant):
        # We can compute orbital parameters yearly because they change very
        # slowly
        year = time.year
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
            years_since_vernal_equinox(time),
            fractional_day(time),
            lat,
            lon,
            solar_constant,
        )


def years_since_vernal_equinox(dt):
    """Fractional years since last March 20, noon UTC (assumed time of
    vernal equinox)."""
    year_start = type(dt)(dt.year, 3, 20, 12)
    year_end = type(dt)(dt.year+1, 3, 20, 12)
    return (dt - year_start).total_seconds() / (year_end - year_start).total_seconds()


def fractional_day(dt):
    day_start = type(dt)(dt.year, dt.month, dt.day)
    return (dt - day_start).total_seconds() / (24.*60.*60.)

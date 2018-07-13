from sympl import DiagnosticComponent
import datetime
import numpy as np


class Instellation(DiagnosticComponent):
    """
    Calculates the zenith angle and star-planet correction
    factor given orbital parameters. Currently useful only
    for Earth-sun system.
    """

    input_properties = {
        'latitude': {
            'dims': ['*'],
            'units': 'degrees_north',
        },
        'longitude': {
            'dims': ['*'],
            'units': 'degrees_east',
        },
    }

    diagnostic_properties = {
        'zenith_angle': {
            'dims': ['*'],
            'units': 'radians',
        }
    }

    def __init__(self, **kwargs):
        super(Instellation, self).__init__(**kwargs)

    def array_call(self, state):
        """
        Calculate zenith angle.

        Args:

            state (dict):
                state dictionary

        """
        lat_radians = np.deg2rad(state['latitude'])
        lon_radians = np.deg2rad(state['longitude'])
        zen_angle = sun_zenith_angle(state['time'], lon=lon_radians, lat=lat_radians)
        zen_angle[zen_angle > np.pi/2] = np.pi/2
        zen_angle[zen_angle < -np.pi/2] = -np.pi/2
        return {'zenith_angle': zen_angle}


def days_from_2000(model_time):
    """Get the days since year 2000.
    """
    return total_days(model_time - datetime.datetime(2000, 1, 1, 12, 0))


def total_days(time_diff):
    """
    Total time in units of days
    """
    return (time_diff.days +
            (time_diff.seconds +
             time_diff.microseconds / (1000000.0)) / (24 * 3600.0))


def greenwich_mean_sidereal_time(model_time):
    """
    Greenwich mean sidereal time, in radians.

    Reference:
        The AIAA 2006 implementation:
            http://www.celestrak.com/publications/AIAA/2006-6753/
    """
    jul_centuries = days_from_2000(model_time) / 36525.0
    theta = 67310.54841 + jul_centuries * (876600 * 3600 + 8640184.812866 + jul_centuries *
                                           (0.093104 - jul_centuries * 6.2 * 10e-6))

    theta_radians = np.deg2rad(theta / 240.0) % (2 * np.pi)

    if theta_radians < 0:
        theta_radians += 2*np.pi

    return theta_radians


def local_mean_sidereal_time(model_time, longitude):
    """
    Local mean sidereal time. requires longitude in radians.
    Ref:
        http://www.setileague.org/askdr/lmst.htm
    """
    return greenwich_mean_sidereal_time(model_time) + longitude


def sun_ecliptic_longitude(model_time):
    """
    Ecliptic longitude of the sun.

    Reference:
        http://www.geoastro.de/elevaz/basics/meeus.htm
    """
    julian_centuries = days_from_2000(model_time) / 36525.0

    # mean anomaly calculation
    mean_anomaly = np.deg2rad(357.52910 +
                              35999.05030*julian_centuries -
                              0.0001559*julian_centuries*julian_centuries -
                              0.00000048*julian_centuries*julian_centuries*julian_centuries)

    # mean longitude
    mean_longitude = np.deg2rad(280.46645 +
                                36000.76983*julian_centuries +
                                0.0003032*(julian_centuries**2))

    d_l = np.deg2rad((1.914600 - 0.004817*julian_centuries -
                      0.000014*(julian_centuries**2))*np.sin(mean_anomaly) +
                     (0.019993 - 0.000101*julian_centuries)*np.sin(2*mean_anomaly) +
                     0.000290*np.sin(3*mean_anomaly))

    # true longitude
    return mean_longitude + d_l


def obliquity_star(julian_centuries):
    """
    return obliquity of the sun
    Use 5th order equation from
    https://en.wikipedia.org/wiki/Ecliptic#Obliquity_of_the_ecliptic
    """
    return np.deg2rad(
        23.0 + 26.0/60 + 21.406/3600.0 -
        (46.836769*julian_centuries - 0.0001831*(julian_centuries**2) +
         0.00200340*(julian_centuries**3) - 0.576e-6*(julian_centuries**4) -
         4.34e-8*(julian_centuries**5))/3600.)


def right_ascension_declination(model_time):
    """
    Right ascension and declination of the sun.
    Ref:
        http://www.geoastro.de/elevaz/basics/meeus.htm
    """
    julian_centuries = days_from_2000(model_time) / 36525.0
    eps = obliquity_star(julian_centuries)

    eclon = sun_ecliptic_longitude(model_time)
    x = np.cos(eclon)
    y = np.cos(eps) * np.sin(eclon)
    z = np.sin(eps) * np.sin(eclon)
    r = np.sqrt(1.0 - z * z)
    # sun declination
    declination = np.arctan2(z, r)
    # right ascension
    right_ascension = 2 * np.arctan2(y, (x + r))
    return right_ascension, declination


def local_hour_angle(model_time, longitude, right_ascension):
    """
    Hour angle at model_time for the given longitude and right_ascension
    longitude in radians

    Ref:
        https://en.wikipedia.org/wiki/Hour_angle#Relation_with_the_right_ascension
    """
    return local_mean_sidereal_time(model_time, longitude) - right_ascension


def star_zenith_azimuth(model_time, lon, lat):
    """
    Return star Zenith and azimuth
    lon,lat in radians
    Ref:
        Azimuth:
            https://en.wikipedia.org/wiki/Solar_azimuth_angle#Formulas
        Zenith:
            https://en.wikipedia.org/wiki/Solar_zenith_angle

    """

    ra, dec = right_ascension_declination(model_time)
    h_angle = local_hour_angle(model_time, lon, ra)

    zenith = np.arccos(np.sin(lat)*np.sin(dec) +
                       np.cos(lat) * np.cos(dec) * np.cos(h_angle))

    azimuth = np.arctan2(-np.sin(h_angle), (np.cos(lat)*np.tan(dec) -
                                            np.sin(lat)*np.cos(h_angle)))

    return zenith, azimuth


def sun_zenith_angle(model_time, lon, lat):
    """
    Sun-zenith angle for lon, lat at model_time.
    lon,lat in radians.
    The angle returned is in radians
    """
    zenith, azimuth = star_zenith_azimuth(model_time, lon, lat)
    return zenith

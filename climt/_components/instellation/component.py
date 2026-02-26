from sympl import DiagnosticComponent
import datetime
import numpy as np
from ..._core.backend import get_array_namespace, jit_compile, prange

class Instellation(DiagnosticComponent):
    """
    Calculates the zenith angle and star-planet correction
    factor given orbital parameters. Currently useful only
    for Earth-sun system.
    """

    input_properties = {
        "latitude": {
            "dims": ["*"],
            "units": "degrees_north",
        },
        "longitude": {
            "dims": ["*"],
            "units": "degrees_east",
        },
    }

    diagnostic_properties = {
        "zenith_angle": {
            "dims": ["*"],
            "units": "radians",
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
        lat = state["latitude"]
        lon = state["longitude"]
        xp = get_array_namespace(lat)
        
        # Flatten inputs
        lat_flat = xp.reshape(lat, (-1,))
        lon_flat = xp.reshape(lon, (-1,))
        
        # All time-based params are scalar for a single array_call
        julian_centuries = days_from_2000(state["time"]) / 36525.0
        fractional_day_val = fractional_day(state["time"])
        
        if xp is np:
            zen_angle = _instellation_kernel_np(lat_flat, lon_flat, julian_centuries, fractional_day_val)
        else:
            zen_angle = _instellation_kernel_jax(lat_flat, lon_flat, julian_centuries, fractional_day_val)
            
        return {"zenith_angle": xp.reshape(zen_angle, lat.shape)}


def days_from_2000(model_time):
    """Get the days since year 2000."""
    return total_days(model_time - datetime.datetime(2000, 1, 1, 12, 0))


def total_days(time_diff):
    """
    Total time in units of days
    """
    return time_diff.days + (
        time_diff.seconds + time_diff.microseconds / (1000000.0)
    ) / (24 * 3600.0)

def fractional_day(dt):
    day_start = type(dt)(dt.year, dt.month, dt.day)
    return (dt - day_start).total_seconds() / (24.0 * 60.0 * 60.0)

@jit_compile(backend=np, parallel=True)
def _instellation_kernel_np(lat_deg, lon_deg, julian_centuries, frac_day):
    ncol = lat_deg.size
    zenith = np.zeros(ncol)
    
    # 1. Earth orbit params
    eps = _obliquity_star_jit(julian_centuries)
    eclon = _sun_ecliptic_longitude_jit(julian_centuries)
    
    # 2. Right ascension / declination
    x = np.cos(eclon)
    y = np.cos(eps) * np.sin(eclon)
    z = np.sin(eps) * np.sin(eclon)
    r = np.sqrt(1.0 - z * z)
    declination = np.arctan2(z, r)
    right_ascension = 2.0 * np.arctan2(y, (x + r))
    
    # 3. Greenwich sidereal time
    gmst = _gmst_jit(julian_centuries)
    
    sin_lat = np.sin(np.deg2rad(lat_deg))
    cos_lat = np.cos(np.deg2rad(lat_deg))
    sin_dec = np.sin(declination)
    cos_dec = np.cos(declination)
    
    pi2 = 2.0 * np.pi
    deg_to_rad = np.pi / 180.0
    
    for i in prange(ncol):
        lmst = gmst + lon_deg[i] * deg_to_rad
        h_angle = lmst - right_ascension
        
        cos_mu = sin_lat[i] * sin_dec + cos_lat[i] * cos_dec * np.cos(h_angle)
        
        # Clamp cos_mu to [-1, 1]
        if cos_mu > 1.0: cos_mu = 1.0
        elif cos_mu < -1.0: cos_mu = -1.0
        
        z_angle = np.arccos(cos_mu)
        
        # Clamp zenith angle to [-PI/2, PI/2]
        if z_angle > np.pi / 2.0: z_angle = np.pi / 2.0
        elif z_angle < -np.pi / 2.0: z_angle = -np.pi / 2.0
        
        zenith[i] = z_angle
        
    return zenith

def _instellation_kernel_jax(lat_deg, lon_deg, julian_centuries, frac_day):
    import jax.numpy as jnp
    
    eps = _obliquity_star_jit(julian_centuries)
    eclon = _sun_ecliptic_longitude_jit(julian_centuries)
    
    x = jnp.cos(eclon); y = jnp.cos(eps) * jnp.sin(eclon)
    z = jnp.sin(eps) * jnp.sin(eclon); r = jnp.sqrt(1.0 - z * z)
    declination = jnp.arctan2(z, r); right_ascension = 2.0 * jnp.arctan2(y, (x + r))
    
    gmst = _gmst_jit(julian_centuries)
    
    deg_to_rad = jnp.pi / 180.0
    lat_rad = lat_deg * deg_to_rad
    lon_rad = lon_deg * deg_to_rad
    
    lmst = gmst + lon_rad
    h_angle = lmst - right_ascension
    
    cos_mu = jnp.sin(lat_rad) * jnp.sin(declination) + jnp.cos(lat_rad) * jnp.cos(declination) * jnp.cos(h_angle)
    zenith = jnp.arccos(jnp.clip(cos_mu, -1.0, 1.0))
    
    return jnp.clip(zenith, -jnp.pi/2.0, jnp.pi/2.0)

@jit_compile
def _obliquity_star_jit(julian_centuries):
    return np.deg2rad(
        23.0
        + 26.0 / 60
        + 21.406 / 3600.0
        - (
            46.836769 * julian_centuries
            - 0.0001831 * (julian_centuries**2)
            + 0.00200340 * (julian_centuries**3)
            - 0.576e-6 * (julian_centuries**4)
            - 4.34e-8 * (julian_centuries**5)
        )
        / 3600.0
    )

@jit_compile
def _sun_ecliptic_longitude_jit(julian_centuries):
    # mean anomaly calculation
    mean_anomaly = np.deg2rad(
        357.52910
        + 35999.05030 * julian_centuries
        - 0.0001559 * julian_centuries * julian_centuries
        - 0.00000048 * julian_centuries * julian_centuries * julian_centuries
    )

    # mean longitude
    mean_longitude = np.deg2rad(
        280.46645 + 36000.76983 * julian_centuries + 0.0003032 * (julian_centuries**2)
    )

    d_l = np.deg2rad(
        (1.914600 - 0.004817 * julian_centuries - 0.000014 * (julian_centuries**2))
        * np.sin(mean_anomaly)
        + (0.019993 - 0.000101 * julian_centuries) * np.sin(2 * mean_anomaly)
        + 0.000290 * np.sin(3 * mean_anomaly)
    )

    # true longitude
    return mean_longitude + d_l

@jit_compile
def _gmst_jit(julian_centuries):
    theta = 67310.54841 + julian_centuries * (
        876600 * 3600
        + 8640184.812866
        + julian_centuries * (0.093104 - julian_centuries * 6.2 * 10e-6)
    )
    theta_radians = np.deg2rad(theta / 240.0) % (2.0 * np.pi)
    # Handle negative result of modulo
    if theta_radians < 0:
        theta_radians += 2.0 * np.pi
    return theta_radians

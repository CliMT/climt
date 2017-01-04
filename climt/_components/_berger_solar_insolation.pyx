from libc.math cimport sin, cos, sqrt, asin, acos, atan, floor
from libc.math cimport M_PI as PI
import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# Abbreviations/term names used are as in Berger 1978. Values are taken from
# CAM 3.0, shr_orb_mod.f90.
# Also useful is http://www.cesm.ucar.edu/models/atm-cam/docs/description/node22.html

# amplitudes for obliquity cos series
cdef DTYPE_t *A = [
    -2462.2214466, -857.3232075, -629.3231835,
    -414.2804924, -311.7632587, 308.9408604,
    -162.5533601, -116.1077911, 101.1189923,
    -67.6856209, 24.9079067, 22.5811241,
    -21.1648355, -15.6549876, 15.3936813,
    14.6660938, -11.7273029, 10.2742696,
    6.4914588, 5.8539148, -5.4872205,
    -5.4290191, 5.1609570, 5.0786314,
    -4.0735782, 3.7227167, 3.3971932,
    -2.8347004, -2.6550721, -2.5717867,
    -2.4712188, 2.4625410, 2.2464112,
    -2.0755511, -1.9713669,-1.8813061,
    -1.8468785, 1.8186742, 1.7601888,
    -1.5428851, 1.4738838, -1.4593669,
    1.4192259, -1.1818980, 1.1756474,
    -1.1316126, 1.0896928
]

cdef int A_length = 47

# rates for obliquity cosine series
cdef DTYPE_t *f = [
    31.609974, 32.620504, 24.172203,
    31.983787, 44.828336, 30.973257,
    43.668246, 32.246691, 30.599444,
    42.681324, 43.836462, 47.439436,
    63.219948, 64.230478, 1.010530,
    7.437771, 55.782177, 0.373813,
    13.218362, 62.583231, 63.593761,
    76.438310, 45.815258, 8.448301,
    56.792707, 49.747842, 12.058272,
    75.278220, 65.241008, 64.604291,
    1.647247, 7.811584, 12.207832,
    63.856665, 56.155990, 77.448840,
    6.801054, 62.209418, 20.656133,
    48.344406, 55.145460, 69.000539,
    11.071350, 74.291298, 11.047742,
    0.636717, 12.844549
]

# phases for obliquity cosine series
cdef DTYPE_t *delta = [
    251.9025, 280.8325, 128.3057,
    292.7252, 15.3747, 263.7951,
    308.4258, 240.0099, 222.9725,
    268.7809, 316.7998, 319.6024,
    143.8050, 172.7351, 28.9300,
    123.5968, 20.2082, 40.8226,
    123.4722, 155.6977, 184.6277,
    267.2772, 55.0196, 152.5268,
    49.1382, 204.6609, 56.5233,
    200.3284, 201.6651, 213.5577,
    17.0374, 164.4194, 94.5422,
    131.9124, 61.0309, 296.2073,
    135.4894, 114.8750, 247.0691,
    256.6114, 32.1008, 143.6804,
    16.8784, 160.6835, 27.5932,
    348.1074, 82.6496
]

# ampl for eccen/fvelp cos/sin series (e cos/sin(PI))
cdef DTYPE_t *P = [
    0.01860798, 0.01627522, -0.01300660,
    0.00988829, -0.00336700, 0.00333077,
    -0.00235400, 0.00140015, 0.00100700,
    0.00085700, 0.00064990, 0.00059900,
    0.00037800, -0.00033700, 0.00027600,
    0.00018200, -0.00017400, -0.00012400,
    0.00001250
]

cdef int P_length=19

# rates for eccen/fvelp cos/sin series (e cos/sin(PI))
cdef DTYPE_t *alpha = [
    4.2072050, 7.3460910, 17.8572630,
    17.2205460, 16.8467330, 5.1990790,
    18.2310760, 26.2167580, 6.3591690,
    16.2100160, 3.0651810, 16.5838290,
    18.4939800, 6.1909530, 18.8677930,
    17.4255670, 6.1860010, 18.4174410,
    0.6678630
]

# phases for eccen/fvelp cos/sin series (e cos/sin(PI))
cdef DTYPE_t *zeta = [
    28.620089, 193.788772, 308.307024,
    320.199637, 279.376984, 87.195000,
    349.129677, 128.443387, 154.143880,
    291.269597, 114.860583, 332.092251,
    296.414411, 145.769910, 337.237063,
    152.092288, 126.839891, 210.667199,
    72.108838,
]

# amplitudes for mvelp sine series
cdef DTYPE_t *F = [
    7391.0225890, 2555.1526947, 2022.7629188,
    -1973.6517951, 1240.2321818, 953.8679112,
    -931.7537108, 872.3795383, 606.3544732,
    -496.0274038, 456.9608039, 346.9462320,
    -305.8412902, 249.6173246, -199.1027200,
    191.0560889, -175.2936572, 165.9068833,
    161.1285917, 139.7878093, -133.5228399,
    117.0673811, 104.6907281, 95.3227476,
    86.7824524, 86.0857729, 70.5893698,
    -69.9719343, -62.5817473, 61.5450059,
    -57.9364011, 57.1899832, -57.0236109,
    -54.2119253, 53.2834147, 52.1223575,
    -49.0059908, -48.3118757, -45.4191685,
    -42.2357920, -34.7971099, 34.4623613,
    -33.8356643, 33.6689362, -31.2521586,
    -30.8798701, 28.4640769, -27.1960802,
    27.0860736, -26.3437456, 24.7253740,
    24.6732126, 24.4272733, 24.0127327,
    21.7150294, -21.5375347, 18.1148363,
    -16.9603104, -16.1765215, 15.5567653,
    15.4846529, 15.2150632, 14.5047426,
    -14.3873316, 13.1351419, 12.8776311,
    11.9867234, 11.9385578, 11.7030822,
    11.6018181, -11.2617293, -10.4664199,
    10.4333970, -10.2377466, 10.1934446,
    -10.1280191, 10.0289441, -10.0034259
]

cdef int F_length=78

# rates for mvelp sine series
cdef DTYPE_t *f_prime = [
    31.609974, 32.620504, 24.172203,
    0.636717, 31.983787, 3.138886,
    30.973257, 44.828336, 0.991874,
    0.373813, 43.668246, 32.246691,
    30.599444, 2.147012, 10.511172,
    42.681324, 13.650058, 0.986922,
    9.874455, 13.013341, 0.262904,
    0.004952, 1.142024, 63.219948,
    0.205021, 2.151964, 64.230478,
    43.836462, 47.439436, 1.384343,
    7.437771, 18.829299, 9.500642,
    0.431696, 1.160090, 55.782177,
    12.639528, 1.155138, 0.168216,
    1.647247, 10.884985, 5.610937,
    12.658184, 1.010530, 1.983748,
    14.023871, 0.560178, 1.273434,
    12.021467, 62.583231, 63.593761,
    76.438310, 4.280910, 13.218362,
    17.818769, 8.359495, 56.792707,
    8.448301, 1.978796, 8.863925,
    0.186365, 8.996212, 6.771027,
    45.815258, 12.002811, 75.278220,
    65.241008, 18.870667, 22.009553,
    64.604291, 11.498094, 0.578834,
    9.237738, 49.747842, 2.147012,
    1.196895, 2.133898, 0.173168
]

# phases for mvelp sine series
cdef DTYPE_t *delta_prime = [
    251.9025, 280.8325, 128.3057,
    348.1074, 292.7252, 165.1686,
    263.7951, 15.3747, 58.5749,
    40.8226, 308.4258, 240.0099,
    222.9725, 106.5937, 114.5182,
    268.7809, 279.6869, 39.6448,
    126.4108, 291.5795, 307.2848,
    18.9300, 273.7596, 143.8050,
    191.8927, 125.5237, 172.7351,
    316.7998, 319.6024, 69.7526,
    123.5968, 217.6432, 85.5882,
    156.2147, 66.9489, 20.2082,
    250.7568, 48.0188, 8.3739,
    17.0374, 155.3409, 94.1709,
    221.1120, 28.9300, 117.1498,
    320.5095, 262.3602, 336.2148,
    233.0046, 155.6977, 184.6277,
    267.2772, 78.9281, 123.4722,
    188.7132, 180.1364, 49.1382,
    152.5268, 98.2198, 97.4808,
    221.5376, 168.2438, 161.1199,
    55.0196, 262.6495, 200.3284,
    201.6651, 294.6547, 99.8233,
    213.5577, 154.1631, 232.7153,
    138.3034, 204.6609, 106.5938,
    250.4676, 332.3345, 27.3039
]

cdef DTYPE_t arcsec_to_degree = 1./3600.


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_orbital_parameters(DTYPE_t years_since_jan_1_1950):
    cdef int i
    cdef DTYPE_t t = years_since_jan_1_1950
    # Equation 1
    cdef DTYPE_t obliquity = 23.320556
    for i in range(A_length):
        obliquity += A[i]*arcsec_to_degree*cos(
            (f[i]*arcsec_to_degree*t + delta[i]) * PI/180.)
    obliquity = obliquity * PI/180.
    # Equations 2-3
    cdef DTYPE_t cos_sum = 0.
    cdef DTYPE_t sin_sum = 0
    for i in range(P_length):
        cos_sum += P[i]*cos(alpha[i]*arcsec_to_degree*t + zeta[i])
        sin_sum += P[i]*sin(alpha[i]*arcsec_to_degree*t + zeta[i])
    cdef DTYPE_t eccentricity_squared = cos_sum*cos_sum + sin_sum*sin_sum
    cdef DTYPE_t eccentricity = sqrt(eccentricity_squared)
    cdef DTYPE_t eccentricity_cubed = eccentricity*eccentricity_squared
    # Equation 4
    cdef DTYPE_t pi = get_fixed_vernal_equinox_longitude_of_perihelion(
        cos_sum, sin_sum)
    # Equation 6
    cdef DTYPE_t omega_tilde = (
        pi * 180./PI +
        50.439273*arcsec_to_degree*t +
        3.392506
    )
    for i in range(F_length):
        omega_tilde += F[i]*sin(
            (f_prime[i]*arcsec_to_degree*t + delta_prime[i]) * PI/180.)
    while omega_tilde < 0.:
        omega_tilde += 360.
    while omega_tilde > 360:
        omega_tilde -= 360.
    omega_tilde = omega_tilde * PI/180.

    cdef DTYPE_t beta = sqrt(1 - eccentricity_squared)

    # Bullet points 1-4 on page 2365
    cdef DTYPE_t lambda_m0 = 2.*(
        (0.5*eccentricity + 0.125*eccentricity_cubed)*(1. + beta)*sin(omega_tilde + PI) -
        0.25*eccentricity_squared*(0.5 + beta)*sin(2*(omega_tilde + PI)) +
        0.125*eccentricity_cubed*(1./3. + beta)*sin(3*(omega_tilde + PI))
    )
    return lambda_m0, eccentricity, omega_tilde, obliquity


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef get_solar_parameters(
        DTYPE_t lambda_m0,
        DTYPE_t eccentricity,
        DTYPE_t omega_tilde,
        DTYPE_t obliquity,
        DTYPE_t years_since_vernal_equinox,
        DTYPE_t fractional_day,
        np.ndarray[DTYPE_t, ndim=1] lat,
        np.ndarray[DTYPE_t, ndim=1] lon,
        DTYPE_t solar_constant):
    """
    Vernal equinox is assumed to occur on March 21. If you give a difference
    from a different date, that date is assumed to be the vernal equinox. The
    date used should however be at 0 UTC.

    Args:
        years_since_march_21_1950:
        years_since_vernal_equinox:
        fractional_day:
        lat:
        lon:
        solar_constant:

    Returns:

    """
    assert lat.dtype == DTYPE and lon.dtype == DTYPE
    cdef np.ndarray[DTYPE_t, ndim=1] solar_zenith_angle = np.zeros([lat.shape[0]], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] solar_insolation = np.zeros([lat.shape[0]], dtype=DTYPE)
    cdef int i_max = lat.size
    cdef int i
    cdef DTYPE_t eccentricity_squared = eccentricity*eccentricity

    cdef DTYPE_t lambda_m = lambda_m0 + years_since_vernal_equinox*2.*PI
    cdef DTYPE_t temp = lambda_m - (omega_tilde + PI)
    cdef DTYPE_t sin_temp = sin(temp)
    cdef DTYPE_t lmbda = lambda_m + eccentricity*(
        2.*sin_temp + eccentricity*(
            1.25*sin(2*temp) + eccentricity * (
                (13./12.)*sin(3*temp) - 0.25*sin_temp)))
    cdef DTYPE_t inverse_rho = (
        (1 + eccentricity*cos(lmbda - (omega_tilde + PI))) / (
        1 - eccentricity_squared))
    cdef DTYPE_t rho = 1./inverse_rho
    cdef DTYPE_t inverse_rho_squared = inverse_rho*inverse_rho
    cdef DTYPE_t solar_declination_angle = asin(sin(obliquity)*sin(lmbda))
    cdef DTYPE_t H  # hour angle of sun during the day, with 0. for noon at 180 degrees
    cdef DTYPE_t cos_mu  # mu is solar zenith angle
    cdef DTYPE_t sin_delta = sin(solar_declination_angle)
    cdef DTYPE_t cos_delta = cos(solar_declination_angle)
    for i in range(i_max):
        H = 2*PI*(fractional_day + lon[i]/360.)
        cos_mu = sin(lat[i])*sin_delta - cos(lat[i])*cos_delta*cos(H)
        solar_zenith_angle[i] = acos(cos_mu)
        solar_insolation[i] = solar_constant * inverse_rho_squared * cos_mu
    return solar_insolation, solar_zenith_angle, obliquity, eccentricity, rho


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef DTYPE_t get_fixed_vernal_equinox_longitude_of_perihelion(
        DTYPE_t cos_sum, DTYPE_t sin_sum):
    cdef DTYPE_t pi
    if abs(cos_sum) < 1e-8:
        if sin_sum == 0.:
            pi = 0.
        elif sin_sum < 0.:
            pi = 1.5*PI
        elif sin_sum > 0.:
            pi = 0.5*PI
    elif cos_sum < 0.:
        pi = atan(sin_sum/cos_sum) + PI
    elif cos_sum > 0.:
        if sin_sum < 0.:
            pi = atan(sin_sum/cos_sum) + 2.*PI
        else:
            pi = atan(sin_sum/cos_sum)
    return pi



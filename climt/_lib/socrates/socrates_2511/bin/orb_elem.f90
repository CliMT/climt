PROGRAM orb_elem

USE realtype_rd, ONLY: RealK
USE coord_transforms, ONLY: CartCoord, rotate_x, rotate_z

IMPLICIT NONE

  REAL(RealK), PARAMETER :: pi = 4 * atan(1.e0_RealK)
  REAL(RealK), PARAMETER :: pio180 = pi/180.e0_RealK
  REAL(RealK), PARAMETER :: pio2 = pi/2.e0_RealK
  REAL(RealK), PARAMETER :: twopi = pi*2.e0_RealK
  REAL(RealK), PARAMETER :: ICRF_oblq = 23.43928_RealK*pio180
  REAL(RealK), PARAMETER :: JCy = 36525.0_RealK

! Orbital and rotational elements for planet with respect to standard
! celestial (ICRF) coordinates (deltas are per Julian day from epoch):
  REAL(RealK) :: epoch ! Epoch (Julian days)
  REAL(RealK) :: a, da ! Semi-major axis (AU)
  REAL(RealK) :: e, de ! Eccentricity
  REAL(RealK) :: I, dI ! Inclination of orbit to ecliptic (degrees)
  REAL(RealK) :: L, dL ! Mean longitude w.r.t. Vernal Equinox (degrees)
  REAL(RealK) :: lph, dlph ! Longitude of perihelion (degrees)
  REAL(RealK) :: lan, dlan ! Longitude of ascending node (degrees)
  REAL(RealK) :: ra_pole, dra_pole ! Right Ascension of North Pole (degrees)
  REAL(RealK) :: dec_pole, ddec_pole ! Declination of North Pole (degrees)
  REAL(RealK) :: W, dW ! Location of planets prime meridian (degrees)

! Some derived orbital elements
  REAL(RealK) :: aph, daph ! Argument of perihelion

! Cartesian coordinates of planet north pole in various reference frames
  TYPE (CartCoord) :: icrf         ! International Celestial Reference Frame (Earth Equatorial)
  TYPE (CartCoord) :: ecliptic     ! Ecliptic coordinate system (Earth orbit, Vernal Point)
  TYPE (CartCoord) :: planet_node  ! Ecliptic plane, x-axis is planet ascending node
  TYPE (CartCoord) :: planet_orbit ! Planet orbit plane, x-axis is planet ascending node
  TYPE (CartCoord) :: planet_peri  ! Planet orbit plane, x-axis is planet perihelion

! Orbital elements of planet w.r.t. planets orbit and planet's equivalent
! of the Vernal Equinox (i.e. position of ascending node)
  REAL(RealK) :: planet_epoch
  REAL(RealK) :: planet_a, planet_da
  REAL(RealK) :: planet_e, planet_de
  REAL(RealK) :: planet_M, planet_dM ! Mean anomaly from perihelion (radians)
  REAL(RealK) :: planet_oblq, planet_doblq ! Obliquity to orbit (radians)
  REAL(RealK) :: planet_lph, planet_dlph ! Longitude of perihelion (radians)
  REAL(RealK) :: planet_ha, planet_dha, planet_dha2 ! Planet hour angle of prime meridian (radians)
  REAL(RealK) :: omega ! Rotation rate (radians per second)

  INTEGER :: hours, minutes
  REAL :: seconds

! https://ssd.jpl.nasa.gov/planets/approx_pos.html
!               a           e            I               L          long.peri.    long.node.
!           AU, AU/Cy     ,    /Cy  deg, deg/Cy     deg, deg/Cy    deg, deg/Cy   deg, deg/Cy
!-------------------------------------------------------------------------------------------
!Mercury   0.38709927   0.20563593   7.00497902     252.25032350   77.45779628   48.33076593
!          0.00000037   0.00001906  -0.00594749  149472.67411175    0.16047689   -0.12534081
!Venus     0.72333566   0.00677672   3.39467605     181.97909950  131.60246718   76.67984255
!          0.00000390  -0.00004107  -0.00078890   58517.81538729    0.00268329   -0.27769418
!EM Bary   1.00000261   0.01671123  -0.00001531     100.46457166  102.93768193    0.0
!          0.00000562  -0.00004392  -0.01294668   35999.37244981    0.32327364    0.0
!Mars      1.52371034   0.09339410   1.84969142      -4.55343205  -23.94362959   49.55953891
!          0.00001847   0.00007882  -0.00813131   19140.30268499    0.44441088   -0.29257343
!Jupiter   5.20288700   0.04838624   1.30439695      34.39644051   14.72847983  100.47390909
!         -0.00011607  -0.00013253  -0.00183714    3034.74612775    0.21252668    0.20469106
!Saturn    9.53667594   0.05386179   2.48599187      49.95424423   92.59887831  113.66242448
!         -0.00125060  -0.00050991   0.00193609    1222.49362201   -0.41897216   -0.28867794
!Uranus   19.18916464   0.04725744   0.77263783     313.23810451  170.95427630   74.01692503
!         -0.00196176  -0.00004397  -0.00242939     428.48202785    0.40805281    0.04240589
!Neptune  30.06992276   0.00859048   1.77004347     -55.12002969   44.96476227  131.78422574
!          0.00026291   0.00005105   0.00035372     218.45945325   -0.32241464   -0.00508664

! The direction of the poles and rotational elements are found in reports of the
! "IAU/IAG Working Group on cartographic coordinates and rotational elements":

! 2006 report:

! Table 1 Recommended values for the direction of the north pole of rotation
! and the prime meridian of the Sun and planets
! RA , DEC are ICRF equatorial coordinates at epoch J2000.0.
! T = interval in Julian centuries (of 36525 days) from the standard epoch
! d = interval in days from the standard epoch.
! The standard epoch is JD 2451545.0, i.e. 2000 January 1 12 hours TDB
! Angles in degrees.

! Sun
! RA = 286.13
! DEC = 63.87
! W = 84.176 + 14.1844000d

! Mercury
! RA = 281.01 - 0.033T
! DEC = 61.45 - 0.005T
! W = 329.548 + 6.1385025d

! Venus
! RA = 272.76
! DEC = 67.16
! W = 160.20 - 1.4813688d

! Earth
! RA = 0.00 - 0.641T
! DEC = 90.00 - 0.557T
! W = 190.147 + 360.9856235d

! Mars
! RA = 317.68143 - 0.1061T
! DEC = 52.88650 - 0.0609T
! W = 176.630 + 350.89198226d

! Jupiter
! RA = 268.056595 - 0.006499T + 0 .000117 sin Ja + 0 .000938 sin Jb
! + 0.001432 sin Jc + 0.000030 sin Jd + 0.002150 sin Je
! DEC = 64.495303 + 0.002413T + 0.000050 cos Ja + 0.000404 cos Jb
! + 0.000617 cos Jc - 0.000013 cos Jd + 0.000926 cos Je
! W = 284.95 + 870.5366420d (e)
! where Ja = 99.360714 + 4850.4046T, Jb = 175.895369 + 1191.9605T,
! Jc = 300.323162 + 262.5475T, Jd = 114.012305 + 6070.2476T,
! Je = 49.511251 + 64.3000T

! Saturn
! RA = 40.589 - 0.036T
! DEC = 83.537 - 0.004T
! W = 38.90 + 810.7939024d

! Uranus
! RA = 257.311
! DEC = -15.175
! W = 203.81 - 501.1600928d

! Neptune
! RA = 299.36 + 0.70 sin N
! DEC = 43.46 - 0.51 cos N
! W = 253.18 + 536.3128492d - 0.48 sin N
! N = 357.85 + 52.316T

! Pluto
! RA = 312.993
! DEC = 6.163
! W = 237.305 - 56.3625225d



! Uncomment planet required:

!! Saturn
!  print*,'Saturn:'
!  epoch = 2451545.0 ! (J2000)
!  a = 9.53667594_RealK
!  da = -0.00125060_RealK
!  e = 0.05386179_RealK
!  de = -0.00050991_RealK
!  I = 2.48599187_RealK*pio180
!  dI = 0.00193609_RealK*pio180
!  L = 49.95424423_RealK*pio180
!  dL = 1222.49362201_RealK*pio180
!  lph = 92.59887831_RealK*pio180
!  dlph = -0.41897216_RealK*pio180
!  lan = 113.66242448_RealK*pio180
!  dlan = -0.28867794_RealK*pio180
!  ra_pole = 40.589_RealK*pio180
!  dra_pole = -0.036_RealK*pio180
!  dec_pole = 83.537_RealK*pio180
!  ddec_pole = -0.004_RealK*pio180
!  W = 38.90_RealK*pio180
!  dW = 810.7939024_RealK*pio180

!! Uranus
!  epoch = 2451545.0 ! (J2000)
!  a = 19.18916464_RealK
!  da = -0.00196176_RealK
!  e = 0.04725744_RealK
!  de = -0.00004397_RealK
!  I = 0.77263783_RealK*pio180
!  dI = -0.00242939_RealK*pio180
!  L = 313.23810451_RealK*pio180
!  dL = 428.48202785_RealK*pio180
!  lph = 170.95427630_RealK*pio180
!  dlph = 0.40805281_RealK*pio180
!  lan = 74.01692503_RealK*pio180
!  dlan = 0.04240589_RealK*pio180
!  ra_pole = 77.311_RealK*pio180
!  dra_pole = 0.0_RealK*pio180
!  dec_pole = 15.175_RealK*pio180
!  ddec_pole = 0.0_RealK*pio180
!  W = 203.81_RealK*pio180
!  dW = 501.1600928_RealK*pio180

!! Earth
!  epoch = 2451545.0_RealK ! (J2000)
!  a = 1.00000261_RealK
!  da = 0.00000562_RealK
!  e = 0.01671123_RealK
!  de = -0.00004392_RealK
!  I = -0.00001531_RealK*pio180
!  dI = -0.01294668_RealK*pio180
!  L = 100.46457166_RealK*pio180
!  dL = 35999.37244981_RealK*pio180
!  lph = 102.93768193_RealK*pio180
!  dlph = 0.32327364_RealK*pio180
!  lan = 0.0_RealK
!  dlan = 0.0_RealK
!  ra_pole = 0.0_RealK
!  dra_pole = -0.641_RealK*pio180
!  dec_pole = 90.00_RealK*pio180
!  ddec_pole = -0.557_RealK*pio180
!  W = 190.147_RealK*pio180
!  dW = 360.9856235_RealK*pio180

! Mars
! Note, see: A post-Pathfinder evaluation of aerocentric solar coordinates
! with improved timing recipes for Mars seasonal/diurnal climate studies
! (http://pubs.giss.nasa.gov/docs/2000/2000_Allison_McEwen_1.pdf)
  print*,'Mars:'
  epoch = 2451545.0_RealK ! (J2000)
  a = 1.52368_RealK
  da = 0.0_RealK
  e = 0.09340_RealK
  de = 2.477E-09_RealK*JCy
  I = 1.8497_RealK*pio180
  dI = -2.23E-07_RealK*pio180*JCy
  lph = (336.0602_RealK)*pio180
  dlph = 1.2153E-05_RealK*pio180*JCy
  L = 19.3870_RealK*pio180 + lph
  dL = 0.52402075_RealK*pio180*JCy + dlph
  lan = 49.5581_RealK*pio180
  dlan = -8.077E-06_RealK*pio180*JCy
  ra_pole = 317.68143_RealK*pio180
  dra_pole = -0.1061_RealK*pio180
  dec_pole = 52.88650_RealK*pio180
  ddec_pole = -0.0609_RealK*pio180
  W = 176.630_RealK*pio180
  dW = 350.89198226_RealK*pio180

!! Jupiter
!  epoch = 2451545.0 ! (J2000)
!  a = 5.20288700_RealK
!  da = -0.00011607_RealK
!  e = 0.04838624_RealK
!  de = -0.00013253_RealK
!  I = 1.30439695_RealK*pio180
!  dI = -0.00183714_RealK*pio180
!  L = 34.39644051_RealK*pio180
!  dL = 3034.74612775_RealK*pio180
!  lph = 14.72847983_RealK*pio180
!  dlph = 0.21252668_RealK*pio180
!  lan = 100.47390909_RealK*pio180
!  dlan = 0.20469106_RealK*pio180
!  ra_pole = 268.056595_RealK*pio180
!  dra_pole = -0.006499_RealK*pio180
!  dec_pole = 64.495303_RealK*pio180
!  ddec_pole = 0.002413_RealK*pio180
!  W = 284.95_RealK*pio180
!  dW = 870.536_RealK*pio180


! Some elements are identical or simply calculated:
  planet_epoch = epoch
  planet_a = a
  planet_da = da/JCy
  planet_e = e
  planet_de = de/JCy
  planet_M = L - lph
  planet_dM = (dL - dlph)/JCy
  omega = dW/86400.0_RealK
  aph = lph - lan

! Other elements require transformation from Earth's equatorial international
! celestial reference frame (ICRF) to the planets equatorial reference frame:

  icrf%x = COS(dec_pole)*COS(ra_pole)
  icrf%y = COS(dec_pole)*SIN(ra_pole)
  icrf%z = SIN(dec_pole)

  CALL rotate_x(ICRF_oblq, icrf,         ecliptic    )
  CALL rotate_z(lan,       ecliptic,     planet_node )
  CALL rotate_x(I,         planet_node,  planet_orbit)
  CALL rotate_z(aph,       planet_orbit, planet_peri )

  planet_oblq = ACOS(planet_peri%z)
  planet_lph = ATAN2(planet_peri%x, planet_peri%y)

  daph = aph + dlph - dlan
  dra_pole=ra_pole+dra_pole
  ddec_pole=dec_pole+ddec_pole
  dI=I+dI
  dlan=lan+dlan

  icrf%x = COS(ddec_pole)*COS(dra_pole)
  icrf%y = COS(ddec_pole)*SIN(dra_pole)
  icrf%z = SIN(ddec_pole)

  CALL rotate_x(ICRF_oblq, icrf,         ecliptic    )
  CALL rotate_z(dlan,      ecliptic,     planet_node )
  CALL rotate_x(dI,        planet_node,  planet_orbit)
  CALL rotate_z(daph,      planet_orbit, planet_peri )

  planet_doblq = (ACOS(planet_peri%z) - planet_oblq)/JCy
  planet_dlph = (ATAN2(planet_peri%x, planet_peri%y) - planet_lph)/JCy

! Planet hour angle at epoch only matters if we have longitudinal asymmetry
  planet_ha = 0.0

! approximate (probably good enough):
  planet_dha = dW - planet_dM - planet_dlph

hours = INT(twopi*24.0/planet_dha)
minutes = INT((twopi*24.0/planet_dha - hours)*60.0)
seconds = ((twopi*24.0/planet_dha - hours)*60.0 - minutes)*60.0

print*,'---Readable units (degrees etc)---'
print*,'planet_epoch = ',planet_epoch
print*,'planet_a, planet_da = ', planet_a, planet_da
print*,'planet_e, planet_de = ', planet_e, planet_de
print*,'planet_oblq, planet_doblq = ', planet_oblq/pio180, planet_doblq/pio180
print*,'planet_lph, planet_dlph = ', planet_lph/pio180, planet_dlph/pio180
print*,'planet_M, planet_dM = ', planet_M/pio180, planet_dM/pio180
print*,'planet_ha, planet_dha = ', planet_ha/pio180, planet_dha/pio180
print*,'omega = ', omega
print*,'Length of day =', hours,' hours', minutes,' minutes', seconds,' secs'
print*,''
print*,'---Units for UM (radians etc)-----'
print*,'planet_epoch = ',planet_epoch
print*,'planet_a, planet_da = ', planet_a, planet_da
print*,'planet_e, planet_de = ', planet_e, planet_de
print*,'planet_oblq, planet_doblq = ', planet_oblq, planet_doblq
print*,'planet_lph, planet_dlph = ', planet_lph, planet_dlph
print*,'planet_M, planet_dM = ', planet_M, planet_dM
print*,'planet_ha, planet_dha = ', planet_ha, planet_dha
print*,'omega = ', omega

END PROGRAM orb_elem

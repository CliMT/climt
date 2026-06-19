! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE astro_constants_mod

! Description: Parameters of the Earth's orbit.
!
!     Default values of the orbital elements
!      (currently those for the epoch J2000 which is 1.5d Jan. 2000):
!     The Eccentricity and Longitue of perhelion are recommended by NASA
!      see (http://ssd.jpl.nasa.gov/elem_planets.html)
!     The Obliquity value comes from the Astronomical Almanac for 1984
!      page S26 and is used on several webpages e.g.
!      nedwww.ipac.caltech.edu/help/calc_doc.txt
!      www.stargazing.net/kepler/astrovba2.html
!      http://edhs1.gsfc.nasa.gov/waisdata/docsw/txt/tp4450505.txt
!
!     The data in the series expansions are adapted from Berger 1978.
!      Andre' L. Berger Journ. of Atm. Sci. Volume 35 p. 2362-2367,
!      and is available from the Met Office library.

USE realtype_rd, ONLY: RealK
USE rad_ccf, ONLY: pi

IMPLICIT NONE

! Eccentricity of the orbit
REAL(RealK), PARAMETER :: E_DFLT = 1.6710222E-02_RealK

! Longitude of the perihelion in radians
REAL(RealK), PARAMETER :: LPH_DFLT = 102.94719_RealK * pi / 180.0_RealK

! Obliquity of the orbit - corresponds to 23.43929111 degrees
REAL(RealK), PARAMETER :: OBLQ_DFLT = 0.409092804_RealK

! Reference year for setting the date of the vernal equinox
INTEGER, PARAMETER :: YEAR_REF_VE = 2000

! Date of the vernal equinox in days after the start of the year
! This date is for the year 2000.
REAL(RealK), PARAMETER :: DATE_VE_DFLT = 79.3159_RealK

!     The final parameter required is the time of the perihelion
!     passage, TAU0. For a pure Keplerian orbit, with a specified
!     eccentricity and longitude of the perihelion, this can be
!     deduced from the date of the vernal equinox (as is specified
!     in AMIP-2, for example). In practice it is somewhat more
!     complicated to calculate the time of the perihelion.
!     For simplicity, a mean value for the years 1995-2005 is used
!     here: note that the range of TAU0 in this period is from
!     1.0 to 3.75 and that there is no simple relationship with leap
!     years.

! Time of the perihelion passage in days
REAL(RealK), PARAMETER :: TAU0_DFLT = 2.667_RealK

!     The tropical year is defined as the mean interval between two
!     successive passages of the sun through the vernal equinox.
!     This value is 365.2424 and will be used here.
!     The value 365.2422 most often used is the mean length of
!     the year starting at different points on the ellipse.

! Number of days in the tropical year
REAL(RealK), PARAMETER :: TropYearLength = 365.2424_RealK

!     ------------------------------------------------------------------

!     The parameters used to calculate secular variations of the orbital
!     elements are taken from A. L. Berger (1978), J. Atm. Sci., vol 35,
!     p. 2362.These have been converted so that:
!     amplitudes are in radians,
!     angular frequencies in radians per year and
!     phases in radians
!     Phases have also been converted to be taken relative to
!     J2000 for consistency with the time for the default values of the
!     orbital parameters.
!     Berger's numbers (with time correction) differ slightly from the
!     default values above.
!
!     The obliquity and longitude of the perihelion have been adjusted
!     to agree with the NASA values for J2000, but adjustment of the
!     eccentricity to NASA values is not so easy and has not been done.


! Reference year
INTEGER, PARAMETER :: YEAR_REF = 2000

!  -----------------------------------------------------------------
!     Obliquity (Table 1) from the first 24 terms:
!     (enough for deviations of less than 0.002 degrees)

! Constant term in the obliquity: from the Astr. Almanac for 1984
! The following value corresponds to 23.320870 degrees at J2000
REAL(RealK), PARAMETER :: OBLQ_CNST = 0.40702597_RealK
!
! Number of terms retained in the series for the obliquity
INTEGER, PARAMETER :: N_TERM_OBQ = 24
!
! Amplitude
REAL(RealK), PARAMETER :: A(N_TERM_OBQ) = (/ &
  -1.19372E-02_RealK,-4.15640E-03_RealK,-3.05103E-03_RealK,-2.00849E-03_RealK, &
  -1.51146E-03_RealK, 1.49778E-03_RealK,-7.88065E-04_RealK,-5.62917E-04_RealK, &
   4.90244E-04_RealK,-3.28170E-04_RealK, 1.20767E-04_RealK, 1.09471E-04_RealK, &
  -1.02587E-04_RealK,-7.58733E-05_RealK, 7.46128E-05_RealK, 7.11222E-05_RealK, &
  -5.68686E-05_RealK, 4.97904E-05_RealK, 3.14644E-05_RealK, 2.83616E-05_RealK, &
  -2.66163E-05_RealK,-2.63254E-05_RealK, 2.50164E-05_RealK, 2.46285E-05_RealK /)
! Angular frequency
REAL(RealK), PARAMETER :: F(N_TERM_OBQ) = (/ &
  1.5324946E-04_RealK,  1.5814864E-04_RealK,  1.1719011E-04_RealK, &
  1.5506174E-04_RealK,  2.1733392E-04_RealK,  1.5016256E-04_RealK, &
  2.1170962E-04_RealK,  1.5633636E-04_RealK,  1.4835028E-04_RealK, &
  2.0692488E-04_RealK,  2.1252514E-04_RealK,  2.2999289E-04_RealK, &
  3.0649899E-04_RealK,  3.1139817E-04_RealK,  4.8991877E-06_RealK, &
  3.6059331E-05_RealK,  2.7043965E-04_RealK,  1.8122966E-06_RealK, &
  6.4084427E-05_RealK,  3.0341210E-04_RealK,  3.0831127E-04_RealK, &
  3.7058338E-04_RealK,  2.2211866E-04_RealK,  4.0958519E-05_RealK /)
! Phase in the series
REAL(RealK), PARAMETER :: D(N_TERM_OBQ) = (/ &
  4.4041E+00_RealK,  4.9093E+00_RealK,  2.2451E+00_RealK,  5.1167E+00_RealK, &
  2.7912E-01_RealK,  4.6115E+00_RealK,  5.3935E+00_RealK,  4.1966E+00_RealK, &
  3.8990E+00_RealK,  4.7014E+00_RealK,  5.5397E+00_RealK,  5.5896E+00_RealK, &
  2.5251E+00_RealK,  3.0303E+00_RealK,  5.0517E-01_RealK,  2.1589E+00_RealK, &
  3.6608E-01_RealK,  7.1253E-01_RealK,  2.1582E+00_RealK,  2.7325E+00_RealK, &
  3.2376E+00_RealK,  4.6833E+00_RealK,  9.7121E-01_RealK,  2.6640E+00_RealK /)

!   -----------------------------------------------------------------
!     Eccentricity and longitude of the fixed perihelion (Table 4):

! Number of terms retained in the series for the
! eccentricty and longitude of the perihelion
INTEGER, PARAMETER :: N_TERM_ECN_LPH = 19
!
! Amplitude
REAL(RealK), PARAMETER :: M(N_TERM_ECN_LPH)= (/ &
  1.8607980E-02_RealK,  1.6275220E-02_RealK, -1.3006600E-02_RealK, &
  9.8882900E-03_RealK, -3.3670000E-03_RealK,  3.3307700E-03_RealK, &
  2.3540000E-03_RealK,  1.4001500E-03_RealK,  1.0070000E-03_RealK, &
  8.5700000E-04_RealK,  6.4990000E-04_RealK,  5.9900000E-04_RealK, &
  3.7800000E-04_RealK, -3.3700000E-04_RealK,  2.7600000E-04_RealK, &
  1.8200000E-04_RealK, -1.7400000E-04_RealK, -1.2400000E-04_RealK, &
  1.2500000E-05_RealK /)
! Angular frequency
REAL(RealK), PARAMETER :: G(N_TERM_ECN_LPH)= (/ &
  2.0397105E-05_RealK,  3.5614854E-05_RealK,  8.6574454E-05_RealK, &
  8.3487563E-05_RealK,  8.1675266E-05_RealK,  2.5205846E-05_RealK, &
  8.8386751E-05_RealK,  1.2710243E-04_RealK,  3.0830121E-05_RealK, &
  7.8588375E-05_RealK,  1.4860417E-05_RealK,  8.0400672E-05_RealK, &
  8.9661345E-05_RealK,  3.0014587E-05_RealK,  9.1473642E-05_RealK, &
  8.4481533E-05_RealK,  2.9990579E-05_RealK,  8.9290274E-05_RealK, &
  3.2378912E-06_RealK /)
! Phase in the series
REAL(RealK), PARAMETER :: B(N_TERM_ECN_LPH) = (/ &
  5.0053E-01_RealK,  3.3839E+00_RealK,  5.3852E+00_RealK,  5.5925E+00_RealK, &
  4.8800E+00_RealK,  1.5230E+00_RealK,  6.0977E+00_RealK,  2.2481E+00_RealK, &
  2.6918E+00_RealK,  5.0874E+00_RealK,  2.0054E+00_RealK,  5.8001E+00_RealK, &
  5.1778E+00_RealK,  2.5455E+00_RealK,  5.8903E+00_RealK,  2.6587E+00_RealK, &
  2.2151E+00_RealK,  3.6812E+00_RealK,  1.2585E+00_RealK /)
!   ------------------------------------------------------------------
!     General Precession (Table 5):

! Linear rate of precession!
! The value corresponds to 50.439273 seconds per year -Berger 1979
REAL(RealK), PARAMETER :: LIN_RATE_GN_PRCS = 2.44536496E-04_RealK

! Constant offset to general precession (in seconds pre year),
! corrected for 50 years difference in reference time.
REAL(RealK), PARAMETER :: GN_PRCS_CNST = 7.14372244E-02_RealK

! Number of terms kept in the series for the general precession
INTEGER, PARAMETER :: N_TERM_GN_PRCS = 10

! Amplitude
REAL(RealK), PARAMETER :: C(N_TERM_GN_PRCS) = (/ &
  3.58327E-02_RealK, 1.23877E-02_RealK, 9.80662E-03_RealK,-9.56853E-03_RealK, &
  6.01280E-03_RealK, 4.62449E-03_RealK,-4.51725E-03_RealK, 4.22942E-03_RealK, &
  2.93967E-03_RealK,-2.40482E-03_RealK /)
! Angular frequency
REAL(RealK), PARAMETER :: H(N_TERM_GN_PRCS) = (/ &
  1.5324946E-04_RealK,  1.5814864E-04_RealK,  1.1719011E-04_RealK,  &
  3.0868911E-06_RealK,  1.5506174E-04_RealK,  1.5217749E-05_RealK,  &
  1.5016256E-04_RealK,  2.1733392E-04_RealK,  4.8087409E-06_RealK,  &
  1.8122966E-06_RealK /)
! Phase in the series
REAL(RealK), PARAMETER :: R(N_TERM_GN_PRCS) = (/ &
  4.4041E+00_RealK,  4.9093E+00_RealK,  2.2451E+00_RealK,  6.0756E+00_RealK, &
  5.1167E+00_RealK,  2.8833E+00_RealK,  4.6115E+00_RealK,  2.7912E-01_RealK, &
  1.0225E+00_RealK,  7.1253E-01_RealK /)
!

END MODULE astro_constants_mod

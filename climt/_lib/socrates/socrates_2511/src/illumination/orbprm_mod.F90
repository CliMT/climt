! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
MODULE orbprm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ORBPRM_MOD'
CONTAINS

SUBROUTINE orbprm(day, year, l_sec_var, lcal360, e, gamph, oblq, mean_anomaly)

! Subroutine to calculate the parameters of the Earth's orbit.
!
! Purpose:
!  This routine returns the parameters of the Earth's orbit:
!  the eccentricity, obliquity, supplement of the longitude of
!  perihelion and the mean anomaly at 12Z on the current day.
!
! Method:
!  For long runs there may be an interest in running with secular
!   variations in the astronomy. The orbital constants have
!   been derived from A. L. Berger 1978, J. Atm. Sci, Volume 35
!   2362-2367. A copy of which can be found in the Met Office library.
!  For short current runs, or long control runs it is preferrable
!   not to allow the astronomy to vary, so fixed values are used.

USE realtype_rd, ONLY: RealK
USE rad_ccf, ONLY: pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE astro_constants_mod , ONLY:                                    &
    a, b, c, d, date_ve_dflt, e_dflt, f, g, gn_prcs_cnst, h,       &
    lin_rate_gn_prcs, lph_dflt, m, n_term_ecn_lph, n_term_gn_prcs, &
    n_term_obq, oblq_cnst, oblq_dflt, r, tau0_dflt, year_ref,      &
    year_ref_ve, tropyearlength

IMPLICIT NONE

INTEGER, INTENT(IN) :: day        ! Day-number in the year
INTEGER, INTENT(IN) :: year       ! Calendar year
LOGICAL, INTENT(IN) :: l_sec_var  ! Include secular variations of the orbit
LOGICAL, INTENT(IN) :: lcal360    ! Use a calendar of 360 days

! Parameters of the Earth's orbit:
REAL(RealK), INTENT(OUT) :: e     ! Eccentricity of the orbit
REAL(RealK), INTENT(OUT) :: gamph ! Supplement of the longitude of perihelion
REAL(RealK), INTENT(OUT) :: oblq  ! Obliquity of the orbit
REAL(RealK), INTENT(OUT) :: mean_anomaly ! Mean anomaly at 12Z


! Local Variables
REAL(RealK) :: tau0          ! Time of the perihelion passage in days
REAL(RealK) :: diny          ! Length of the calendar year (in whole days)
REAL(RealK) :: year_offset   ! Offset of the year from the reference year
                             !   when default values apply
REAL(RealK) :: ecn_sn        ! Eccentricity multiplied by
                             !   the sine of the longitude of the perihelion
REAL(RealK) :: ecn_cn        ! Eccentricity multiplied by
                             !   the cosine of the longitude of the perihelion
REAL(RealK) :: lph_fixed_ve  ! Longitude of the perihelion
                             !   relative to a fixed vernal equinox
REAL(RealK) :: gn_prcs       ! General precession
REAL(RealK) :: date_ve       ! Date of the vernal equinox in days into the year
REAL(RealK) :: no_leap_days  ! The number of leap days, used to calculate DATE_VE
REAL(RealK) :: mean_anom_ve  ! Mean anomaly at the vernal equinox

! Synthetic constants
REAL(RealK) :: beta
REAL(RealK) :: ee1
REAL(RealK) :: ee2
REAL(RealK) :: ee3

INTEGER :: i  ! Loop variable

! Mathematical constants:
REAL(RealK), PARAMETER :: twopi = 2.0_RealK * pi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ORBPRM'

! Astronomical Parameters: from astro_constants_mod


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The length of the calendar year may be set for a 360-day calendar
!  (as is often used in climate runs),
!  or for a real Gregorian calendar which has 365 days in
!  non-leap years and 366 in leap years.
IF (lcal360) THEN
  diny=360.0_RealK
ELSE
  ! Is this a leap year?
  IF (MOD(year,4)    ==  0 .AND.                                  &
     (MOD(year,400)  ==  0 .OR. MOD(year,100)  /=  0)) THEN
    diny = 366.0_RealK
  ! Is this a normal year?
  ELSE
    diny = 365.0_RealK
  END IF
END IF

! The orbital elements are normally set to default values, but
! secular variations may be required in some longer climate runs.

IF (l_sec_var) THEN

  year_offset = REAL( year - year_ref, RealK )

  ! Obliquity: (Equation 1 from Berger 1978)
  oblq = oblq_cnst
  DO i=1, n_term_obq
    oblq = oblq+a(i)*COS(f(i)*year_offset+d(i))
  END DO

  ! Eccentricity: this is better computed from its components
  ! than directly from the series.(Equation (4) of Berger 1978).
  ecn_sn = m(1) * SIN (g(1) * year_offset + b(1))
  ecn_cn = m(1) * COS (g(1) * year_offset + b(1))

  DO i=2, n_term_ecn_lph
    ecn_sn = ecn_sn + m(i) * SIN (g(i) * year_offset + b(i))
    ecn_cn = ecn_cn + m(i) * COS (g(i) * year_offset + b(i))
  END DO
  e = SQRT(ecn_sn*ecn_sn+ecn_cn*ecn_cn)

  ! We now obtain the longitude of the perihelion relative to the
  ! fixed equinox.
  lph_fixed_ve = ATAN2 (ecn_sn,ecn_cn)

  ! The longitude of perihelion and
  !  the supplement of the longitude of the perihelion relative to
  !  the actual vernal equinox requires the general precession.

  ! General Precession.
  gn_prcs = lin_rate_gn_prcs * year_offset + gn_prcs_cnst
  DO i=1, n_term_gn_prcs
    gn_prcs = gn_prcs + c(i) * SIN (h(i) * year_offset + r(i))
  END DO

  ! Supplement of the longitude of the perihelion
  gamph = pi - lph_fixed_ve - gn_prcs

  ! Time of perihelion: The time at which an object is at perihelion
  !  (its closest distance to the sun).
  ! The time of perihelion is inferred from the date of
  !  the vernal equinox using the Gregorian calendar.

  ! Calculate the date of the vernal equinox.
  !  First we need to:
  !   Calculate the no of leap days between year & year_ref_ve.
  !   This needs to be corrected when using the Gregorian calendar.
  !    by adding (DINY-366.0) when the year_ref_ve is a leap year or
  !    by adding (DINY-365.0) when the year_ref_ve is a normal year.
  !   This correction is done when the DATE_VE is calculated below!

  !   In the calculation of NO_LEAP_DAYS below, the divisions of type
  !    'YEAR'/x (where x is 4, 100 or 400) are integer computations.
  !    These integers are then subtracted and the resulting integer
  !    is then converted to a real.

  no_leap_days = ( tropyearlength - 365.0_RealK)                  &
    * REAL( year     - year_ref_ve,     RealK )                   &
    - REAL( year/4   - year_ref_ve/4,   RealK )                   &
    + REAL( year/100 - year_ref_ve/100, RealK )                   &
    - REAL( year/400 - year_ref_ve/400, RealK )

  !  Now we can calculate the date of the vernal equinox!
  !  Because the date of the vernal equinox is varying with the year,
  !  we have to keep track of its position in the sky.
  !  In order to accomodate a time varying vernal equinox when using
  !  a 360-day year, we still have to calculate the difference in
  !  the vernal equinox depending on leap years, normal years and the
  !  difference between the length of the tropical year and the
  !  "normal" year and then we adjust this by multiplying the
  !  DATE_VE by 360/(length of tropical year).

  ! Is a 360 day calendar being used?
  IF (lcal360) THEN
    date_ve = date_ve_dflt + no_leap_days
    date_ve = date_ve * diny / tropyearlength
  ! Is a 365 day calendar being used?
  ELSE
    ! Is the epoch reference year a leap year?
    IF (MOD(year_ref_ve,4)    ==  0 .AND.                         &
       (MOD(year_ref_ve,400)  ==  0 .OR.                          &
        MOD(year_ref_ve,100)  /=  0)) THEN
      date_ve = date_ve_dflt + (no_leap_days + (diny - 366.0_RealK))

    ! Is the epoch reference year a normal year?
    ELSE
      date_ve = date_ve_dflt + (no_leap_days + (diny - 365.0_RealK))
    END IF
  END IF

  beta = SQRT(1.0e+00_RealK-e*e)
  ee1  = (0.5_RealK*e + 0.125_RealK*e*e*e)*(1.0_RealK + beta)
  ee2  = -0.25_RealK*e*e* (0.5_RealK + beta)
  ee3  = 0.125_RealK*e*e*e*((1.0_RealK/3.0_RealK) + beta)
  mean_anom_ve = gamph - 2.0e+00_RealK * (                        &
      ee1 * SIN (gamph)                                           &
    + ee2 * SIN (2.0_RealK * gamph)                               &
    + ee3 * SIN (3.0_RealK * gamph)                               &
    )

  tau0 = date_ve - mean_anom_ve * tropyearlength/(twopi)

ELSE

  e     = e_dflt
  oblq  = oblq_dflt
  gamph = pi - lph_dflt
  tau0  = tau0_dflt

END IF

! If using a 360-day calendar the time of the perihelion is adjusted.
IF (lcal360) THEN
  tau0 = tau0*(360.0_RealK/tropyearlength)+0.71_RealK
END IF

! Calculate the mean anomaly at 12Z on the current day.
! The 0.5 accounts for the time in days since mid-night.
! The references are to Smart 1944 (and UMDP23)
! Eq 67 p. 113 and n=2pi/orbital period     (Eq 3.1.1)
IF (lcal360) THEN
  mean_anomaly = (twopi/diny)           * (REAL(day, RealK) - tau0 - 0.5_RealK)
ELSE
  mean_anomaly = (twopi/tropyearlength) * (REAL(day, RealK) - tau0 - 0.5_RealK)
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE orbprm
END MODULE orbprm_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE solpos_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLPOS_MOD'
CONTAINS

SUBROUTINE solpos (orbit, year, day, seconds, timestep, sindec, eqt, scs, &
                   l_observer, l_calendar_360, sindec_obs, eqt_obs, phase_obs)

! Description:
!  Calculations of the planet's orbit from the day of the year and the
!  orbital elements. It calculates the sin of the solar declination and the
!  inverse-square scaling factor for the solar "constant". It is thus
!  intrinsically scalar. Calculations depend on whether l_cal360 is set:
!  this replaces the Julian calendar with the climate-mode 360-day calendar.

USE realtype_rd, ONLY: RealK
USE rad_ccf, ONLY: pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE def_orbit, ONLY: StrOrbit, ip_smart, ip_mueller, ip_elements_user, &
  ip_elements_earth_fixed, ip_elements_earth_secular_variation, &
  ip_spin_user, ip_spin_earth_day, ip_spin_fixed_sun
USE orbprm_mod, ONLY: orbprm

IMPLICIT NONE

TYPE(StrOrbit), INTENT(IN) :: orbit ! Orbital elements

INTEGER, INTENT(IN)     :: year       ! Gregorian calendar year
INTEGER, INTENT(IN)     :: day        ! Day-number in the year
REAL(RealK), INTENT(IN) :: seconds    ! Seconds since midnight GMT
REAL(RealK), INTENT(IN) :: timestep   ! Length of timestep in seconds

REAL(RealK), INTENT(OUT) :: sindec
!       Sin(solar declination)
REAL(RealK), INTENT(OUT) :: eqt
!       The equation of time, specified as an hour angle in radians.
REAL(RealK), INTENT(OUT) :: scs
!       Solar constant scaling factor

LOGICAL, INTENT(IN), OPTIONAL :: l_observer
!       Calculate angles towards distant observer
LOGICAL, INTENT(IN), OPTIONAL :: l_calendar_360
!       Use 360 day calendar

REAL(RealK), INTENT(OUT), OPTIONAL :: eqt_obs
!       The equation of time, specified as an hour angle in radians for the
!       position of the observer rather than the sun.
REAL(RealK), INTENT(OUT), OPTIONAL :: sindec_obs
!       Sin(observer declination)
REAL(RealK), INTENT(OUT), OPTIONAL :: phase_obs
!       Phase angle of orbit as viewed by observer

! Mathematical constants:
REAL(RealK), PARAMETER :: twopi = 2.0_RealK * pi

! Parameters of the Earth's orbit:
REAL(RealK) :: e       ! Eccentricity of the orbit
REAL(RealK) :: gamph   ! Supplement of the longitude of the perihelion
REAL(RealK) :: oblq    ! Obliquity of the orbit

! Derived orbital constants:
REAL(RealK) :: e1, e2, e3, e4, y, y2, p

REAL(RealK) :: m
!       Mean anomaly: positional angle of a "mean" Earth rotating around
!       the sun with a constant angular speed equal to 2pi/T and counted
!       counterclockwise from the perihelion
REAL(RealK) :: v
!       True anomaly: positional angle of Earth in its orbit, counted
!       counterclockwise from perihelion
REAL(RealK) :: ha
!       Hour angle of the planet (in radians) at Earth midnight
REAL(RealK) :: day_number_at_midnight
!       Number of days from epoch to Earth midnight at the beginning of
!       the current Earth day.
REAL(RealK) :: day_number
!       Number of days from epoch to the middle of the current timestep

REAL(RealK), PARAMETER :: jd1 = 1721425.5_RealK
!       Julian day number for the beginning of the first day (year 1) of
!       the Gregorian calendar.

! Cartesian coordinates of distant observer in orbital and equatorial
! reference frames:
REAL(RealK) :: obs_orb_x, obs_orb_y, obs_orb_z
REAL(RealK) :: obs_eqt_x, obs_eqt_y, obs_eqt_z

LOGICAL :: l_sec_var
!       Use secular variations of the Earth's orbit

LOGICAL :: l_cal360
LOGICAL :: l_obs
!       Local logicals for optional inputs

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='SOLPOS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(l_observer)) THEN
  l_obs = l_observer
ELSE
  l_obs = .FALSE.
END IF

IF (PRESENT(l_calendar_360)) THEN
  l_cal360 = l_calendar_360
ELSE
  l_cal360 = .FALSE.
END IF

SELECT CASE (orbit%i_elements)
CASE (ip_elements_user)

  ! Calculate number of days from epoch to the middle of the 
  ! current timestep. The conversion to Julian day number assumes
  ! the Gregorian calendar for the calculation of leap years.
  day_number_at_midnight = jd1 + REAL( 365*(year-1)     &
    + (year-1)/4 - (year-1)/100 + (year-1)/400 , RealK) &
    + REAL(day, RealK) - 1.0_RealK - orbit%epoch
  day_number = day_number_at_midnight                   &
    + (seconds+timestep/2.0_RealK)/86400.0_RealK

  m    = orbit%mean_anomaly + orbit%mean_anomaly_inc * day_number
  e    = orbit%eccentricity + orbit%eccentricity_inc * day_number
  oblq = orbit%obliquity    + orbit%obliquity_inc    * day_number
  gamph = pi - (orbit%arg_periapsis + orbit%arg_periapsis_inc * day_number)

  SELECT CASE (orbit%i_spin)
  CASE (ip_spin_user)
    ha = orbit%hour_angle + orbit%hour_angle_inc * day_number_at_midnight
  CASE (ip_spin_earth_day)
    ha = orbit%hour_angle + twopi * day_number_at_midnight
  END SELECT

CASE (ip_elements_earth_fixed)

  l_sec_var = .FALSE.
  ! Fixed orbital parameters for the Earth
  CALL orbprm(day, year, l_sec_var, l_cal360, e, gamph, oblq, m)

CASE (ip_elements_earth_secular_variation)

  l_sec_var = .TRUE.
  ! Varying orbital parameters for the Earth
  CALL orbprm(day, year, l_sec_var, l_cal360, e, gamph, oblq, m)

END SELECT

! Calculate the coefficients in the equation of the centre and
! thence derive the true anomaly.
e1 = e * ( 2.0_RealK - 0.25_RealK * e*e )
e2 = 1.25_RealK * e*e
e3 = e*e*e * 13.0_RealK / 12.0_RealK

! True anomaly, equation 87 in Smart on p. 120
v  = m + e1*SIN(m) + e2*SIN(2.0_RealK*m) + e3*SIN(3.0_RealK*m)

! Solar constant scaling factor
e4  = ( 1.0_RealK / (1.0_RealK - e*e) )**2
scs = e4 * ( 1.0_RealK + e * COS(v) )**2
IF (orbit%i_elements == ip_elements_user) THEN
  scs = scs / (orbit%semimajor_axis + orbit%semimajor_axis_inc*day_number)**2
END IF

! sin(solar declination)
! The solar declination is related to
!  the true longitude of the earth (lambda) by:
!  sindec = sin(obliquity) * sin(lambda)
! Lambda is counted counterclockwise from the vernal equinox
!  and is related to v (the true anomaly) through
!  lambda = v + (longitude of perihelion)
sindec = SIN(oblq) * SIN (v - gamph)


! Calculate the equation of time.
SELECT CASE (orbit%i_eqt)
CASE (ip_smart)
  ! Use equation (29) on page 149 of Smart (1944).
  ! (Recall the factor of 5/4 in the definition of e2).
  y   = ( TAN ( 0.5_RealK*oblq ) )**2
  eqt = y * SIN(2.0_RealK * ( m - gamph ))                &
        - 2.0_RealK*e * SIN(m)                            &
        + 4.0_RealK*e*y * SIN(m) * COS(2.0_RealK*( m - gamph )) &
        - 0.5_RealK*y*y * SIN(4.0_RealK*( m - gamph ))          &
        - e2 * SIN(2.0_RealK*m)
CASE (ip_mueller)
  ! M. Mueller, Acta Physica Polonica A 88 Supplement, S-49 (1995)
  y   = ( TAN ( 0.5_RealK*oblq ) )**2
  y2  = y*y
  e2  = e*e 
  p   = 0.5_RealK*pi - gamph
  eqt = - y * (1.0_RealK-4.0_RealK*e2) * SIN(2.0_RealK*(m+p)) &
        - 2.0_RealK  * e * SIN(m)                 &
        + 2.0_RealK  * e * y  * SIN(m+2.0_RealK*p)      &
        - 2.0_RealK  * e * y  * SIN(3.0_RealK*m+2.0_RealK*p)  &
        - 0.5_RealK  * y2 * SIN(4.0_RealK*(m+p))        &
        - 1.25_RealK * e2 * SIN(2.0_RealK*m)            &
        + 2.0_RealK  * e * y2 * SIN(3.0_RealK*m+4.0_RealK*p)  &
        - 2.0_RealK  * e * y2 * SIN(5.0_RealK*m+4.0_RealK*p)  &
        - 3.25_RealK * e2* y  * SIN(4.0_RealK*m+2.0_RealK*p)  &
        - y2   * y * SIN(6.0_RealK*(m+p))/3.0_RealK
CASE DEFAULT
  eqt=0.0e+00_RealK
END SELECT

IF (orbit%i_elements == ip_elements_user) THEN
  ! Add on a correction for the planets hour angle at Earth midnight
  eqt = eqt + MODULO(ha + pi, twopi)
END IF


IF (l_obs) THEN
  ! Calculations needed to output emission spectra towards a distant observer
  IF (PRESENT(phase_obs)) THEN
    phase_obs = MODULO(v - gamph + pi - orbit%observer_lon, twopi)
  END IF

  ! The observer declination is calculated using coordinate transformation
  ! from the planet's orbital to the planet's equatorial reference frame:
  obs_orb_x = COS(orbit%observer_lat)*COS(orbit%observer_lon)
  obs_orb_y = COS(orbit%observer_lat)*SIN(orbit%observer_lon)
  obs_orb_z = SIN(orbit%observer_lat)
  obs_eqt_x = obs_orb_x
  obs_eqt_y = COS(oblq)*obs_orb_y - SIN(oblq)*obs_orb_z
  obs_eqt_z = SIN(oblq)*obs_orb_y + COS(oblq)*obs_orb_z
  IF (PRESENT(sindec_obs)) sindec_obs = obs_eqt_z
  
  IF (PRESENT(eqt_obs)) THEN
    ! We add the position of the observer on to the equation of time
    ! noting that the equation of time is in a clockwise sense.
    ! First add the direction of the vernal equinox (the zero point for
    ! the longitude of the observer in the planet's equatorial frame):
    eqt_obs = eqt + ATAN2( COS(oblq) * SIN(v - gamph), COS(v - gamph) )
    ! Then subtract the longitude of the observer in the planet's equatorial
    ! reference frame (reduces to orbit%observer_lon if oblq=0):
    eqt_obs = eqt_obs - ATAN2(obs_eqt_y, obs_eqt_x)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE solpos
END MODULE solpos_mod

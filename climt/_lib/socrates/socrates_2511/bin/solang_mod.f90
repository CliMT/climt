! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE solang_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLANG_MOD'
CONTAINS

! Calculation of solar zenith/azimuth angles and sunlit fraction.
!
! Purpose :
!  Calculations of the earth's orbit described in the "Calculation
!  of incoming insolation" section of the technical guide, i.e.
!  from the sin of the solar declination, the position of each point
!  and the time limits it calculates how much sunlight, if any, it
!  receives.
!
!------------------------------------------------------------------------------
SUBROUTINE solang(orbit, t, dt, sindec, eqt, lat, lon, k,                      &
     lit, cosz, sol_azimuth, cosz_beg, cosz_end)

USE realtype_rd, ONLY: RealK
USE rad_ccf, ONLY: pi, seconds_per_day
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE def_orbit, ONLY: StrOrbit, &
  ip_spin_user, ip_spin_earth_day, ip_spin_fixed_sun

IMPLICIT NONE


TYPE(StrOrbit), INTENT(IN) :: orbit ! Orbital elements

INTEGER, INTENT(IN) :: k
!   Number of points

REAL(RealK), INTENT(IN) ::                                                     &
  sindec,                                                                      &
!   Sin(solar declination)
  t, dt,                                                                       &
!   Start time (GMT) & timestep
  lat(k), lon(k),                                                              &
!   latitude & longitude of each point
  eqt
!   The equation of time (the difference between true and mean solar time)

REAL(RealK), INTENT(OUT) ::                                                    &
  lit(k),                                                                      &
!   Sunlit fraction of the timestep
  cosz(k)
!   Mean cos(solar zenith angle) during the sunlit fraction

REAL(RealK), INTENT(OUT), OPTIONAL ::                                          &
  sol_azimuth(k),                                                              &
!   Mean solar azimuth angle (radians clockwise from grid north) during the
!   sunlit fraction
  cosz_beg(k), cosz_end(k)
!   Cosine of the solar zenith angle at the beginning and end of the
!   period over which cosz is integrated


INTEGER :: j
!   Loop counter over points

REAL(RealK) ::                                                                 &
  mean_omega(k),                                                               &
!   Mean hour angle over the timestep measured in radians west of local noon
  delta_hour_angle,                                                            &
!   Increment to hour angle per Earth day
  seconds_to_radians,                                                          &
!   Convert seconds to the angle through which the planet has spun
  sinsin, coscos,                                                              &
!   Products of the sines and cosines of solar declination and latitude
  hld, coshld,                                                                 &
!   Half-length of the day in radians (equal to the hour-angle of sunset,
!   and minus the hour-angle of sunrise) & its cosine.
  hat,                                                                         &
!   Local hour angle at the start time.
  omegab, omegae, omega1, omega2, omegas,                                      &
!   Local hour angle at the beginning and end of the timestep and
!   of the period over which cosz is integrated, and sunset - all measured in
!   radians after local sunrise, not from local noon as the true hour angle is.
  difsin, diftim,                                                              &
!   A difference-of-sines intermediate value and the corresponding time period
  trad, dtrad,                                                                 &
!   These are the start-time and length of the timestep (T & DT)
!   converted to radians after midday GMT, or equivalently, hour
!   angle of the mean sun on the Greenwich meridian.
  sinlat(k), dec
!   Working variables for bearing calculation.

REAL(RealK), PARAMETER :: twopi = 2.0_RealK*pi
REAL(RealK), PARAMETER :: eps = EPSILON(1.0_RealK)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLANG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (orbit%i_spin)
CASE (ip_spin_user)
  delta_hour_angle = orbit%hour_angle_inc
CASE (ip_spin_earth_day)
  delta_hour_angle = twopi
END SELECT

IF (orbit%i_spin == ip_spin_fixed_sun) THEN

  ! Fix the sun at a particular zenith and azimuth angle for all points
  ! for idealised tests.
  cosz = COS(orbit%fixed_zenith_angle)
  WHERE (cosz < eps)
    cosz = 0.0_RealK
    lit = 0.0_RealK
  ELSEWHERE
    lit = 1.0_RealK
  END WHERE
  IF (PRESENT(cosz_beg)) cosz_beg = cosz
  IF (PRESENT(cosz_end)) cosz_end = cosz
  IF (PRESENT(sol_azimuth)) sol_azimuth = orbit%fixed_azimuth_angle

ELSE IF (ABS(delta_hour_angle) < eps) THEN

  ! Planet is tidally locked.
  ! The longitude of the overhead sun may vary with the equation of time
  ! to allow for an eccentric orbit and a constant rotation rate.
  mean_omega = MODULO(lon + eqt - pi, twopi)

  ! The solar declination is allowed to vary with the orbital parameters.
  dec = ASIN(sindec)
  cosz = sindec * SIN(lat) + COS(dec)*COS(lat)*COS(mean_omega)
  WHERE (cosz < eps)
    cosz = 0.0_RealK
    lit = 0.0_RealK
  ELSEWHERE
    lit = 1.0_RealK
  END WHERE
  IF (PRESENT(cosz_beg)) cosz_beg = cosz
  IF (PRESENT(cosz_end)) cosz_end = cosz
  IF (PRESENT(sol_azimuth)) THEN
    sol_azimuth = MODULO(ATAN2( -COS(dec)*SIN(mean_omega),                     &
      COS(lat)*sindec - SIN(lat)*COS(dec)*COS(mean_omega) ), twopi)
  END IF

ELSE

  seconds_to_radians = delta_hour_angle/seconds_per_day
  trad = t * seconds_to_radians - pi
  dtrad = dt * seconds_to_radians
  !DIR$ IVDEP
  DO j = 1, k
    coshld = 0.0_RealK ! Logically unnecessary statements without which
    hld = 0.0_RealK    ! the CRAY compiler will not vectorize this code.
    sinlat(j) = SIN(lat(j))
    sinsin = sindec * sinlat(j)
    coscos = SQRT( (1.0_RealK-sindec**2) * (1.0_RealK-sinlat(j)**2) )
    IF ( sinsin < -coscos ) THEN
      ! Perpetual night
      lit(j) = 0.0_RealK
      cosz(j) = 0.0_RealK
      IF (PRESENT(cosz_beg)) cosz_beg(j) = 0.0_RealK
      IF (PRESENT(cosz_end)) cosz_end(j) = 0.0_RealK
      mean_omega(j) = 0.0_RealK
    ELSE
      ! Since T and DT represent mean time and all solar calculations
      ! are done using the variable HAT, it suffices to add the
      ! equation of time on to HAT.
      hat = lon(j) + trad + eqt
      IF ( sinsin > coscos ) THEN
        ! Perpetual day: hour angles are start and end of timestep
        omega1 = hat
        omega2 = hat + dtrad
      ELSE
        ! At this latitude some points are sunlit, some not. Different ones
        ! need different treatment.
        coshld = sinsin / coscos
        hld = ACOS(-coshld)
        ! The logic seems simplest if one takes all "times" - actually hour
        ! angles - relative to sunrise (or sunset), but they must be kept in the
        ! range 0 to 2pi for the tests on their orders to work.
        omegab = MODULO(hat + hld, twopi)
        omegae = MODULO(omegab + dtrad, twopi)
        omegas = 2.0_RealK * hld
        ! Now that the start-time, end-time and sunset are set in terms of
        ! hour angle, we can set the limiting hour-angles for integration.
        ! The simple cases are start-to-end-of-timestep, start-to-sunset,
        ! sunrise-to-end and sunrise-to-sunset, but two other cases exist
        ! and need special treatment.
        IF (omegab <= omegas .OR. omegab < omegae) THEN
          omega1 = omegab - hld
        ELSE
          omega1 = - hld
        END IF
        IF (omegae <= omegas) THEN
          omega2 = omegae - hld
        ELSE
          omega2 = omegas - hld
        END IF
        IF (omegae > omegab .AND. omegab > omegas) omega2=omega1
        ! Put in an arbitrary marker for the case when the sun does not rise
        ! during the timestep (though it is up elsewhere at this latitude).
      END IF ! This finishes the ELSE (perpetual day) block
      difsin = SIN(omega2) - SIN(omega1)
      diftim = omega2 - omega1
      mean_omega(j) = (omega1 + omega2)/2.0_RealK
      ! Next, deal with the case where the sun sets and then rises again
      ! within the timestep. There the integration has actually been done
      ! backwards over the night, and the resulting negative DIFSIN and DIFTIM
      ! must be combined with positive values representing the whole of the
      ! timestep to get the right answer, which combines contributions from
      ! the two separate daylit periods. A simple analytic expression for the
      ! total sun throughout the day is used. (This could of course be used
      ! alone at points where the sun rises and then sets within the timestep)
      IF (diftim < 0.0_RealK) THEN
        difsin = difsin + 2.0_RealK * SQRT(1.0_RealK-coshld**2)
        diftim = diftim + 2.0_RealK * hld
        mean_omega(j) = mean_omega(j) + Pi
      END IF
      IF (mean_omega(j) > pi) mean_omega(j)=mean_omega(j)-twopi
      IF (diftim == 0.0_RealK) THEN
        ! Pick up the arbitrary marker for night points at a partly-lit latitude
        cosz(j) = 0.0_RealK
        IF (PRESENT(cosz_beg)) cosz_beg(j) = 0.0_RealK
        IF (PRESENT(cosz_end)) cosz_end(j) = 0.0_RealK
        lit(j) = 0.0_RealK
      ELSE
        cosz(j) = MIN(MAX(difsin*coscos/diftim + sinsin, 0.0_RealK), 1.0_RealK)
        IF (PRESENT(cosz_beg)) cosz_beg(j) = coscos*COS(omega1) + sinsin
        IF (PRESENT(cosz_end)) cosz_end(j) = coscos*COS(omega2) + sinsin
        lit(j) = diftim / dtrad
      END IF
    END IF ! This finishes the ELSE (perpetual night) block
  END DO

  dec = ASIN(sindec)
  IF (PRESENT(sol_azimuth)) THEN
    sol_azimuth = MODULO(ATAN2( -COS(dec)*SIN(mean_omega),                     &
      COS(lat)*sindec - sinlat*COS(dec)*COS(mean_omega) ), twopi)
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE solang
END MODULE solang_mod

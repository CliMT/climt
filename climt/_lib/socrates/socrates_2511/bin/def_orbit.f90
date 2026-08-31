! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for orbital elements.
!
! Description:
!   This module contains the declaration of the structure
!   used to store orbital elements in the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_orbit

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrOrbit

  ! Selector for use of orbital elements
  INTEGER :: i_elements

  ! Selector for treatment of stellar motion from planet's surface
  INTEGER :: i_spin

  ! Selector for formulation of equation of time
  INTEGER :: i_eqt

  ! Epoch in Julian Days
  REAL(RealK) :: epoch

  ! Eccentricity
  REAL(RealK) :: eccentricity

  ! Increment to eccentricity per day number from epoch
  REAL(RealK) :: eccentricity_inc

  ! Longitude of perihelion in radians
  REAL(RealK) :: arg_periapsis

  ! Increment to longitude of perihelion per day number from epoch
  REAL(RealK) :: arg_periapsis_inc

  ! Obliquity of the orbit in radians
  REAL(RealK) :: obliquity

  ! Increment to obliquity of the orbit per day number from epoch
  REAL(RealK) :: obliquity_inc

  ! Semi-major axis in AU
  REAL(RealK) :: semimajor_axis

  ! Increment to semi-major axis per day number from epoch
  REAL(RealK) :: semimajor_axis_inc

  ! Mean anomaly at epoch in radians
  REAL(RealK) :: mean_anomaly

  ! Increment to mean anomaly per day number from epoch
  REAL(RealK) :: mean_anomaly_inc

  ! Hour angle at epoch in radians
  REAL(RealK) :: hour_angle

  ! Increment to hour angle per day number from epoch
  REAL(RealK) :: hour_angle_inc

  ! Fixed stellar zenith angle for all points in radians
  REAL(RealK) :: fixed_zenith_angle

  ! Fixed stellar azimuth angle measured clockwise from grid north in radians
  REAL(RealK) :: fixed_azimuth_angle

  ! Orbital longitude of observer
  REAL(RealK) :: observer_lon

  ! Orbital latitude of observer
  REAL(RealK) :: observer_lat

END TYPE StrOrbit

! Identifiers for use of orbital elements
INTEGER, PARAMETER :: ip_elements_user = 1
!   User defined orbit
INTEGER, PARAMETER :: ip_elements_earth_fixed = 2
!   Fixed Earth orbit
INTEGER, PARAMETER :: ip_elements_earth_secular_variation = 3
!   Secular variations of the Earth's orbit

! Identifiers for treatment of stellar motion from planet's surface
INTEGER, PARAMETER :: ip_spin_user = 1
!   User defined hour angle increments
INTEGER, PARAMETER :: ip_spin_earth_day = 2
!   Diurnal cycle of exactly one Earth day
INTEGER, PARAMETER :: ip_spin_fixed_sun = 3
!   Fixed stellar position for all points

! Identifiers for equation of time formulation:
INTEGER, PARAMETER :: ip_smart   = 1
!   Equation 29 on page 149 of Smart (1944).
INTEGER, PARAMETER :: ip_mueller = 2 
!   M. Mueller, Acta Physica Polonica A 88 Supplement, S-49 (1995)


CONTAINS


SUBROUTINE set_orbit(orbit, &
  i_elements, i_spin, i_eqt, &
  epoch, eccentricity, eccentricity_inc, arg_periapsis, arg_periapsis_inc, &
  obliquity, obliquity_inc, semimajor_axis, semimajor_axis_inc, &
  mean_anomaly, mean_anomaly_inc, hour_angle, hour_angle_inc, &
  fixed_zenith_angle, fixed_azimuth_angle, observer_lon, observer_lat)

USE realtype_rd, ONLY: RealExt

IMPLICIT NONE

TYPE (StrOrbit), INTENT(OUT) :: orbit

INTEGER, INTENT(IN), OPTIONAL :: i_elements, i_spin, i_eqt

REAL(RealExt), INTENT(IN), OPTIONAL :: epoch
REAL(RealExt), INTENT(IN), OPTIONAL :: eccentricity, eccentricity_inc
REAL(RealExt), INTENT(IN), OPTIONAL :: arg_periapsis, arg_periapsis_inc
REAL(RealExt), INTENT(IN), OPTIONAL :: obliquity, obliquity_inc
REAL(RealExt), INTENT(IN), OPTIONAL :: semimajor_axis, semimajor_axis_inc
REAL(RealExt), INTENT(IN), OPTIONAL :: mean_anomaly, mean_anomaly_inc
REAL(RealExt), INTENT(IN), OPTIONAL :: hour_angle, hour_angle_inc
REAL(RealExt), INTENT(IN), OPTIONAL :: fixed_zenith_angle, fixed_azimuth_angle
REAL(RealExt), INTENT(IN), OPTIONAL :: observer_lon, observer_lat

! Defaults to use when arguments are not present
INTEGER, PARAMETER :: default_i_elements = ip_elements_user
INTEGER, PARAMETER :: default_i_spin = ip_spin_user
REAL(RealK), PARAMETER :: default_epoch = 2451545.0_RealK ! (J2000)
REAL(RealK), PARAMETER :: default_a = 1.0_RealK
REAL(RealK), PARAMETER :: zero_element = 0.0_RealK

IF (PRESENT(i_elements)) THEN
  orbit%i_elements = i_elements
ELSE
  orbit%i_elements = default_i_elements
END IF

IF (PRESENT(i_spin)) THEN
  orbit%i_spin = i_spin
ELSE
  orbit%i_spin = default_i_spin
END IF

! Set orbital elements
SELECT CASE (orbit%i_elements)
CASE (ip_elements_user)
  IF (PRESENT(i_eqt)) THEN
    orbit%i_eqt = i_eqt
  ELSE
    orbit%i_eqt = ip_mueller
  END IF
  IF (PRESENT(epoch)) THEN
    orbit%epoch = real(epoch, RealK)
  ELSE
    orbit%epoch = default_epoch
  END IF
  IF (PRESENT(eccentricity)) THEN
    orbit%eccentricity = real(eccentricity, RealK)
  ELSE
    orbit%eccentricity = zero_element
  END IF
  IF (PRESENT(eccentricity_inc)) THEN
    orbit%eccentricity_inc = real(eccentricity_inc, RealK)
  ELSE
    orbit%eccentricity_inc = zero_element
  END IF
  IF (PRESENT(arg_periapsis)) THEN
    orbit%arg_periapsis = real(arg_periapsis, RealK)
  ELSE
    orbit%arg_periapsis = zero_element
  END IF
  IF (PRESENT(arg_periapsis_inc)) THEN
    orbit%arg_periapsis_inc = real(arg_periapsis_inc, RealK)
  ELSE
    orbit%arg_periapsis_inc = zero_element
  END IF
  IF (PRESENT(obliquity)) THEN
    orbit%obliquity = real(obliquity, RealK)
  ELSE
    orbit%obliquity = zero_element
  END IF
  IF (PRESENT(obliquity_inc)) THEN
    orbit%obliquity_inc = real(obliquity_inc, RealK)
  ELSE
    orbit%obliquity_inc = zero_element
  END IF
  IF (PRESENT(semimajor_axis)) THEN
    orbit%semimajor_axis = real(semimajor_axis, RealK)
  ELSE
    orbit%semimajor_axis = default_a
  END IF
  IF (PRESENT(semimajor_axis_inc)) THEN
    orbit%semimajor_axis_inc = real(semimajor_axis_inc, RealK)
  ELSE
    orbit%semimajor_axis_inc = zero_element
  END IF
  IF (PRESENT(mean_anomaly)) THEN
    orbit%mean_anomaly = real(mean_anomaly, RealK)
  ELSE
    orbit%mean_anomaly = zero_element
  END IF
  IF (PRESENT(mean_anomaly_inc)) THEN
    orbit%mean_anomaly_inc = real(mean_anomaly_inc, RealK)
  ELSE
    orbit%mean_anomaly_inc = zero_element
  END IF
  IF (PRESENT(hour_angle)) THEN
    orbit%hour_angle = real(hour_angle, RealK)
  ELSE
    orbit%hour_angle = zero_element
  END IF
CASE (ip_elements_earth_fixed)
  IF (PRESENT(i_eqt)) THEN
    orbit%i_eqt = i_eqt
  ELSE
    orbit%i_eqt = ip_smart
  END IF
CASE (ip_elements_earth_secular_variation)
  IF (PRESENT(i_eqt)) THEN
    orbit%i_eqt = i_eqt
  ELSE
    orbit%i_eqt = ip_smart
  END IF
END SELECT

! Set motion of sun across the sky
SELECT CASE (orbit%i_spin)
CASE (ip_spin_user)
  IF (PRESENT(hour_angle_inc)) THEN
    orbit%hour_angle_inc = real(hour_angle_inc, RealK)
  ELSE
    orbit%hour_angle_inc = zero_element
  END IF
CASE (ip_spin_earth_day)
CASE (ip_spin_fixed_sun)
  IF (PRESENT(fixed_zenith_angle)) THEN
    orbit%fixed_zenith_angle = real(fixed_zenith_angle, RealK)
  ELSE
    orbit%fixed_zenith_angle = zero_element
  END IF
  IF (PRESENT(fixed_azimuth_angle)) THEN
    orbit%fixed_azimuth_angle = real(fixed_azimuth_angle, RealK)
  ELSE
    orbit%fixed_azimuth_angle = zero_element
  END IF
END SELECT

IF (PRESENT(observer_lon)) THEN
  orbit%observer_lon = real(observer_lon, RealK)
ELSE
  orbit%observer_lon = zero_element
END IF

IF (PRESENT(observer_lat)) THEN
  orbit%observer_lat = real(observer_lat, RealK)
ELSE
  orbit%observer_lat = zero_element
END IF

END SUBROUTINE set_orbit

END MODULE def_orbit

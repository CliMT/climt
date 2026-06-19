! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to declare a structure for solar spectral data.

MODULE def_solarspec

  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: solar_t_effective, solar_radius

  IMPLICIT NONE


  TYPE StrSolarSpec

    INTEGER :: n_points = 0
!     Number of points in the spectrum
    REAL (RealK), POINTER :: wavelength(:)
!     Wavelength at which the spectral irradiance is specified
    REAL (RealK), POINTER :: irrad(:)
!     Solar spectral irradiance in units of Wm-2.m-1
    REAL (RealK), POINTER :: bandsize(:)
!     Band size
    REAL (RealK), POINTER :: bandbnds(:,:)
!     Bounds of each band
    REAL (RealK) :: t_effective = solar_t_effective
!     Effective solar temperature
    REAL (RealK) :: radius = solar_radius
!     Radius at the photosphere
    LOGICAL :: l_binned = .FALSE.
!     Spectrum is either binned (.TRUE.) or point-values (.FALSE.)

  END TYPE StrSolarSpec

END MODULE def_solarspec

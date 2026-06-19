! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the solar irradiance at given wavenumbers.
!
SUBROUTINE solar_intensity_90(nu, SolarSpec, solar_intensity)
!
!
! Description:
!   This subroutine returns the solar irradiance per unit wavenumber
!   at given wavenumbers.
!
! Method:
!   The wavelength is bracketed between two points of the solar
!   spectrum and linear interpolation is used in the interval.
!   At long wavelngths a Rayleigh-Jeans tail is used.


  USE realtype_rd
  USE def_solarspec
  USE dimensions_pp_ucf
  USE rad_ccf, ONLY: astronomical_unit

  IMPLICIT NONE

  
! Dummy arguments
  REAL  (RealK), Intent(IN) :: nu(:)
!   Wavenumbers supplied to the routine
  TYPE (StrSolarSpec), Intent(IN) :: SolarSpec
!   Spectral solar irradiance at the top of the atmosphere
  REAL  (RealK), Intent(OUT) :: solar_intensity(:)
!   Returned intensity

! Local variables.
  INTEGER :: n
!   Size of supplied array of wavenumbers
  INTEGER :: i
!   Loop variable
  INTEGER :: i_short
!   Array point just below lambda
  INTEGER :: i_long
!   Array point just above lambda
  REAL  (RealK) :: lambda
!   Wavelength
  REAL  (RealK) :: s(1)
!   Array used for conformance in the function call
  REAL  (RealK) :: fraction
!   Fraction of interval covered


! Subroutines called:
  INTERFACE
!
    SUBROUTINE planck_90(nu, t, b)
!     Subroutine to calculate the Planckian radiance
!
      USE realtype_rd
!
      REAL  (RealK), Intent(IN) :: nu(:)
      REAL  (RealK), Intent(IN) :: t
!
      REAL  (RealK), Intent(OUT) :: b(:)
!
    END SUBROUTINE planck_90
!
  END INTERFACE


! Solar spectra are defined by wavelength.
  n=SIZE(nu)

  DO i=1, n
!
    lambda=1.0_RealK/nu(i)
    CALL point_bracket(lambda, SolarSpec%n_points, SolarSpec%wavelength, &
      i_short, i_long)
!
!   If lambda lies within the range covered we interpolate: if not
!   use a simple extrapolation
    IF (i_short == 0) THEN
!
!     Lambda is at a wavelength shorter than any in the spectrum.
      solar_intensity(i)=0.0_RealK
!
    ELSE IF (i_long >  SolarSpec%n_points) THEN
!
!     Use a black body fit at the effective temperature.
      CALL planck_90( (/ nu(i) /), SolarSpec%t_effective, s)
      solar_intensity(i) = s(1) * (SolarSpec%radius / astronomical_unit)**2
!
    ELSE
!
!     Within the spectrum use linear interpolation.
      fraction = (lambda - SolarSpec%wavelength(i_short)) / &
                 (SolarSpec%wavelength(i_long) - &
                  SolarSpec%wavelength(i_short))
      s = (fraction * SolarSpec%irrad(i_long) + &
        (1.0_RealK - fraction) * SolarSpec%irrad(i_short) ) * &
        lambda**2
      solar_intensity(i) = s(1)
!
    ENDIF
!
  ENDDO

END SUBROUTINE solar_intensity_90

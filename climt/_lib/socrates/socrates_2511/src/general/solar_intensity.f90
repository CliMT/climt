! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the solar irradiance at a given wavelength.
!
! Method:
!     The wavelength is bracketed between two points of the solar
!     spectrum and linear interpolation is used in the interval.
!     At long wavelngths a Rayleigh-Jeans tail is used.
!
!- ---------------------------------------------------------------------
FUNCTION solar_intensity(lambda, Sol)

  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: astronomical_unit
  USE def_solarspec, ONLY: StrSolarSpec

  IMPLICIT NONE


  REAL (RealK) :: solar_intensity
!   Returned intensity

  REAL (RealK), Intent(IN) :: lambda
!   Wavelength
  TYPE (StrSolarSpec), INTENT(IN) :: Sol
!   Solar spectrum

! Local variables.
  INTEGER :: i_short
!   Array point just below lambda
  INTEGER :: i_long
!   Array point just above lambda
  REAL  (RealK) :: fraction
!   Fraction of interval covered
!

! Functions called:
  REAL (RealK), EXTERNAL :: planck
!   Planckian function


! Find the points of the wavelength array bracketting lambda.
  CALL point_bracket(lambda, Sol%n_points, Sol%wavelength, i_short, i_long)

! If lambda lies within the range covered we interpolate: if not
! use a simple extrapolation
  IF (i_short == 0) THEN
!   Lambda is at a wavelength shorter than any in the spectrum.
    solar_intensity=0.0_RealK
  ELSE IF (i_long >  Sol%n_points) THEN
!   Use a black body fit at the effective temperature.
    solar_intensity=planck(Sol%t_effective, lambda) &
      *(Sol%radius/astronomical_unit)**2
  ELSE
!   Within the spectrum use linear interpolation.
    fraction=(lambda-Sol%wavelength(i_short)) &
      /(Sol%wavelength(i_long)-Sol%wavelength(i_short))
    solar_intensity=fraction*Sol%irrad(i_long) &
      +(1.0_RealK-fraction)*Sol%irrad(i_short)
  ENDIF

END FUNCTION solar_intensity

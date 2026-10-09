! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set the possible types of line profile corrections.

MODULE line_prof_corr_mod

! Description:
!   Corrections of Voigt line profiles are usually parameterised in
!   terms of a correction factor. This module defines the different
!   types of correction factors supported.

USE realtype_rd, ONLY: RealK

IMPLICIT NONE

! Definitions of line profile correction types
  INTEGER, PARAMETER :: ip_lpc_unity      = 1
!   No correction applied
  INTEGER, PARAMETER :: ip_lpc_co2_ph1989 = 2
!   Corrections to CO2 line profiles from Perrin & Hartmann (1989)

! Various constants used in line profile correction calculations
  REAL (RealK) :: &
    b_co2(3), b_n2(3)

CONTAINS

! Main function to calculate line profile correction
REAL (RealK) FUNCTION line_prof_corr(nu, nu_line, alpha_lorentz_air &
  , alpha_lorentz_self, i_line_prof_corr)

! Input variables
  REAL (RealK), INTENT(IN) :: &
      nu &
!       Current wavenumber
    , nu_line &
!       Wavenumber of line center
    , alpha_lorentz_air &
!       Lorentz HWHM for broadening by air
    , alpha_lorentz_self
!       Lorentz HWHM for broadening by 
  INTEGER, INTENT(IN) :: &
      i_line_prof_corr
!       Line correction type

! Local variables
  REAL (RealK) :: &
      delta_nu &
!       Distance from line center in cm-1
    , line_prof_corr_air &
!       Correction of line profile due to broadening by air
    , line_prof_corr_self
!       Correction of line profile due to self-broadening

  SELECT CASE (i_line_prof_corr)

  CASE (ip_lpc_unity)
    line_prof_corr = 1.0_RealK

  CASE (ip_lpc_co2_ph1989)
    delta_nu = ABS(nu - nu_line)*0.01_RealK

!   Correction of N2/air-broadened line profile
    IF (delta_nu < 3.0_RealK) THEN
      line_prof_corr_air = 1.0_RealK
    ELSE IF (delta_nu < 10.0_RealK) THEN
      line_prof_corr_air = EXP(-b_n2(1)*(delta_nu - 3.0_RealK))
    ELSE IF (delta_nu < 70.0_RealK) THEN
      line_prof_corr_air = EXP(-b_n2(1)*(10.0_RealK - 3.0_RealK) &
                               -b_n2(2)*(delta_nu - 10.0_RealK))
    ELSE
      line_prof_corr_air = EXP(-b_n2(1)*(10.0_RealK - 3.0_RealK) &
                               -b_n2(2)*(70.0_RealK - 10.0_RealK) &
                               -b_n2(3)*(delta_nu - 70.0_RealK))
    END IF

!   Correction of CO2-broadened line profile
    IF (delta_nu < 3.0_RealK) THEN
      line_prof_corr_self = 1.0_RealK
    ELSE IF (delta_nu < 30.0_RealK) THEN
      line_prof_corr_self = EXP(-b_co2(1)*(delta_nu - 3.0_RealK))
    ELSE IF (delta_nu < 120.0_RealK) THEN
      line_prof_corr_self = EXP(-b_co2(1)*(30.0_RealK - 3.0_RealK) &
                                -b_co2(2)*(delta_nu - 30.0_RealK))
    ELSE
      line_prof_corr_self = EXP(-b_co2(1)*(30.0_RealK - 3.0_RealK) &
                                -b_co2(2)*(120.0_RealK - 30.0_RealK) &
                                -b_co2(3)*(delta_nu - 120.0_RealK))
    END IF

!   Total line correction is a weighted mean
    line_prof_corr = (alpha_lorentz_air*line_prof_corr_air &
                      +alpha_lorentz_self*line_prof_corr_self)/ &
                     (alpha_lorentz_air + alpha_lorentz_self)

  CASE DEFAULT
    WRITE(*,'(a)') &
      '***Error: Invalid line profile type.'
    STOP

  END SELECT

END FUNCTION line_prof_corr

! Set constants needed in calculation of line profile corrections
SUBROUTINE set_line_prof_corr_cnst(t_calc, i_line_prof_corr)

  REAL (RealK), INTENT(IN) :: t_calc
!   Temperature
  INTEGER, INTENT(IN) :: i_line_prof_corr
!   Line profile correction type

  SELECT CASE (i_line_prof_corr)

  CASE (ip_lpc_co2_ph1989)
    b_co2(1) = ph1989_b(0.0888_RealK, -0.16_RealK,   0.0041_RealK,  t_calc)
    b_co2(2) = ph1989_b(0.0_RealK,     0.0526_RealK, 0.00152_RealK, t_calc)
    b_co2(3) = ph1989_b(0.0232_RealK,  0.0_RealK,    0.0_RealK,     t_calc)
    b_n2(1)  = ph1989_b(0.416_RealK,  -0.354_RealK,  0.00386_RealK, t_calc)
    b_n2(2)  = ph1989_b(0.00167_RealK, 0.0421_RealK, 0.00248_RealK, t_calc)
    b_n2(3)  = ph1989_b(0.02_RealK,    0.0_RealK,    0.0_RealK,     t_calc)
  
  END SELECT

END SUBROUTINE set_line_prof_corr_cnst

! Function to calculate temperature dependent B1, B2 and B3 factors in
! Perrin & Hartmann (1989)
REAL (RealK) FUNCTION ph1989_b(alpha, beta, epsilon, t_calc)

! Input variables
  REAL (RealK), INTENT(IN) :: alpha, beta, epsilon
!   Constants in parametrisation
  REAL (RealK), INTENT(IN) :: t_calc
!   Temperature

! Local variables
  REAL (RealK) :: t_calc_loc
!   Temperature restricted to be within validity range of parametrisation

  t_calc_loc = MIN(800.0_RealK, MAX(190.0_RealK, t_calc))
  ph1989_b = alpha + beta*EXP(-epsilon*t_calc_loc)

END FUNCTION ph1989_b

END MODULE line_prof_corr_mod

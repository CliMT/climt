! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Evaluate Planck flux for a given spectral band
!
!------------------------------------------------------------------------------
MODULE planck_flux_band_mod

USE realtype_rd, ONLY: RealK

IMPLICIT NONE

CONTAINS

REAL(RealK) FUNCTION planck_flux_band(spectrum, i_band, temperature)

  USE def_spectrum, ONLY: StrSpecData
  USE interp1d_mod, ONLY: interp1d

  IMPLICIT NONE

  TYPE (StrSpecData), INTENT(IN) :: spectrum
  INTEGER, INTENT(IN) :: i_band
  REAL(RealK), INTENT(IN) :: temperature

  INTEGER :: j
  REAL(RealK) :: t_ratio

  
  t_ratio=temperature/spectrum%planck%t_ref_planck

  IF (spectrum%planck%l_planck_tbl) THEN
    ! Use a look-up table to the Planck function
    planck_flux_band=interp1d(spectrum%planck%theta_planck_tbl, &
      spectrum%planck%thermal_coeff(:,i_band), &
      t_ratio, spectrum%dim%nd_thermal_coeff)
  ELSE
    ! Use the polynomial fit to the Planck function
    planck_flux_band &
      = spectrum%planck%thermal_coeff(spectrum%planck%n_deg_fit, i_band)
    DO j=spectrum%planck%n_deg_fit-1, 0, -1
      planck_flux_band = planck_flux_band * t_ratio &
        + spectrum%planck%thermal_coeff(j, i_band)
    END DO
  END IF

END FUNCTION planck_flux_band
END MODULE planck_flux_band_mod

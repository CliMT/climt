! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the weightings for correlated-k.
!
SUBROUTINE rad_weight_90 &
!
(i_weight, nu_wgt, SolarSpec, t, wgt)
!
!
! Description:
!   Weightings at the frequencies supplied are calculated. nu_wgt and wgt
!   should have the same lengths.
!
!
!
! Modules used:
  USE realtype_rd
  USE def_solarspec
  USE weighting_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  INTEGER, Intent(IN) :: i_weight
!   Method of weighting
  TYPE (StrSolarSpec), Intent(IN) :: SolarSpec
!   Solar spectral irradiance data
  REAL  (RealK), Intent(IN) :: t
!   Temperatures for Planckian weighting
  REAL  (RealK), Intent(IN) :: nu_wgt(:)
!   Frequencies where weighting is applied
  REAL  (RealK), Intent(OUT) :: wgt(:)
!   Calculated weightings
!
! Local variables
!
!
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
!
    SUBROUTINE d_planck_90(nu, t, db)
!     Subroutine to calculate the differential Planckian radiance
!
      USE realtype_rd
!
      REAL  (RealK), Intent(IN) :: nu(:)
      REAL  (RealK), Intent(IN) :: t
!
      REAL  (RealK), Intent(OUT) :: db(:)
!
    END SUBROUTINE d_planck_90
!
!
    SUBROUTINE solar_intensity_90(nu, SolarSpec, solar_irrad)
!     Subroutine to calculate the solar intensity
!
      USE realtype_rd
      USE def_solarspec
!
      REAL  (RealK), Intent(IN)        :: nu(:)
      TYPE  (StrSolarSpec), Intent(IN) :: SolarSpec
!
      REAL  (RealK), Intent(OUT) :: solar_irrad(:)
!
    END SUBROUTINE solar_intensity_90
!
  END INTERFACE
!
!
  SELECT CASE (i_weight)
    CASE (IP_weight_planck)
      CALL planck_90(nu_wgt, t, wgt)
    CASE (IP_weight_d_planck)
      CALL d_planck_90(nu_wgt, t, wgt)
    CASE (IP_weight_solar)
      CALL solar_intensity_90(nu_wgt, SolarSpec, wgt)
    CASE (IP_weight_uniform)
      wgt(:)=1.0_RealK
  END SELECT
!
!
!
END SUBROUTINE rad_weight_90

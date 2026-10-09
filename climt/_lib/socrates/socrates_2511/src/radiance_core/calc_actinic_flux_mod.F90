! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the actinic flux
!
!------------------------------------------------------------------------------
MODULE calc_actinic_flux_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE calc_actinic_flux(control, sph_comp, sph_common, &
  n_profile, n_layer, tau_abs, flux_total, flux_direct, sec_0, &
  l_scale_solar, adjust_solar_ke, &
  actinic_flux, &
  nd_profile, nd_layer)

  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_spherical_geometry, ONLY: StrSphComp, StrSphCommon

  IMPLICIT NONE

! Control options:
  TYPE(StrCtrl),    INTENT(IN) :: control

! Spherical geometry fields for component (either all-sky or clear-sky)
  TYPE(StrSphComp), INTENT(IN) :: sph_comp

! Spherical geometry fields common to all components
  TYPE(StrSphCommon), INTENT(IN) :: sph_common

  INTEGER, INTENT(IN) :: nd_profile
!   Size allocated for atmospheric profiles
  INTEGER, INTENT(IN) :: nd_layer
!   Size allocated for atmospheric layers

  INTEGER, INTENT(IN) :: n_profile
!   Number of profiles
  INTEGER, INTENT(IN) :: n_layer
!   Number of layers

  REAL(RealK), INTENT(IN) :: tau_abs(nd_profile, nd_layer)
!   Optical depth due to absorption
  REAL(RealK), INTENT(IN) :: flux_total(nd_profile, 2*nd_layer+2)
!   Total flux
  REAL(RealK), INTENT(IN) :: flux_direct(nd_profile, 0: nd_layer)
!   Direct flux (plane-parallel geometry)
  REAL(RealK), INTENT(IN) :: sec_0(nd_profile)
!   Secants of the solar zenith angles

  LOGICAL, INTENT(IN) :: l_scale_solar
!   Apply scaling to solar flux
  REAL (RealK), INTENT(IN) :: adjust_solar_ke(nd_profile, nd_layer)
!   Adjustment of solar beam with equivalent extinction

  REAL(RealK), INTENT(OUT) :: actinic_flux(nd_profile, nd_layer)
!   Actinic flux

  INTEGER :: i, l
!   Loop variables
  REAL(RealK), PARAMETER :: tol = SQRT(EPSILON(1.0_RealK))
!   Tolerance for weak absorption
  REAL(RealK), PARAMETER :: tol_ke = 0.1_RealK
!   Tolerance for change in absorption from equivalent extinction
  REAL(RealK), PARAMETER :: tiny_real = TINY(1.0_RealK)
!   Minimum value to prevent LOG(0)
  REAL(RealK) :: adjust_tau
!   Adjustment to tau with equivalent extinction


! Calculate the actinic flux from the flux divergence
! (or the incoming fluxes where the absorption is very weak)
  IF (control%l_spherical_solar) THEN
    IF (l_scale_solar) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          adjust_tau = LOG(MAX(adjust_solar_ke(l, i), tiny_real)) &
                     / sph_common%path_div(l, i)
          IF (tau_abs(l, i) - adjust_tau > tol) THEN
            IF (ABS(adjust_tau)/tau_abs(l, i) < tol_ke) THEN
              ! Here the adjustment to the total absorption from equivalent
              ! extinction is within the tolerance, so calculation of the
              ! actinic flux from the flux divergence is valid.
              actinic_flux(l, i) = &
                ( sph_comp%flux_direct_div(l, i) &
                + flux_total(l, 2*i+1) - flux_total(l, 2*i-1) &
                + flux_total(l, 2*i) - flux_total(l, 2*i+2) ) &
                / ( tau_abs(l, i) - adjust_tau )
            ELSE
              ! Where the adjustment due to equivalent extinction is large
              ! we fall back to a calculation from the incoming fluxes.
              actinic_flux(l, i) = &
                2.0_RealK * ( flux_total(l, 2*i+1) + flux_total(l, 2*i) ) &
                + sph_comp%flux_direct(l, i)*( 1.0_RealK-adjust_solar_ke(l, i) &
                * EXP(-tau_abs(l, i)*sph_common%path_div(l, i)) ) &
                / (sph_common%path_div(l, i) * ( tau_abs(l, i) - adjust_tau ))
            END IF
          ELSE
            actinic_flux(l, i) = 2.0_RealK*flux_total(l, 2*i+1) &
                               + 2.0_RealK*flux_total(l, 2*i) &
                               + sph_comp%flux_direct(l, i)
          END IF
        END DO
      END DO
    ELSE
      DO i=1, n_layer
        DO l=1, n_profile
          IF (tau_abs(l, i) > tol) THEN
            actinic_flux(l, i) = &
              ( sph_comp%flux_direct_div(l, i) &
              + flux_total(l, 2*i+1) - flux_total(l, 2*i-1) &
              + flux_total(l, 2*i) - flux_total(l, 2*i+2) ) &
              / tau_abs(l, i)
          ELSE
            actinic_flux(l, i) = 2.0_RealK*flux_total(l, 2*i+1) &
                               + 2.0_RealK*flux_total(l, 2*i) &
                               + sph_comp%flux_direct(l, i)
          END IF
        END DO
      END DO
    END IF
  ELSE
    IF (l_scale_solar) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          adjust_tau = LOG(MAX(adjust_solar_ke(l, i), tiny_real)) / sec_0(l)
          IF (tau_abs(l, i) - adjust_tau > tol) THEN
            IF (ABS(adjust_tau)/tau_abs(l, i) < tol_ke) THEN
              actinic_flux(l, i) = &
                ( flux_total(l, 2*i+1) - flux_total(l, 2*i-1) &
                + flux_total(l, 2*i) - flux_total(l, 2*i+2) ) &
                / ( tau_abs(l, i) - adjust_tau )
            ELSE
              actinic_flux(l, i) = 2.0_RealK * ( flux_total(l, 2*i+1) &
                + flux_total(l, 2*i) - flux_direct(l, i-1) ) &
                + flux_direct(l, i-1) * ( 1.0_RealK - adjust_solar_ke(l, i) &
                * EXP(-tau_abs(l, i)*sec_0(l)) ) &
                / ( tau_abs(l, i) - adjust_tau )
            END IF
          ELSE
            actinic_flux(l, i) = 2.0_RealK*flux_total(l, 2*i+1) &
                               + 2.0_RealK*flux_total(l, 2*i) &
                               + flux_direct(l, i-1)*(sec_0(l) - 2.0_RealK)
          END IF
        END DO
      END DO
    ELSE
      DO i=1, n_layer
        DO l=1, n_profile
          IF (tau_abs(l, i) > tol) THEN
            actinic_flux(l, i) = &
              ( flux_total(l, 2*i+1) - flux_total(l, 2*i-1) &
              + flux_total(l, 2*i) - flux_total(l, 2*i+2) ) &
              / tau_abs(l, i)
          ELSE
            actinic_flux(l, i) = 2.0_RealK*flux_total(l, 2*i+1) &
                               + 2.0_RealK*flux_total(l, 2*i) &
                               + flux_direct(l, i-1)*(sec_0(l) - 2.0_RealK)
          END IF
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE calc_actinic_flux
END MODULE calc_actinic_flux_mod

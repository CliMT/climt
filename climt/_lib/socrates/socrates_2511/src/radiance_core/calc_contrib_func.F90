! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the contribution function.
!
! Method:
!   Calculates the contribution, or weighting, function per decade in 
!   pressure.
!
!------------------------------------------------------------------------------
MODULE calc_contrib_func_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_CONTRIB_FUNC_MOD'
CONTAINS
SUBROUTINE calc_contrib_func(n_profile, n_layer, n_cloud_top                   &
      , plev, planck_flux_band, ss_prop, contrib_funci_part                    &
      , contrib_funcf_part, nd_profile, nd_layer)

  USE rad_ccf, ONLY: pi
  USE diffusivity_factor, ONLY: diffusivity_factor_cf
  USE realtype_rd, ONLY: RealK
  USE def_ss_prop, ONLY: str_ss_prop
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE
  
! Dummy arguments.

! Input
  INTEGER, INTENT(IN) ::                                                       &
      n_profile, n_layer, n_cloud_top, nd_profile, nd_layer
!       Dimensions
  REAL (RealK), INTENT(IN) ::                                                  &
      plev(nd_profile, 0:nd_layer)
!       Pressure      
  REAL (RealK), INTENT(IN) ::                                                  &
      planck_flux_band(nd_profile, nd_layer)
!       Planckian
  REAL (RealK), INTENT(OUT) ::                                                 &
      contrib_funci_part(nd_profile, nd_layer)
!       Contribution function (intensity) increment
  REAL (RealK), INTENT(OUT) ::                                                 &
      contrib_funcf_part(nd_profile, nd_layer)
!       Contribution function (flux) increment
  TYPE(str_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

! Local variables
  INTEGER :: i, l
!       Loop counters
  REAL (RealK) ::                                                              &
      tau_clr_f(nd_profile, nd_layer)
!       Clear-sky optical depths for the whole column
  REAL (RealK) ::                                                              &
      tau_clr_sum(nd_profile, 0:nd_layer)
!       Clear-sky optical depths for the whole column from TOA

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CALC_CONTRIB_FUNC'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise optical depth at top-of-atmosphere
  tau_clr_sum(1:n_profile, 0) = 0.0_RealK

! Calculate clear sky optical depth
  DO i=1, n_cloud_top-1
    DO l=1, n_profile
      tau_clr_f(l, i)=ss_prop%tau_clr(l, i)
    END DO
  END DO
  DO i=n_cloud_top, n_layer
    DO l=1, n_profile
      tau_clr_f(l, i)=ss_prop%tau(l, i, 0)
    END DO
  END DO

! Calculate contribution function on each level
  DO i=1, n_layer
    DO l=1, n_profile
!     Optical depth from TOA
      tau_clr_sum(l, i) = SUM(tau_clr_f(l, 1:i))

!     Intensity contribution function
      contrib_funci_part(l, i) = (exp(-tau_clr_sum(l, i-1))                    &
        -exp(-tau_clr_sum(l, i)))*planck_flux_band(l, i)/pi                    &
        /(log10(plev(l,i))-log10(plev(l,i-1)))

!     Flux contribution function
      contrib_funcf_part(l, i) =                                               &
        2.*(exp(-diffusivity_factor_cf*tau_clr_sum(l, i-1))                    &
        -exp(-diffusivity_factor_cf*tau_clr_sum(l, i)))                        &
        *planck_flux_band(l, i)/diffusivity_factor_cf                          &
        /(log10(plev(l,i))-log10(plev(l,i-1)))
    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_contrib_func
END MODULE calc_contrib_func_mod

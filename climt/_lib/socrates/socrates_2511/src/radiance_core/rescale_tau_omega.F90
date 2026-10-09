! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to rescale optical depth and albedo.
!
! Method:
!   The standard rescaling formulae are applied.
!
!- ---------------------------------------------------------------------
MODULE rescale_tau_omega_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'RESCALE_TAU_OMEGA_MOD'
CONTAINS
SUBROUTINE rescale_tau_omega(n_profile                                  &
    , i_layer_first, i_layer_last                                       &
    , tau, omega, forward_scatter                                       &
    , nd_profile, nd_layer, id_1                                        &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , id_1
!       Topmost declared layer for optical properties

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to rescale
    , i_layer_last
!       First layer to rescale
  REAL (RealK), INTENT(IN) ::                                           &
      forward_scatter(nd_profile, id_1: nd_layer)
!       Forward scattering
  REAL (RealK), INTENT(INOUT) ::                                        &
      tau(nd_profile, id_1: nd_layer)                                   &
!       Optical depth
    , omega(nd_profile, id_1: nd_layer)
!       Albedo of single scattering

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='RESCALE_TAU_OMEGA'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile
      tau(l, i)=tau(l, i)*(1.0e+00_RealK                                &
        -omega(l, i)*forward_scatter(l, i))
      omega(l, i)=omega(l, i)*(1.0e+00_RealK-forward_scatter(l, i))     &
        /(1.0e+00_RealK-omega(l, i)*forward_scatter(l, i))
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE rescale_tau_omega
END MODULE rescale_tau_omega_mod

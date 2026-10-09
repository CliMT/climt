! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate fluxes in a column ignoring scattering.
!
!- ---------------------------------------------------------------------
MODULE solver_no_scat_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLVER_NO_SCAT_MOD'
CONTAINS
SUBROUTINE solver_no_scat(n_profile, n_layer                            &
    , trans, s_down, s_up                                               &
    , diffuse_albedo, flux_inc_down, d_planck_flux_surface              &
    , flux_total                                                        &
    , nd_profile, nd_layer)


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer
!       Size allocated for atmospheric layers

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                           &
      trans(nd_profile, nd_layer)                                       &
!       Transmission coefficient
    , s_down(nd_profile, nd_layer)                                      &
!       Downward diffuse source
    , s_up(nd_profile, nd_layer)                                        &
!       Upward diffuse source
    , diffuse_albedo(nd_profile)                                        &
!       Diffuse surface albedo
    , d_planck_flux_surface(nd_profile)                                 &
!       Difference between the Planckian flux at the surface
!       temperature and that of the overlaying air
    , flux_inc_down(nd_profile)
!       Incident total flux

  REAL (RealK), INTENT(OUT) ::                                          &
      flux_total(nd_profile, 2*nd_layer+2)
!       Total flux

! Declaration of local variables.
  INTEGER ::                                                            &
      i, l
!       Loop variables

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLVER_NO_SCAT'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Downward fluxes.
  DO l=1, n_profile
    flux_total(l, 2)=flux_inc_down(l)
  END DO
  DO i=1, n_layer
    DO l=1, n_profile
      flux_total(l, 2*i+2)=s_down(l, i)+trans(l, i)*flux_total(l, 2*i)
    END DO
  END DO

! Upward fluxes.
  DO l=1, n_profile
    flux_total(l, 2*n_layer+1) = &
      (1.0_RealK - diffuse_albedo(l))*d_planck_flux_surface(l)          &
      + diffuse_albedo(l)*flux_total(l, 2*n_layer+2)
  END DO
  DO i=n_layer, 1, -1
    DO l=1, n_profile
      flux_total(l, 2*i-1)=s_up(l, i)+trans(l, i)*flux_total(l, 2*i+1)
    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solver_no_scat
END MODULE solver_no_scat_mod

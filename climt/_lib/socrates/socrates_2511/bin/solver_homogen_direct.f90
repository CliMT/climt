! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate fluxes in a homogeneous column directly.
!
!- ---------------------------------------------------------------------
MODULE solver_homogen_direct_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLVER_HOMOGEN_DIRECT_MOD'
CONTAINS
SUBROUTINE solver_homogen_direct(n_profile, n_layer                     &
    , trans, reflect                                                    &
    , s_down, s_up                                                      &
    , diffuse_albedo                                                    &
    , flux_inc_down                                                     &
    , source_ground                                                     &
    , flux_total                                                        &
    , nd_profile, nd_layer                                              &
    )


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
    , reflect(nd_profile, nd_layer)                                     &
!       Reflection coefficient
    , s_down(nd_profile, nd_layer)                                      &
!       Downward diffuse source
    , s_up(nd_profile, nd_layer)                                        &
!       Upward diffuse source
    , diffuse_albedo(nd_profile)                                        &
!       Diffuse surface albedo
    , flux_inc_down(nd_profile)                                         &
!       Incident total flux
    , source_ground(nd_profile)
!       Source from ground

  REAL (RealK), INTENT(OUT) ::                                          &
      flux_total(nd_profile, 2*nd_layer+2)
!       Total flux

! Declaration of local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) ::                                                       &
      alpha(nd_profile)                                                 &
!       Combined albedo of lower layers
    , beta(nd_profile, nd_layer)                                        &
!       Working array
    , gamma(nd_profile, nd_layer)                                       &
!       Working array
    , h(nd_profile, nd_layer)                                           &
!       Working array
    , s_up_prime(nd_profile)
!       Modified upward source function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLVER_HOMOGEN_DIRECT'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Eliminating loop:
  DO i=n_layer+1, 1, -1
    IF ( i < n_layer+1 ) THEN
      DO l=1, n_profile
        beta(l, i)=1.0e+00_RealK/(1.0e+00_RealK-alpha(l)*reflect(l, i))
        gamma(l, i)=alpha(l)*trans(l, i)
        h(l, i)=s_up_prime(l)+alpha(l)*s_down(l, i)
        alpha(l)=reflect(l, i)+beta(l, i)*gamma(l, i)*trans(l, i)
        s_up_prime(l)=s_up(l, i)+beta(l, i)*trans(l, i)*h(l, i)
      END DO
    ELSE
!     Initialization at the bottom for upward elimination:
      DO l=1, n_profile
        alpha(l)=diffuse_albedo(l)
        s_up_prime(l)=source_ground(l)
      END DO
    END IF
  END DO

! Backward substitution:
  DO i=0, n_layer
    IF ( i > 0 ) THEN
      DO l=1, n_profile
!       Upward flux
        flux_total(l, 2*i+1)                                              &
          =beta(l, i)*(h(l, i)+gamma(l, i)*flux_total(l, 2*i))
!       Downward flux
        flux_total(l, 2*i+2)=s_down(l, i)                                 &
          +trans(l, i)*flux_total(l, 2*i)                                 &
          +reflect(l, i)*flux_total(l, 2*i+1)
      END DO
    ELSE
!     Initialize for backward substitution.
      DO l=1, n_profile
        flux_total(l, 2)=flux_inc_down(l)
        flux_total(l, 1)=alpha(l)*flux_total(l, 2)+s_up_prime(l)
      END DO
    END IF
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solver_homogen_direct
END MODULE solver_homogen_direct_mod

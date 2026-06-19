! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the pentadiagonal matrix for the fluxes.
!
!- ---------------------------------------------------------------------
MODULE set_matrix_pentadiagonal_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_MATRIX_PENTADIAGONAL_MOD'
CONTAINS
SUBROUTINE set_matrix_pentadiagonal(n_profile, n_layer                  &
     , trans, reflect                                                   &
     , s_down, s_up                                                     &
     , diffuse_albedo                                                   &
     , flux_inc_down                                                    &
     , source_ground                                                    &
     , a5, b                                                            &
     , nd_profile, nd_layer                                             &
     )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       SIze allocated for atmospheric profiles
    , nd_layer
!       SIze allocated for atmospheric layers


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
      a5(nd_profile, 5, 2*nd_layer+2)                                   &
!       Pentadiagonal matrix
    , b(nd_profile, 2*nd_layer+2)
!       Source terms for equations

! Declaration of local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_MATRIX_PENTADIAGONAL'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The top boundary condition:
  DO l=1, n_profile
    a5(l, 4, 2)=0.0e+00_RealK
    a5(l, 3, 2)=1.0e+00_RealK
    a5(l, 2, 2)=0.0e+00_RealK
    a5(l, 1, 2)=0.0e+00_RealK
    b(l, 2)=flux_inc_down(l)
  END DO

! Interior rows: odd and even rows correspond to different boundary
! conditions.
  DO i=1, n_layer
    DO l=1, n_profile

      a5(l, 5, 2*i-1)=0.0e+00_RealK
      a5(l, 4, 2*i-1)=0.0e+00_RealK
      a5(l, 3, 2*i-1)=-1.0e+00_RealK
      a5(l, 2, 2*i-1)=reflect(l, i)
      a5(l, 1, 2*i-1)=trans(l, i)
      b(l, 2*i-1)=-s_up(l, i)

      a5(l, 5, 2*i+2)=trans(l, i)
      a5(l, 4, 2*i+2)=reflect(l, i)
      a5(l, 3, 2*i+2)=-1.0e+00_RealK
      a5(l, 2, 2*i+2)=0.0e+00_RealK
      a5(l, 1, 2*i+2)=0.0e+00_RealK
      b(l, 2*i+2)=-s_down(l, i)

    END DO
  END DO

! The surface boundary condition:
  DO l=1, n_profile
    a5(l, 5, 2*n_layer+1)=0.0e+00_RealK
    a5(l, 4, 2*n_layer+1)=0.0e+00_RealK
    a5(l, 3, 2*n_layer+1)=1.0e+00_RealK
    a5(l, 2, 2*n_layer+1)=-diffuse_albedo(l)
    b(l, 2*n_layer+1)=source_ground(l)
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_matrix_pentadiagonal
END MODULE set_matrix_pentadiagonal_mod

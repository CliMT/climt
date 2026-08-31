! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set clear-sky solar phase function.
!
! Method:
!   The arrays of clear-sky forward scattering and the solar
!   phase function at the top of the column and of these
!   same properties from the total list lower down
!   are combined to give unified arrays of clear-sky
!   optical properties.
!
!- ---------------------------------------------------------------------
MODULE copy_clr_sol_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'COPY_CLR_SOL_MOD'
CONTAINS
SUBROUTINE copy_clr_sol(n_profile, n_layer, n_cloud_top                 &
    , n_direction, l_rescale                                            &
    , forward_scatter_clr, phase_fnc_solar_clr                          &
    , forward_scatter, phase_fnc_solar                                  &
    , forward_scatter_clr_f                                             &
    , phase_fnc_solar_clr_f                                             &
!                   Sizes of arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_direction           &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Sizes of arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_layer_clr                                                      &
!       Size allocated for totally clear layers
    , id_ct                                                             &
!       Topmost declared layer for cloudy optical properties
    , nd_direction
!       SIze allocated for viewing directions

!                 Atmospheric properties
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of atmospheric layers
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_direction
!       Number of terms in the phase function
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale
!       Flag for rescaling of the optical properties

!                 Optical properties
  REAL (RealK), INTENT(IN) ::                                           &
      forward_scatter_clr(nd_profile, nd_layer_clr)                     &
!       Forward scattering in the totally clear region
    , phase_fnc_solar_clr(nd_profile, nd_layer_clr, nd_direction)
!       Phase function in totally clear region
  REAL (RealK), INTENT(IN) ::                                           &
      forward_scatter(nd_profile, id_ct: nd_layer)                      &
!       Forward scattering in the cloudy regions
    , phase_fnc_solar(nd_profile, id_ct: nd_layer, nd_direction)
!       Phase function restricted to clear-sky regions

!                 Single scattering properties
  REAL (RealK), INTENT(OUT) ::                                          &
      forward_scatter_clr_f(nd_profile, nd_layer)                       &
!       Forward scattering expanded to the whole column
    , phase_fnc_solar_clr_f(nd_profile, nd_layer, nd_direction)
!       Phase function expanded to the whole column



! Local variables.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , k
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_CLR_SOL'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Above cloud top.
  DO i=1, n_cloud_top-1
    IF (l_rescale) THEN
      DO l=1, n_profile
        forward_scatter_clr_f(l, i)=forward_scatter_clr(l, i)
      END DO
    END IF
    DO k=1, n_direction
      DO l=1, n_profile
        phase_fnc_solar_clr_f(l, i, k)=phase_fnc_solar_clr(l, i, k)
      END DO
    END DO
  END DO

! Below cloud top.
  DO i=n_cloud_top, n_layer
    IF (l_rescale) THEN
      DO l=1, n_profile
        forward_scatter_clr_f(l, i)=forward_scatter(l, i)
      END DO
    END IF
    DO k=1, n_direction
      DO l=1, n_profile
        phase_fnc_solar_clr_f(l, i, k)=phase_fnc_solar(l, i, k)
      END DO
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE copy_clr_sol
END MODULE copy_clr_sol_mod

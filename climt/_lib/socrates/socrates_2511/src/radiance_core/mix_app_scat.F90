! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve for fluxes treating scattering approximately.
!
! Method:
!   The routine is applicable in the infra-red. downward
!   differential fluxes are calculated first assuming that the
!   upward differential fluxes are 0. Upward fluxes are then
!   calculated using the previously calculated downward fluxes
!   in the reflected terms.
!
!- ---------------------------------------------------------------------
MODULE mix_app_scat_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MIX_APP_SCAT_MOD'
CONTAINS
SUBROUTINE mix_app_scat(n_profile, n_layer, n_cloud_top                 &
     , t_free, r_free, s_down_free, s_up_free                           &
     , t_cloud, r_cloud, s_down_cloud, s_up_cloud                       &
     , g_ff, g_fc, g_cf, g_cc                                           &
     , b_ff, b_fc, b_cf, b_cc                                           &
     , flux_inc_down                                                    &
     , source_ground, albedo_surface_diff                               &
     , flux_diffuse                                                     &
     , nd_profile, nd_layer, id_ct                                      &
     )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_ct
!       Topmost declared cloudy layer


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top
!       Topmost cloudy layer
  REAL (RealK), INTENT(IN) ::                                           &
      t_free(nd_profile, nd_layer)                                      &
!       Free transmission
    , r_free(nd_profile, nd_layer)                                      &
!       Free reflection
    , s_down_free(nd_profile, nd_layer)                                 &
!       Free downward source function
    , s_up_free(nd_profile, nd_layer)                                   &
!       Free upward source function
    , t_cloud(nd_profile, nd_layer)                                     &
!       Cloudy transmission
    , r_cloud(nd_profile, nd_layer)                                     &
!       Cloudy reflection
    , s_down_cloud(nd_profile, nd_layer)                                &
!       Downward cloudy source function
    , s_up_cloud(nd_profile, nd_layer)
!       Upward cloudy source function
  REAL (RealK), INTENT(IN) ::                                           &
      g_ff(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , g_fc(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , g_cf(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , g_cc(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , b_ff(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , b_fc(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , b_cf(nd_profile, id_ct-1: nd_layer)                               &
!       Energy transfer coefficient
    , b_cc(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient
  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_down(nd_profile)                                         &
!       Incident diffuse flux
    , source_ground(nd_profile)                                         &
!       Source from ground
    , albedo_surface_diff(nd_profile)
!       Diffuse albedo
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_diffuse(nd_profile, 2*nd_layer+2)
!       Diffuse flux

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) ::                                                       &
      flux_down(nd_profile, 0: nd_layer)                                &
!       Downward fluxes outside clouds just below i'th level
    , flux_down_cloud(nd_profile, 0: nd_layer)                          &
!       Downward fluxes inside clouds just below i'th level
    , flux_up(nd_profile, 0: nd_layer)                                  &
!       Upward fluxes outside clouds just above i'th level
    , flux_up_cloud(nd_profile, 0: nd_layer)                            &
!       Upward fluxes inside clouds just above i'th level
    , flux_propagated                                                   &
!       Temporary propagated flux outside cloud
    , flux_propagated_cloud                                             &
!       Temporary propagated flux inside cloud
    , flux_cloud_top(nd_profile)
!       Total downward flux at top of cloud

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MIX_APP_SCAT'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The arrays flux_down and flux_up will eventually contain the total
! fluxes, but initially they are used for the clear fluxes.
! Note that downward fluxes refer to values just below the interface
! and upward fluxes to values just above it.


! Downward flux:

! Region above clouds:
  DO l=1, n_profile
    flux_down(l, 0)=flux_inc_down(l)
  END DO
  DO i=1, n_cloud_top-1
    DO l=1, n_profile
      flux_down(l, i)=t_free(l, i)*flux_down(l, i-1)                    &
        +s_down_free(l, i)
    END DO
  END DO
  DO l=1, n_profile
    flux_cloud_top(l)=flux_down(l, n_cloud_top-1)
  END DO

! Region of clouds:
  DO l=1, n_profile
    flux_down(l, n_cloud_top-1)                                         &
      =g_ff(l, n_cloud_top-1)*flux_cloud_top(l)
    flux_down_cloud(l, n_cloud_top-1)                                   &
      =g_fc(l, n_cloud_top-1)*flux_cloud_top(l)
  END DO

  DO i=n_cloud_top, n_layer-1
    DO l=1, n_profile

!     Propagate downward fluxes through the layer.
      flux_propagated=t_free(l, i)*flux_down(l, i-1)                    &
        +s_down_free(l, i)
      flux_propagated_cloud=t_cloud(l, i)*flux_down_cloud(l, i-1)       &
        +s_down_cloud(l, i)
!     Transfer downward fluxes across the interface.
      flux_down(l, i)                                                   &
        =g_ff(l, i)*flux_propagated                                     &
        +g_cf(l, i)*flux_propagated_cloud
      flux_down_cloud(l, i)                                             &
        =g_cc(l, i)*flux_propagated_cloud                               &
        +g_fc(l, i)*flux_propagated

    END DO
  END DO

! Propagate across the bottom layer, but without transferring
! across the surface and form the reflected beams.
  DO l=1, n_profile
!   Propagate downward fluxes through the layer.
    flux_down(l, n_layer)                                               &
      =t_free(l, n_layer)*flux_down(l, n_layer-1)                       &
      +s_down_free(l, n_layer)
    flux_down_cloud(l, n_layer)                                         &
      =t_cloud(l, n_layer)*flux_down_cloud(l, n_layer-1)                &
      +s_down_cloud(l, n_layer)
    flux_up(l, n_layer)                                                 &
      =albedo_surface_diff(l)*flux_down(l, n_layer)                     &
      +b_ff(l, n_layer)*source_ground(l)
    flux_up_cloud(l, n_layer)                                           &
      =albedo_surface_diff(l)*flux_down_cloud(l, n_layer)               &
      +b_cf(l, n_layer)*source_ground(l)
  END DO


! Calculate the upward fluxes using the previous downward fluxes
! to approximate the scattering term.
  DO i=n_layer, n_cloud_top, -1
    DO l=1, n_profile

!     Propagate upward fluxes through the layer.
      flux_propagated=t_free(l, i)*flux_up(l, i)+s_up_free(l, i)        &
        +r_free(l, i)*flux_down(l, i-1)
      flux_propagated_cloud=t_cloud(l, i)*flux_up_cloud(l, i)           &
        +s_up_cloud(l, i)+r_cloud(l, i)*flux_down_cloud(l, i-1)
!     Transfer upward fluxes across the interface.
      flux_up(l, i-1)=b_ff(l, i-1)*flux_propagated                      &
        +b_fc(l, i-1)*flux_propagated_cloud
      flux_up_cloud(l, i-1)=b_cc(l, i-1)*flux_propagated_cloud          &
        +b_cf(l, i-1)*flux_propagated

    END DO
  END DO

! Continue through the region above clouds.
  DO i=n_cloud_top-1, 1, -1
    DO l=1, n_profile
      flux_up(l, i-1)=t_free(l, i)*flux_up(l,i)+s_up_free(l, i)         &
        +r_free(l, i)*flux_down(l, i-1)
    END DO
  END DO



! Calculate the overall flux.
  DO i=0, n_cloud_top-2
    DO l=1, n_profile
      flux_diffuse(l, 2*i+1)=flux_up(l, i)
      flux_diffuse(l, 2*i+2)=flux_down(l, i)
    END DO
  END DO
  DO i=n_cloud_top-1, n_layer
    DO l=1, n_profile
      flux_diffuse(l, 2*i+1)=flux_up(l, i)+flux_up_cloud(l, i)
      flux_diffuse(l, 2*i+2)=flux_down(l, i)                            &
        +flux_down_cloud(l, i)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE mix_app_scat
END MODULE mix_app_scat_mod

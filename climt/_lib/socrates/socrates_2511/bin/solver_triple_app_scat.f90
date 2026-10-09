! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve for triple overlaps with approximate scattering.
!
! Method:
!   The flux is propagated downwards, ignoring reflection terms.
!   since the routine uses differential fluxes, this effectively
!   treats the upward flux as Planckian at this point. Upward
!   fluxes are calculated using the newly available approximate
!   downward fluxes in the reflected terms.
!
!- ---------------------------------------------------------------------
MODULE solver_triple_app_scat_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLVER_TRIPLE_APP_SCAT_MOD'
CONTAINS
SUBROUTINE solver_triple_app_scat(n_profile, n_layer, n_cloud_top       &
     , t, r, s_down, s_up                                               &
     , t_strat, r_strat, s_down_strat, s_up_strat                       &
     , t_conv, r_conv, s_down_conv, s_up_conv                           &
     , v11, v12, v13, v21, v22, v23, v31, v32, v33                      &
     , u11, u12, u13, u21, u22, u23, u31, u32, u33                      &
     , flux_inc_down                                                    &
     , source_ground_free, source_ground_strat                          &
     , source_ground_conv, albedo_surface_diff                          &
     , flux_total                                                       &
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
      t(nd_profile, nd_layer)                                           &
!       Clear-sky transmission
    , r(nd_profile, nd_layer)                                           &
!       Clear-sky reflection
    , s_down(nd_profile, nd_layer)                                      &
!       Clear-sky downward source function
    , s_up(nd_profile, nd_layer)                                        &
!       Clear-sky upward source function
    , t_strat(nd_profile, nd_layer)                                     &
!       Stratfiform transmission
    , r_strat(nd_profile, nd_layer)                                     &
!       Stratfiform reflection
    , s_down_strat(nd_profile, nd_layer)                                &
!       Downward stratfiform source function
    , s_up_strat(nd_profile, nd_layer)                                  &
!       Upward stratfiform source function
    , t_conv(nd_profile, nd_layer)                                      &
!       Convective transmission
    , r_conv(nd_profile, nd_layer)                                      &
!       Convective reflection
    , s_down_conv(nd_profile, nd_layer)                                 &
!       Downward convective source function
    , s_up_conv(nd_profile, nd_layer)
!       Upward convective source function
  REAL (RealK), INTENT(IN) ::                                           &
      v11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v13(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v22(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v23(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v31(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v32(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for downward radiation
    , v33(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient for downward radiation
  REAL (RealK) ::                                                       &
      u11(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u12(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u13(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u21(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u22(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u23(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u31(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u32(nd_profile, id_ct-1: nd_layer)                                &
!       Energy transfer coefficient for upward radiation
    , u33(nd_profile, id_ct-1: nd_layer)
!       Energy transfer coefficient for upward radiation
  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_down(nd_profile)                                         &
!       Incident flux
    , source_ground_free(nd_profile)                                    &
!       Source from ground (clear sky)
    , source_ground_strat(nd_profile)                                   &
!       Source from ground (cloudy region)
    , source_ground_conv(nd_profile)                                    &
!       Source from ground (cloudy region)
    , albedo_surface_diff(nd_profile)
!       Diffuse albedo
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_total(nd_profile, 2*nd_layer+2)
!       Total flux

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable


! Temporary fluxes
  REAL (RealK) ::                                                       &
      flux_down_1(nd_profile, 0: nd_layer)                              &
!       Downward fluxes outside clouds just below i'th level
    , flux_down_2(nd_profile, 0: nd_layer)                              &
!       Downward fluxes inside clouds just below i'th level
    , flux_down_3(nd_profile, 0: nd_layer)                              &
!       Downward fluxes inside clouds just below i'th level
    , flux_up_1(nd_profile, 0: nd_layer)                                &
!       Upward fluxes outside clouds just above i'th level
    , flux_up_2(nd_profile, 0: nd_layer)                                &
!       Upward fluxes inside clouds just above i'th level
    , flux_up_3(nd_profile, 0: nd_layer)                                &
!       Upward fluxes inside clouds just above i'th level
    , flux_propag_1(nd_profile)                                         &
!       Temporary fluxes for propagation across layers
    , flux_propag_2(nd_profile)                                         &
!       Temporary fluxes for propagation across layers
    , flux_propag_3(nd_profile)
!       Temporary fluxes for propagation across layers

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLVER_TRIPLE_APP_SCAT'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The arrays flux_down and flux_up will eventually contain the total
! fluxes, but initially they are used for the clear fluxes.
! Note that downward fluxes refer to values just below the interface
! and upward fluxes to values just above it.


! Downward flux:

! Region above clouds:
  DO l=1, n_profile
    flux_total(l, 2)=flux_inc_down(l)
  END DO
  DO i=1, n_cloud_top-1
    DO l=1, n_profile
      flux_total(l, 2*i+2)=t(l, i)*flux_total(l, 2*i)                   &
        +s_down(l, i)
    END DO
  END DO

! Pass into the cloudy region. here, downward fluxes hold values
! just below the level and upward fluxes the values just above it.
! Thus the fluxes impinging on the layer are held.
  i=n_cloud_top-1
  DO l=1, n_profile
    flux_down_1(l, i)=v11(l, i)*flux_total(l, 2*i+2)
    flux_down_2(l, i)=v21(l, i)*flux_total(l, 2*i+2)
    flux_down_3(l, i)=v31(l, i)*flux_total(l, 2*i+2)
  END DO

  DO i=n_cloud_top, n_layer-1
    DO l=1, n_profile

!     Propagte the flux across the layer.
      flux_propag_1(l)=t(l, i)*flux_down_1(l, i-1)                      &
        +s_down(l, i)
      flux_propag_2(l)=t_strat(l, i)*flux_down_2(l, i-1)                &
        +s_down_strat(l, i)
      flux_propag_3(l)=t_conv(l, i)*flux_down_3(l, i-1)                 &
        +s_down_conv(l, i)

!     Transfer across the interface.
      flux_down_1(l, i)=v11(l, i)*flux_propag_1(l)                      &
        +v12(l, i)*flux_propag_2(l)                                     &
        +v13(l, i)*flux_propag_3(l)
      flux_down_2(l, i)=v21(l, i)*flux_propag_1(l)                      &
        +v22(l, i)*flux_propag_2(l)                                     &
        +v23(l, i)*flux_propag_3(l)
      flux_down_3(l, i)=v31(l, i)*flux_propag_1(l)                      &
        +v32(l, i)*flux_propag_2(l)                                     &
        +v33(l, i)*flux_propag_3(l)

    END DO
  END DO

! Propagate across the bottom layer and form the reflected beam.
! We do not transfer fluxes across the bottom interface, so as
! to make the reflection consistent between regions.
  DO l=1, n_profile

!   Propagte the flux through the layer.
    flux_down_1(l, n_layer)                                             &
      =t(l, n_layer)*flux_down_1(l, n_layer-1)                          &
      +s_down(l, n_layer)
    flux_down_2(l, n_layer)                                             &
      =t_strat(l, n_layer)*flux_down_2(l, n_layer-1)                    &
      +s_down_strat(l, n_layer)
    flux_down_3(l, n_layer)                                             &
      =t_conv(l, n_layer)*flux_down_3(l, n_layer-1)                     &
      +s_down_conv(l, n_layer)

!   Reflect from the surface.
    flux_up_1(l, n_layer)                                               &
      =albedo_surface_diff(l)*flux_down_1(l, n_layer)                   &
      +source_ground_free(l)
    flux_up_2(l, i)                                                     &
      =albedo_surface_diff(l)*flux_down_2(l, n_layer)                   &
      +source_ground_strat(l)
    flux_up_3(l, i)                                                     &
      =albedo_surface_diff(l)*flux_down_3(l, n_layer)                   &
      +source_ground_conv(l)

!   Propagate across the bottom layer.
    flux_propag_1(l)                                                    &
      =t(l, n_layer)*flux_up_1(l, n_layer)+s_up(l, n_layer)             &
      +r(l, n_layer)*flux_down_1(l, n_layer-1)
    flux_propag_2(l)                                                    &
      =t_strat(l, n_layer)*flux_up_2(l, n_layer)                        &
      +s_up_strat(l, n_layer)                                           &
      +r_strat(l, n_layer)*flux_down_2(l, n_layer-1)
    flux_propag_3(l)                                                    &
      =t_conv(l, n_layer)*flux_up_3(l, n_layer)                         &
      +s_up_conv(l, n_layer)                                            &
      +r_conv(l, n_layer)*flux_down_3(l, n_layer-1)

  END DO



! Work back up through the column assigning the upward fluxes.
  DO i=n_layer-1, n_cloud_top, -1
    DO l=1, n_profile

      flux_up_1(l, i)=u11(l, i)*flux_propag_1(l)                        &
        +u12(l, i)*flux_propag_2(l)                                     &
        +u13(l, i)*flux_propag_3(l)
      flux_up_2(l, i)=u21(l, i)*flux_propag_1(l)                        &
        +u22(l, i)*flux_propag_2(l)                                     &
        +u23(l, i)*flux_propag_3(l)
      flux_up_3(l, i)=u31(l, i)*flux_propag_1(l)                        &
        +u32(l, i)*flux_propag_2(l)                                     &
        +u33(l, i)*flux_propag_3(l)

      flux_propag_1(l)=t(l, i)*flux_up_1(l, i)+s_up(l, i)               &
        +r(l, i)*flux_down_1(l, i-1)
      flux_propag_2(l)=t_strat(l, i)*flux_up_2(l, i)                    &
        +s_up_strat(l, i)+r_strat(l, i)*flux_down_2(l, i-1)
      flux_propag_3(l)=t_conv(l, i)*flux_up_3(l, i)                     &
        +s_up_conv(l, i)+r_conv(l, i)*flux_down_3(l, i-1)

    END DO
  END DO

! Propagate into the cloud-free region.
  i=n_cloud_top-1
  DO l=1, n_profile
    flux_total(l, 2*i+1)=flux_propag_1(l)+flux_propag_2(l)              &
      +flux_propag_3(l)
  END DO

! Continue through the upper cloudy layers.
  DO i=n_cloud_top-1, 1, -1
    DO l=1, n_profile
      flux_total(l, 2*i-1)=t(l, i)*flux_total(l, 2*i+1)                 &
        +r(l, i)*flux_total(l, 2*i)+s_up(l, i)
    END DO
  END DO

! Assign the total fluxes on the intermediate cloudy layers.
  DO i=n_cloud_top, n_layer
    DO l=1, n_profile
      flux_total(l, 2*i+1)=flux_up_1(l, i)+flux_up_2(l, i)              &
        +flux_up_3(l, i)
      flux_total(l, 2*i+2)=flux_down_1(l, i)+flux_down_2(l, i)          &
        +flux_down_3(l, i)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solver_triple_app_scat
END MODULE solver_triple_app_scat_mod

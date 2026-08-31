! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate fluxes using gaussian quadrature.
!
! Method:
!   Fluxes are calculated by using gaussian quadrature for
!   the angular integration. This is not a full implementation
!   of gaussian quadrature for multiple scattering, but is
!   intended only for non-scattering calculations in the
!   infra-red. In this case, the fluxes can be calculated as
!   a weighted sum of two-stream fluxes where the diffusivity
!   factors for the two-stream approximations are determined
!   from the gaussian points.
!
!- ---------------------------------------------------------------------
MODULE gauss_angle_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'GAUSS_ANGLE_MOD'
CONTAINS
SUBROUTINE gauss_angle(n_profile, n_layer                               &
     , n_order_gauss                                                    &
     , tau                                                              &
     , flux_inc_down                                                    &
     , diff_planck, source_ground, albedo_surface_diff                  &
     , flux_diffuse                                                     &
     , nd_profile, nd_layer                                             &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE gaussian_weight_pcf, ONLY : gauss_weight, gauss_point
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE monochromatic_ir_radiance_mod, ONLY: monochromatic_ir_radiance

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer
!       Maximum number of layers


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_order_gauss
!       Order of gaussian integration
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)                                         &
!       Optical depth
    , albedo_surface_diff(nd_profile)                                   &
!       Diffuse albedo
    , flux_inc_down(nd_profile)                                         &
!       Incident total flux
    , diff_planck(nd_profile, nd_layer)                                 &
!       Difference in pi*Planckian function
    , source_ground(nd_profile)
!       Ground source function
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_diffuse(nd_profile, 2*nd_layer+2)
!       Diffuse fluxes

! Local variabales.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , k
!       Loop variable
  REAL (RealK) ::                                                       &
      flux_stream(nd_profile, 2*nd_layer+2)                             &
!       Flux in stream
    , secant_ray                                                        &
!       Secant of angle with vertical
    , diff_planck_rad(nd_profile, nd_layer)                             &
!       Difference in pi*Planckian function
    , source_ground_rad(nd_profile)                                     &
!       Ground source function
    , radiance_inc(nd_profile)                                          &
!       Incidnet radiance
    , weight_stream
!       Weighting for stream

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='GAUSS_ANGLE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the source function.
  DO l=1, n_profile
    source_ground_rad(l)=source_ground(l)/pi
    radiance_inc(l)=flux_inc_down(l)/pi
  END DO
  DO i=1, n_layer
    DO l=1, n_profile
      diff_planck_rad(l, i)=diff_planck(l, i)/pi
    END DO
  END DO
  DO i=1, 2*n_layer+2
    DO l=1, n_profile
      flux_diffuse(l, i)=0.0
    END DO
  END DO

! Calculate the fluxes with a number of diffusivity factors
! and sum the results.
  DO k=1, n_order_gauss
    secant_ray=2.0e+00_RealK                                            &
      /(gauss_point(k, n_order_gauss)+1.0e+00_RealK)

!   Calculate the radiance at this angle.
    CALL monochromatic_ir_radiance(n_profile, n_layer                   &
      , tau                                                             &
      , radiance_inc                                                    &
      , diff_planck_rad, source_ground_rad, albedo_surface_diff         &
      , secant_ray                                                      &
      , flux_stream                                                     &
      , nd_profile, nd_layer                                            &
      )

!   Augment the flux by the amount in this stream.
    weight_stream=5.0e-01_RealK*pi*gauss_weight(k, n_order_gauss)       &
      *(gauss_point(k, n_order_gauss)+1.0e+00_RealK)
    DO i=1, 2*n_layer+2
      DO l=1, n_profile
        flux_diffuse(l, i)=flux_diffuse(l, i)                           &
          +weight_stream*flux_stream(l, i)
      END DO
    END DO

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE gauss_angle
END MODULE gauss_angle_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to increment upward fluxes on a tiled surface.
!
! Method:
!   The arrays holding the local cumulative fluxes or radiances
!   on each tile are incremented by the variables suffixed
!   with _INCR, multiplied by appropriate weights. The routine
!   can be called to initialize fluxes.
!
!- ---------------------------------------------------------------------
MODULE augment_tiled_radiance_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'AUGMENT_TILED_RADIANCE_MOD'
CONTAINS
SUBROUTINE augment_tiled_radiance(control, spectrum, radout             &
    , i_band, iex, iex_minor                                            &
    , n_point_tile, n_tile, list_tile                                   &
    , l_initial_channel_tile                                            &
    , weight_incr, weight_blue_incr, weight_sub_band_incr               &
!                 Surface characteristics
    , rho_alb                                                           &
!                 Increments to radiances
    , flux_direct_incr, flux_down_incr                                  &
    , planck_flux_tile, planck_flux_air                                 &
!                 Dimensions
    , nd_flux_profile, nd_point_tile, nd_tile                           &
    , nd_brdf_basis_fnc, nd_channel, nd_abs                             &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_spectrum, ONLY: StrSpecData
  USE def_out,     ONLY: StrOut
  USE rad_pcf, ONLY: i_err_fatal, ip_infra_red, ip_ir_gauss, ip_solar,  &
                     ip_spherical_harmonic, ip_surf_alb_diff,           &
                     ip_surf_alb_dir, ip_two_stream
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl), INTENT(IN)     :: control

! Spectral data:
  TYPE(StrSpecData), INTENT(IN) :: spectrum

! Output fields:
  TYPE(StrOut),  INTENT(INOUT)  :: radout

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_flux_profile                                                   &
!       Size allocated for points where fluxes are calculated
    , nd_point_tile                                                     &
!       Size allocated for points where the surface is tiled
    , nd_tile                                                           &
!       Size allocated for surface tiles
    , nd_brdf_basis_fnc                                                 &
!       Size allocated for BRDF basis functions
    , nd_channel                                                        &
!       Size allocated for output channels
    , nd_abs
!       Size allocated for absorbers

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      i_band, iex                                                       &
!       Band and k-term being considered
    , iex_minor(nd_abs)
!       k-term of minor gases (only used with exact_major overlap)

  INTEGER, INTENT(IN) ::                                                &
      n_point_tile                                                      &
!       Number of points where the surface is tiled
    , n_tile                                                            &
!       Number of tiles used
    , list_tile(nd_point_tile)
!       List of tiled points

  LOGICAL, INTENT(INOUT) ::                                             &
      l_initial_channel_tile(nd_channel)
!       Perform initialisation instead of incrementing

  REAL (RealK), INTENT(IN) ::                                           &
      weight_incr, weight_blue_incr, weight_sub_band_incr
!       Weight to apply to increments

!                 Surface Characteristics
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb(nd_point_tile, nd_brdf_basis_fnc, nd_tile)
!       Weighting functions for BRDFs

!                 Increments to Fluxes
  REAL (RealK), INTENT(IN) ::                                           &
      flux_direct_incr(nd_flux_profile)                                 &
!       Increment to mean direct flux
    , flux_down_incr(nd_flux_profile)
!       Increment to total downward flux

!                 Planckian Fluxes
  REAL (RealK), INTENT(IN) ::                                           &
      planck_flux_tile(nd_point_tile, nd_tile)                          &
!       Local Planckian flux emitted from each tile
    , planck_flux_air(nd_flux_profile)
!       Hemispheric Planckian flux at the temperature of the air


! Local arguments.
  INTEGER :: l, ll, k, i_sub_k, i_sub, i_channel, i_abs, i_gas, i_k_sub, i_k
!       Loop variables
  REAL (RealK) :: weight_channel_incr
!       Weight to apply to channel increments

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  INTEGER :: ierr
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'AUGMENT_TILED_RADIANCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (control%l_map_sub_bands) THEN
    ! Increment the fluxes with sub-bands mapping to channels
    IF (iex == 0) THEN
      ! The k-term is not set so assume a single sub-band for the band
      i_sub = spectrum%map%list_sub_band_k(1, 1, i_band)
      i_channel = control%map_channel(i_sub)
      weight_channel_incr = weight_sub_band_incr
      CALL augment_channel_tile()
    ELSE
      ! Increment sub-bands for this k-term
      DO i_sub_k=1, spectrum%map%n_sub_band_k(iex, i_band)
        i_sub = spectrum%map%list_sub_band_k(i_sub_k, iex, i_band)
        i_channel = control%map_channel(i_sub)
        weight_channel_incr = weight_sub_band_incr &
          * spectrum%map%weight_sub_band_k(i_sub_k, iex, i_band)
        DO i_abs=2, spectrum%gas%n_band_absorb(i_band)
          IF (iex_minor(i_abs) > 0) THEN
            ! For the exact_major overlap method the minor gas k-term may
            ! be specified allowing the use of the sub-band weights for
            ! exact overlap with the major gas.
            i_gas = spectrum%gas%index_absorb(i_abs, i_band)
            DO i_k_sub=1, spectrum%map%n_k_sub_band(i_gas, i_sub)
              i_k = spectrum%map%list_k_sub_band(i_k_sub, i_gas, i_sub)
              IF (i_k == iex_minor(i_abs)) THEN
                weight_channel_incr = weight_channel_incr &
                  * spectrum%map%weight_k_sub_band(i_k_sub, i_gas, i_sub) &
                  / spectrum%map%weight_k_major(i_k, i_abs, iex, i_band)
              END IF
            END DO
          END IF
        END DO
        CALL augment_channel_tile()
      END DO
    END IF
  ELSE
    ! Increment the fluxes with bands mapping to channels
    i_channel = control%map_channel(i_band)
    weight_channel_incr = weight_incr
    CALL augment_channel_tile()
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


CONTAINS


! Increment the tiled fluxes for each channel.
! An internal subroutine is used to allow multiple channels to be updated
! where sub-bands map to channels.
SUBROUTINE augment_channel_tile()
  IMPLICIT NONE

  IF (l_initial_channel_tile(i_channel)) THEN
!   Initialization of the radiance field takes place here.
    IF ( (control%i_angular_integration == ip_two_stream).OR. &
         (control%i_angular_integration == ip_ir_gauss) ) THEN

!     Initialize the actual fluxes.
      IF (control%isolir == ip_solar) THEN

        IF (control%l_spherical_solar) THEN
          DO k=1, n_tile
            DO ll=1, n_point_tile
              l=list_tile(ll)
              radout%flux_up_tile(ll, k, i_channel) &
                = weight_channel_incr &
                *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                + rho_alb(ll, ip_surf_alb_dir, k)*flux_direct_incr(l))
            END DO
          END DO
          IF (control%l_blue_flux_surf .AND. .NOT.control%l_map_sub_bands) THEN
            DO k=1, n_tile
              DO ll=1, n_point_tile
                l=list_tile(ll)
                radout%flux_up_blue_tile(ll, k, i_channel) &
                  = weight_blue_incr &
                  *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                  + rho_alb(ll, ip_surf_alb_dir, k)*flux_direct_incr(l))
              END DO
            END DO
          END IF
        ELSE
          DO k=1, n_tile
            DO ll=1, n_point_tile
              l=list_tile(ll)
              radout%flux_up_tile(ll, k, i_channel) &
                = weight_channel_incr &
                *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                +(rho_alb(ll, ip_surf_alb_dir, k) &
                - rho_alb(ll, ip_surf_alb_diff, k)) &
                * flux_direct_incr(l))
            END DO
          END DO
          IF (control%l_blue_flux_surf .AND. .NOT.control%l_map_sub_bands) THEN
            DO k=1, n_tile
              DO ll=1, n_point_tile
                l=list_tile(ll)
                radout%flux_up_blue_tile(ll, k, i_channel) &
                  = weight_blue_incr &
                  *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                  +(rho_alb(ll, ip_surf_alb_dir, k) &
                  - rho_alb(ll, ip_surf_alb_diff, k)) &
                  * flux_direct_incr(l))
              END DO
            END DO
          END IF
        END IF

      ELSE IF (control%isolir == ip_infra_red) THEN

        DO k=1, n_tile
          DO ll=1, n_point_tile
            l=list_tile(ll)
            radout%flux_up_tile(ll, k, i_channel) &
              = weight_channel_incr*(planck_flux_tile(ll, k) &
              + rho_alb(ll, ip_surf_alb_diff, k) &
              *(flux_down_incr(l) &
              + planck_flux_air(l)-planck_flux_tile(ll, k)))
          END DO
        END DO

      END IF

    ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

      cmessage = '*** Error: Tiled surfaces have not yet been ' // &
        'implemented with the spherical harmonic solver.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)

    END IF

    l_initial_channel_tile(i_channel) = .FALSE.

  ELSE
!   Most commonly, this routine will be called to increment
!   rather than to initialize fluxes.
    IF ( (control%i_angular_integration == ip_two_stream).OR. &
         (control%i_angular_integration == ip_ir_gauss) ) THEN

!     Increment the actual fluxes.
      IF (control%isolir == ip_solar) THEN

        IF (control%l_spherical_solar) THEN
          DO k=1, n_tile
            DO ll=1, n_point_tile
              l=list_tile(ll)
              radout%flux_up_tile(ll, k, i_channel) &
                = radout%flux_up_tile(ll, k, i_channel) &
                + weight_channel_incr &
                *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                + rho_alb(ll, ip_surf_alb_dir, k)*flux_direct_incr(l))
            END DO
          END DO  
          IF (control%l_blue_flux_surf.AND. .NOT.control%l_map_sub_bands) THEN
            DO k=1, n_tile
              DO ll=1, n_point_tile
                l=list_tile(ll)
                radout%flux_up_blue_tile(ll, k, i_channel) &
                  = radout%flux_up_blue_tile(ll, k, i_channel) &
                  + weight_blue_incr &
                  *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                  + rho_alb(ll, ip_surf_alb_dir, k)*flux_direct_incr(l))
              END DO
            END DO
          END IF
        ELSE
          DO k=1, n_tile
            DO ll=1, n_point_tile
              l=list_tile(ll)
              radout%flux_up_tile(ll, k, i_channel) &
                = radout%flux_up_tile(ll, k, i_channel) &
                + weight_channel_incr &
                *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                +(rho_alb(ll, ip_surf_alb_dir, k) &
                - rho_alb(ll, ip_surf_alb_diff, k)) &
                * flux_direct_incr(l))
            END DO
          END DO
          IF (control%l_blue_flux_surf .AND. .NOT.control%l_map_sub_bands) THEN
            DO k=1, n_tile
              DO ll=1, n_point_tile
                l=list_tile(ll)
                radout%flux_up_blue_tile(ll, k, i_channel) &
                  = radout%flux_up_blue_tile(ll, k, i_channel) &
                  + weight_blue_incr &
                  *(rho_alb(ll, ip_surf_alb_diff, k)*flux_down_incr(l) &
                  +(rho_alb(ll, ip_surf_alb_dir, k) &
                  - rho_alb(ll, ip_surf_alb_diff, k)) &
                  * flux_direct_incr(l))
              END DO
            END DO
          END IF
        END IF

      ELSE IF (control%isolir == ip_infra_red) THEN

        DO k=1, n_tile
          DO ll=1, n_point_tile
            l=list_tile(ll)
            radout%flux_up_tile(ll, k, i_channel) &
              = radout%flux_up_tile(ll, k, i_channel) &
              + weight_channel_incr*(planck_flux_tile(ll, k) &
              + rho_alb(ll, ip_surf_alb_diff, k) &
              *(flux_down_incr(l) &
              + planck_flux_air(l)-planck_flux_tile(ll, k)))
          END DO
        END DO
      END IF

    ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

      cmessage = '*** Error: Tiled surfaces have not yet been ' // &
        'implemented with the spherical harmonic solver.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)

    END IF
  END IF

END SUBROUTINE augment_channel_tile

END SUBROUTINE augment_tiled_radiance
END MODULE augment_tiled_radiance_mod

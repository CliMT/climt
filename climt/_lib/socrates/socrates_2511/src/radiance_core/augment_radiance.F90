! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to increment radiances or fluxes.
!
! Method:
!   The arrays holding the summed fluxes or radiances are
!   incremented by a weighted sum of the variables suffixed
!   with _INCR. Arguments specify which arrays are to be
!   incremented.
!
!- ---------------------------------------------------------------------
MODULE augment_radiance_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'AUGMENT_RADIANCE_MOD'
CONTAINS
SUBROUTINE augment_radiance(control, sp, atm, bound, radout             &
    , i_band, iex, iex_minor                                            &
    , n_profile, n_layer, n_viewing_level, n_direction                  &
    , l_clear, l_initial, l_initial_band, l_initial_channel             &
    , weight_incr, weight_blue_incr, weight_sub_band_incr               &
!                 Actual radiances
    , i_direct                                                          &
!                 Increments to radiances
    , flux_direct_incr, flux_total_incr, actinic_flux_incr              &
    , i_direct_incr, radiance_incr, photolysis_incr                     &
    , flux_direct_incr_clear, flux_total_incr_clear                     &
    , actinic_flux_incr_clear, k_abs_layer                              &
    , photol, sph, contrib_funci_incr, contrib_funcf_incr               &
!                 Dimensions
    , nd_profile, nd_flux_profile, nd_radiance_profile, nd_j_profile    &
    , nd_layer, nd_viewing_level, nd_direction, nd_channel              &
    , nd_abs, nd_esft_term                                              &
    )


  USE realtype_rd, ONLY: RealK
  USE augment_channel_mod, ONLY: augment_channel
  USE calc_photolysis_incr_mod, ONLY: calc_photolysis_incr
  USE def_control, ONLY: StrCtrl
  USE def_spectrum, ONLY: StrSpecData
  USE def_atm, ONLY: StrAtm
  USE def_bound, ONLY: StrBound
  USE def_out, ONLY: StrOut
  USE def_qy, ONLY: StrQy
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE finalise_photol_incr_mod, ONLY:  finalise_photol_incr
  USE rad_pcf, ONLY: ip_solar, ip_spherical_harmonic, ip_two_stream,    &
                     ip_ir_gauss, ip_sph_mode_flux, ip_sph_mode_rad
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl), INTENT(IN)     :: control

! Spectral data:
  TYPE(StrSpecData), INTENT(IN) :: sp

! Atmospheric properties:
  TYPE(StrAtm), INTENT(IN)      :: atm

! Boundary conditions:
  TYPE(StrBound), INTENT(IN)    :: bound

! Output fields:
  TYPE(StrOut), INTENT(INOUT)   :: radout

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_flux_profile                                                   &
!       Size allocated for points where fluxes are calculated
    , nd_radiance_profile                                               &
!       Size allocated for points where radiances are calculated
    , nd_j_profile                                                      &
!       Size allocated for points where photolysis is calculated
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_channel                                                        &
!       Size allocated for output channels
    , nd_abs                                                            &
!       Size allocated for absorbers
    , nd_esft_term
!       Size allocated for ESFT terms

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      i_band, iex
!       Band and k-term being considered
  INTEGER, INTENT(IN) ::                                                &
      iex_minor(nd_abs)
!       k-term of minor gases (only used with exact_major overlap)
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , n_direction
!       Number of viewing directions
  LOGICAL, INTENT(IN) ::                                                &
      l_clear
!       Clear fluxes calculated
  LOGICAL, INTENT(INOUT) ::                                             &
      l_initial                                                         &
!       Perform initialisation rather than incrementing total fluxes
    , l_initial_band(sp%dim%nd_band)                              &
!       Perform initialisation rather than incrementing band fluxes
    , l_initial_channel(nd_channel)
!       Perform initialisation rather than incrementing channel fluxes

  REAL (RealK), INTENT(IN) ::                                           &
      weight_incr, weight_blue_incr, weight_sub_band_incr
!       Weights to apply to incrementing fluxes

!                 Increments to Fluxes
  REAL (RealK), INTENT(IN) ::                                           &
      flux_direct_incr(nd_flux_profile, 0: nd_layer)                    &
!       Increment to direct flux
    , flux_total_incr(nd_flux_profile, 2*nd_layer+2)                    &
!       Increment to total flux
    , actinic_flux_incr(nd_flux_profile, nd_layer)                      &
!       Increment to actinic flux
    , flux_direct_incr_clear(nd_flux_profile, 0: nd_layer)              &
!       Increment to clear direct flux
    , flux_total_incr_clear(nd_flux_profile, 2*nd_layer+2)              &
!       Increment to clear total flux
    , actinic_flux_incr_clear(nd_flux_profile, nd_layer)
!       Increment to clear actinic flux

  TYPE(StrQy), INTENT(IN) :: photol(sp%photol%n_pathway)
!   Photolysis quantum yields interpolated to model grid temperatures

  TYPE(StrSphGeo), INTENT(IN) :: sph
!   Spherical geometry fields

  REAL (RealK), INTENT(IN) :: &
    k_abs_layer(nd_profile, nd_layer, nd_esft_term, nd_abs)
!       Scaled absorption terms

!                 Increments to Radiances
  REAL (RealK), INTENT(IN) ::                                           &
      i_direct_incr(nd_radiance_profile, 0: nd_layer)                   &
!       Increments to the solar irradiance
    , radiance_incr(nd_radiance_profile, nd_viewing_level               &
        , nd_direction)
!       Increments to the radiance
!                 Increments to Rates of photolysis
  REAL (RealK), INTENT(IN) ::                                           &
      photolysis_incr(nd_j_profile, nd_viewing_level)
!       Increments to the rates of photolysis
  REAL (RealK), INTENT(IN) ::                                           &
      contrib_funci_incr(nd_flux_profile, nd_layer)
!       Contribution function (intensity)
  REAL (RealK), INTENT(IN) ::                                           &
      contrib_funcf_incr(nd_flux_profile, nd_layer)
!       Contribution function (flux)

!                 Total Radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_radiance_profile, 0: nd_layer)
!       Solar irradiance


! Local arguments.
  INTEGER :: i, l, i_sub_k, i_sub, i_abs, i_k_sub, i_k
  INTEGER :: i_channel, i_path, i_gas
!       Loop variables
  REAL (RealK) :: weight_channel_incr
!       Weight to apply to channel increments
  REAL (RealK) :: weight_band_incr
!       Weight to apply to channel increments mapped from bands
  REAL (RealK) :: photolysis_div_incr(nd_flux_profile, nd_layer,        &
                                      sp%dim%nd_pathway)
!       Flux divergence for photolysis increment for the sub-band
  REAL (RealK) :: photolysis_rate_incr(nd_flux_profile, nd_layer,       &
                                       sp%dim%nd_pathway)
!       Photolysis rate increment for the sub-band
  LOGICAL :: l_path(sp%dim%nd_pathway)
!       Flag to calculate pathway increments for these sub-bands
  LOGICAL :: l_weight_set
!       Flag to indicate a weight has been set for this absorber
  LOGICAL :: l_calc_sub_band
!       Flag to calculate increments for a particular sub-band

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='AUGMENT_RADIANCE'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (control%l_map_sub_bands) THEN
    ! Increment the fluxes with sub-bands mapping to channels
    IF (iex == 0) THEN
      ! The k-term is not set so assume a single sub-band for the band
      i_sub = sp%map%list_sub_band_k(1, 1, i_band)
      i_channel = control%map_channel(i_sub)
      weight_channel_incr = weight_sub_band_incr
      photolysis_rate_incr = 0.0_RealK
      photolysis_div_incr = 0.0_RealK
      CALL augment_channel(                                                    &
        control, sp, bound, radout, l_initial_channel, l_clear,                &
        i_channel, n_profile, n_layer, n_viewing_level, n_direction,           &
        weight_channel_incr, flux_direct_incr, flux_total_incr,                &
        actinic_flux_incr, photolysis_div_incr, photolysis_rate_incr,          &
        radiance_incr, photolysis_incr,                                        &
        flux_direct_incr_clear, flux_total_incr_clear, actinic_flux_incr_clear,&
        sph, contrib_funci_incr, contrib_funcf_incr,                           &
        nd_flux_profile, nd_radiance_profile, nd_j_profile,                    &
        nd_layer, nd_viewing_level, nd_direction, nd_channel)
    ELSE
      ! Increment sub-bands for this k-term
      DO i_sub_k=1, sp%map%n_sub_band_k(iex, i_band)
        i_sub = sp%map%list_sub_band_k(i_sub_k, iex, i_band)
        i_channel = control%map_channel(i_sub)
        weight_channel_incr = weight_sub_band_incr &
          * sp%map%weight_sub_band_k(i_sub_k, iex, i_band)
        l_calc_sub_band = .TRUE.
        DO i_abs=2, sp%gas%n_band_absorb(i_band)
          IF (iex_minor(i_abs) > 0) THEN
            ! For the exact_major overlap method the minor gas k-term may
            ! be specified allowing the use of the sub-band weights for
            ! exact overlap with the major gas.
            i_gas = sp%gas%index_absorb(i_abs, i_band)
            l_weight_set = .FALSE.
            DO i_k_sub=1, sp%map%n_k_sub_band(i_gas, i_sub)
              i_k = sp%map%list_k_sub_band(i_k_sub, i_gas, i_sub)
              IF (i_k == iex_minor(i_abs)) THEN
                weight_channel_incr = weight_channel_incr &
                  * sp%map%weight_k_sub_band(i_k_sub, i_gas, i_sub) &
                  / sp%map%weight_k_major(i_k, i_abs, iex, i_band)
                l_weight_set = .TRUE.
              END IF
            END DO
            ! If the minor gas k-term does not contribute to the sub-band
            ! the weight will not have been set (and should be zero) so
            ! no increments need to be calculated for the sub-band.
            l_calc_sub_band = l_calc_sub_band .AND. l_weight_set
          END IF
        END DO
        IF (l_calc_sub_band) THEN
          photolysis_rate_incr = 0.0_RealK
          photolysis_div_incr = 0.0_RealK
          l_path=.FALSE.
          CALL calc_photolysis_incr(                                           &
            sp, photolysis_rate_incr, nd_profile, nd_flux_profile, nd_layer,   &
            nd_esft_term, nd_abs, weight_channel_incr, i_band, l_path,         &
            iex_minor, i_sub, n_layer, n_profile, k_abs_layer, photol )
          CALL finalise_photol_incr(                                           &
            control, sp, atm, l_path,                                          &
            nd_flux_profile, nd_layer, n_profile, n_layer,                     &
            photolysis_div_incr, photolysis_rate_incr, actinic_flux_incr)
          CALL augment_channel(                                                &
            control, sp, bound, radout, l_initial_channel, l_clear,            &
            i_channel, n_profile, n_layer, n_viewing_level, n_direction,       &
            weight_channel_incr, flux_direct_incr, flux_total_incr,            &
            actinic_flux_incr, photolysis_div_incr, photolysis_rate_incr,      &
            radiance_incr, photolysis_incr, flux_direct_incr_clear,            &
            flux_total_incr_clear, actinic_flux_incr_clear,                    &
            sph, contrib_funci_incr, contrib_funcf_incr,                       &
            nd_flux_profile, nd_radiance_profile, nd_j_profile,                &
            nd_layer, nd_viewing_level, nd_direction, nd_channel)
        END IF
      END DO
    END IF
  ELSE
    ! Photolysis increments are always calculated using sub-bands
    photolysis_rate_incr = 0.0_RealK
    photolysis_div_incr = 0.0_RealK
    IF (iex > 0) THEN
      l_path=.FALSE.
      weight_band_incr = 0.0_RealK
      DO i_sub_k=1, sp%map%n_sub_band_k(iex, i_band)
        i_sub = sp%map%list_sub_band_k(i_sub_k, iex, i_band)
        weight_channel_incr = weight_sub_band_incr &
          * sp%map%weight_sub_band_k(i_sub_k, iex, i_band)
        l_calc_sub_band = .TRUE.
        DO i_abs=2, sp%gas%n_band_absorb(i_band)
          IF (iex_minor(i_abs) > 0) THEN
            i_gas = sp%gas%index_absorb(i_abs, i_band)
            l_weight_set = .FALSE.
            DO i_k_sub=1, sp%map%n_k_sub_band(i_gas, i_sub)
              i_k = sp%map%list_k_sub_band(i_k_sub, i_gas, i_sub)
              IF (i_k == iex_minor(i_abs)) THEN
                weight_channel_incr = weight_channel_incr &
                  * sp%map%weight_k_sub_band(i_k_sub, i_gas, i_sub) &
                  / sp%map%weight_k_major(i_k, i_abs, iex, i_band)
                l_weight_set = .TRUE.
              END IF
            END DO
            l_calc_sub_band = l_calc_sub_band .AND. l_weight_set
          END IF
        END DO
        IF (l_calc_sub_band) THEN
          CALL calc_photolysis_incr(                                           &
            sp, photolysis_rate_incr, nd_profile, nd_flux_profile, nd_layer,   &
            nd_esft_term, nd_abs, weight_channel_incr, i_band, l_path,         &
            iex_minor, i_sub, n_layer, n_profile, k_abs_layer, photol )
          weight_band_incr = weight_band_incr + weight_channel_incr
        END IF
      END DO
      CALL finalise_photol_incr(                                               &
        control, sp, atm, l_path,                                              &
        nd_flux_profile, nd_layer, n_profile, n_layer,                         &
        photolysis_div_incr, photolysis_rate_incr, actinic_flux_incr)

      IF (sp%map%n_sub_band_k(iex, i_band) > 1) THEN
        ! Where there is more than one sub-band in the band the sum
        ! of the sub-band weights are used. This is more accurate than
        ! than weight_incr in the case of exact major overlap.
        weight_channel_incr = weight_band_incr
      ELSE
        weight_channel_incr = weight_incr
      END IF
    ELSE
      weight_channel_incr = weight_incr
    END IF
    ! Increment the fluxes with bands mapping to channels
    i_channel = control%map_channel(i_band)
    CALL augment_channel(                                                      &
      control, sp, bound, radout, l_initial_channel, l_clear,                  &
      i_channel, n_profile, n_layer, n_viewing_level, n_direction,             &
      weight_channel_incr, flux_direct_incr, flux_total_incr,                  &
      actinic_flux_incr, photolysis_div_incr, photolysis_rate_incr,            &
      radiance_incr, photolysis_incr,                                          &
      flux_direct_incr_clear, flux_total_incr_clear, actinic_flux_incr_clear,  &
      sph, contrib_funci_incr, contrib_funcf_incr,                             &
      nd_flux_profile, nd_radiance_profile, nd_j_profile,                      &
      nd_layer, nd_viewing_level, nd_direction, nd_channel)

  END IF

  IF (l_initial) THEN
    ! Initialise diagnostic fields that aren't in channels.
    IF ( (control%i_angular_integration == ip_two_stream).OR. &
         (control%i_angular_integration == ip_ir_gauss).OR. &
       ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
          (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

      IF (control%isolir == ip_solar) THEN
        IF (control%l_spherical_solar) THEN
          IF (control%l_blue_flux_surf) THEN
            DO l=1, n_profile
              radout%flux_direct_blue_surf(l) &
                = weight_blue_incr*sph%allsky%flux_direct(l, n_layer+1)
            END DO
          END IF
        ELSE
          IF (control%l_blue_flux_surf) THEN
            DO l=1, n_profile
              radout%flux_direct_blue_surf(l) &
                = weight_blue_incr*flux_direct_incr(l, n_layer)
            END DO
          END IF
        END IF
        IF (control%l_blue_flux_surf) THEN
          DO l=1, n_profile
            radout%flux_up_blue_surf(l) &
              = weight_blue_incr*flux_total_incr(l, 2*n_layer+1)
            radout%flux_down_blue_surf(l) &
              = weight_blue_incr*flux_total_incr(l, 2*n_layer+2)
          END DO
        END IF
      END IF

    ELSE IF ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
              (control%i_sph_mode == ip_sph_mode_rad) ) THEN

      IF (control%isolir == ip_solar) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            i_direct(l, i)=weight_incr*i_direct_incr(l, i)
          END DO
        END DO
      END IF

    END IF
    l_initial = .FALSE.
  ELSE
    IF ( (control%i_angular_integration == ip_two_stream).OR. &
         (control%i_angular_integration == ip_ir_gauss).OR. &
       ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
         (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

      IF (control%isolir == ip_solar) THEN
        IF (control%l_spherical_solar) THEN
          IF (control%l_blue_flux_surf) THEN
            DO l=1, n_profile
              radout%flux_direct_blue_surf(l) &
                = radout%flux_direct_blue_surf(l) &
                + weight_blue_incr*sph%allsky%flux_direct(l, n_layer+1)
            END DO
          END IF
        ELSE
          IF (control%l_blue_flux_surf) THEN
            DO l=1, n_profile
              radout%flux_direct_blue_surf(l) &
                = radout%flux_direct_blue_surf(l) &
                + weight_blue_incr*flux_direct_incr(l, n_layer)
            END DO
          END IF
        END IF
        IF (control%l_blue_flux_surf) THEN
          DO l=1, n_profile
            radout%flux_up_blue_surf(l) &
              = radout%flux_up_blue_surf(l) &
              + weight_blue_incr*flux_total_incr(l, 2*n_layer+1)
            radout%flux_down_blue_surf(l) &
              = radout%flux_down_blue_surf(l) &
              + weight_blue_incr*flux_total_incr(l, 2*n_layer+2)
          END DO
        END IF
      END IF

    ELSE IF ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
              (control%i_sph_mode == ip_sph_mode_rad) ) THEN

      IF (control%isolir == ip_solar) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            i_direct(l, i)=i_direct(l, i) &
              +weight_incr*i_direct_incr(l, i)
          END DO
        END DO
      END IF

    END IF
  END IF


  IF (l_initial_band(i_band)) THEN

!   Initialise the band-by-band fluxes
    IF (control%l_flux_direct_band) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_direct_band(l, i, i_band) &
            = weight_incr*flux_direct_incr(l, i)
        END DO
      END DO
    END IF
    IF (control%l_flux_direct_div_band .AND. &
        control%l_spherical_solar) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          radout%flux_direct_div_band(l, i, i_band) &
            = weight_incr*sph%allsky%flux_direct_div(l, i)
        END DO
      END DO
    END IF
    IF (control%l_flux_direct_sph_band .AND. &
        control%l_spherical_solar) THEN
      DO i=0, n_layer+1
        DO l=1, n_profile
          radout%flux_direct_sph_band(l, i, i_band) &
            = weight_incr*sph%allsky%flux_direct(l, i)
        END DO
      END DO
    END IF
    IF (control%l_flux_down_band) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_down_band(l, i, i_band) &
            = weight_incr*flux_total_incr(l, 2*i+2)
        END DO
      END DO
    END IF
    IF (control%l_flux_up_band) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_up_band(l, i, i_band) &
            = weight_incr*flux_total_incr(l, 2*i+1)
        END DO
      END DO
    END IF
    IF (control%l_flux_div_band .AND. .NOT.control%l_map_sub_bands) THEN
      IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_band(l, i, i_band) &
              = weight_incr * &
              ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
              + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) &
              + sph%allsky%flux_direct_div(l, i) )
          END DO
        END DO
      ELSE
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_band(l, i, i_band) &
              = weight_incr * &
              ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
              + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) )
          END DO
        END DO
      END IF
      IF (control%isolir == ip_solar .AND. control%l_orog .AND. &
        .NOT. control%l_spherical_solar) THEN
        DO l=1, n_profile
          radout%flux_div_band(l, n_layer, i_band) &
            = radout%flux_div_band(l, n_layer, i_band) &
            + weight_incr*flux_direct_incr(l, n_layer) &
            * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
        END DO
      END IF
      DO i_path=1, sp%photol%n_pathway
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_band(l, i, i_band) &
              = radout%flux_div_band(l, i, i_band) &
              - photolysis_div_incr(l, i, i_path)
          END DO
        END DO
      END DO
    END IF
    IF (control%l_actinic_flux_band) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          radout%actinic_flux_band(l, i, i_band) &
            = weight_incr*actinic_flux_incr(l, i)
        END DO
      END DO
    END IF

    IF (l_clear) THEN

      IF (control%l_flux_direct_clear_band .OR. &
           (.NOT.control%l_spherical_solar .AND. &
             ( control%l_cloud_extinction .OR. &
               control%l_ls_cloud_extinction .OR. &
               control%l_cnv_cloud_extinction ) ) ) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_direct_clear_band(l, i, i_band) &
              = weight_incr*flux_direct_incr_clear(l, i)
          END DO
        END DO
      END IF
      IF (control%l_flux_direct_clear_div_band .AND. &
          control%l_spherical_solar) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_direct_clear_div_band(l, i, i_band) &
              = weight_incr*sph%clear%flux_direct_div(l, i)
          END DO
        END DO
      END IF
      IF (control%l_spherical_solar .AND. &
           (control%l_flux_direct_clear_sph_band .OR. &
            control%l_cloud_extinction .OR. &
            control%l_ls_cloud_extinction .OR. &
            control%l_cnv_cloud_extinction)) THEN
        DO i=0, n_layer+1
          DO l=1, n_profile
            radout%flux_direct_clear_sph_band(l, i, i_band) &
              = weight_incr*sph%clear%flux_direct(l, i)
          END DO
        END DO
      END IF
      IF (control%l_flux_down_clear_band) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_down_clear_band(l, i, i_band) &
              = weight_incr*flux_total_incr_clear(l, 2*i+2)
          END DO
        END DO
      END IF
      IF (control%l_flux_up_clear_band .OR. &
          control%l_cloud_absorptivity .OR. &
          control%l_ls_cloud_absorptivity .OR. &
          control%l_cnv_cloud_absorptivity) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_up_clear_band(l, i, i_band) &
              = weight_incr*flux_total_incr_clear(l, 2*i+1)
          END DO
        END DO
      END IF
      IF (control%l_flux_div_clear_band .AND. .NOT.control%l_map_sub_bands) THEN
        IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear_band(l, i, i_band) &
                = weight_incr * &
                ( flux_total_incr_clear(l, 2*i+1) &
                - flux_total_incr_clear(l, 2*i-1) &
                + flux_total_incr_clear(l, 2*i) &
                - flux_total_incr_clear(l, 2*i+2) &
                + sph%clear%flux_direct_div(l, i) )
            END DO
          END DO
        ELSE
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear_band(l, i, i_band) &
                = weight_incr * &
                ( flux_total_incr_clear(l, 2*i+1) &
                - flux_total_incr_clear(l, 2*i-1) &
                + flux_total_incr_clear(l, 2*i) &
                - flux_total_incr_clear(l, 2*i+2) )
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar .AND. control%l_orog .AND. &
          .NOT. control%l_spherical_solar) THEN
          DO l=1, n_profile
            radout%flux_div_clear_band(l, n_layer, i_band) &
              = radout%flux_div_clear_band(l, n_layer, i_band) &
              + weight_incr*flux_direct_incr_clear(l, n_layer) &
              * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
          END DO
        END IF
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear_band(l, i, i_band) &
                = radout%flux_div_clear_band(l, i, i_band) &
                - photolysis_div_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF
      IF (control%l_actinic_flux_clear_band) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%actinic_flux_clear_band(l, i, i_band) &
              = weight_incr*actinic_flux_incr_clear(l, i)
          END DO
        END DO
      END IF

    ELSE ! .NOT. l_clear_band

      IF (control%l_flux_direct_clear_band .OR. &
           (.NOT.control%l_spherical_solar .AND. &
             ( control%l_cloud_extinction .OR. &
               control%l_ls_cloud_extinction .OR. &
               control%l_cnv_cloud_extinction ) ) ) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_direct_clear_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF
      IF (control%l_flux_direct_clear_div_band .AND. &
          control%l_spherical_solar) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_direct_clear_div_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF
      IF (control%l_spherical_solar .AND. &
           (control%l_flux_direct_clear_sph_band .OR. &
            control%l_cloud_extinction .OR. &
            control%l_ls_cloud_extinction .OR. &
            control%l_cnv_cloud_extinction)) THEN
        DO i=0, n_layer+1
          DO l=1, n_profile
            radout%flux_direct_clear_sph_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF
      IF (control%l_flux_down_clear_band) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_down_clear_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF
      IF (control%l_flux_up_clear_band .OR. &
          control%l_cloud_absorptivity .OR. &
          control%l_ls_cloud_absorptivity .OR. &
          control%l_cnv_cloud_absorptivity) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_up_clear_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF
      IF (control%l_flux_div_clear_band) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_clear_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF
      IF (control%l_actinic_flux_clear_band) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%actinic_flux_clear_band(l, i, i_band) = 0.0_RealK
          END DO
        END DO
      END IF      

    END IF

    IF (control%l_contrib_func_band) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          radout%contrib_funci_band(l, i, i_band) &
            = weight_incr*contrib_funci_incr(l, i)
          radout%contrib_funcf_band(l, i, i_band) &
            = weight_incr*contrib_funcf_incr(l, i)
        END DO
      END DO
    END IF

    l_initial_band(i_band) = .FALSE.

  ELSE

!   Increment the band-by-band fluxes
    IF (control%l_flux_direct_band) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_direct_band(l, i, i_band) &
            = radout%flux_direct_band(l, i, i_band) &
            + weight_incr*flux_direct_incr(l, i)
        END DO
      END DO
    END IF
    IF (control%l_flux_direct_div_band .AND. &
        control%l_spherical_solar) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          radout%flux_direct_div_band(l, i, i_band) &
            = radout%flux_direct_div_band(l, i, i_band) &
            + weight_incr*sph%allsky%flux_direct_div(l, i)
        END DO
      END DO
    END IF
    IF (control%l_flux_direct_sph_band .AND. &
        control%l_spherical_solar) THEN
      DO i=0, n_layer+1
        DO l=1, n_profile
          radout%flux_direct_sph_band(l, i, i_band) &
            = radout%flux_direct_sph_band(l, i, i_band) &
            + weight_incr*sph%allsky%flux_direct(l, i)
        END DO
      END DO
    END IF
    IF (control%l_flux_down_band) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_down_band(l, i, i_band) &
            = radout%flux_down_band(l, i, i_band) &
            + weight_incr*flux_total_incr(l, 2*i+2)
        END DO
      END DO
    END IF
    IF (control%l_flux_up_band) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_up_band(l, i, i_band) &
            = radout%flux_up_band(l, i, i_band) &
            + weight_incr*flux_total_incr(l, 2*i+1)
        END DO
      END DO
    END IF
    IF (control%l_flux_div_band .AND. .NOT.control%l_map_sub_bands) THEN
      IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_band(l, i, i_band) &
              = radout%flux_div_band(l, i, i_band) + weight_incr * &
              ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
              + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) &
              + sph%allsky%flux_direct_div(l, i) )
          END DO
        END DO
      ELSE
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_band(l, i, i_band) &
              = radout%flux_div_band(l, i, i_band) + weight_incr * &
              ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
              + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) )
          END DO
        END DO
      END IF
      IF (control%isolir == ip_solar .AND. control%l_orog .AND. &
        .NOT. control%l_spherical_solar) THEN
        DO l=1, n_profile
          radout%flux_div_band(l, n_layer, i_band) &
            = radout%flux_div_band(l, n_layer, i_band) &
            + weight_incr*flux_direct_incr(l, n_layer) &
            * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
        END DO
      END IF
      DO i_path=1, sp%photol%n_pathway
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_div_band(l, i, i_band) &
              = radout%flux_div_band(l, i, i_band) &
              - photolysis_div_incr(l, i, i_path)
          END DO
        END DO
      END DO
    END IF
    IF (control%l_actinic_flux_band) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          radout%actinic_flux_band(l, i, i_band) &
            = radout%actinic_flux_band(l, i, i_band) &
            + weight_incr*actinic_flux_incr(l, i)
        END DO
      END DO
    END IF

    IF (l_clear) THEN
      IF (control%l_flux_direct_clear_band .OR. &
           (.NOT.control%l_spherical_solar .AND. &
             ( control%l_cloud_extinction .OR. &
               control%l_ls_cloud_extinction .OR. &
               control%l_cnv_cloud_extinction ) ) ) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_direct_clear_band(l, i, i_band) &
              = radout%flux_direct_clear_band(l, i, i_band) &
              + weight_incr*flux_direct_incr_clear(l, i)
          END DO
        END DO
      END IF
      IF (control%l_flux_direct_clear_div_band .AND. &
          control%l_spherical_solar) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%flux_direct_clear_div_band(l, i, i_band) &
              = radout%flux_direct_clear_div_band(l, i, i_band) &
              + weight_incr*sph%clear%flux_direct_div(l, i)
          END DO
        END DO
      END IF
      IF (control%l_spherical_solar .AND. &
           (control%l_flux_direct_clear_sph_band .OR. &
            control%l_cloud_extinction .OR. &
            control%l_ls_cloud_extinction .OR. &
            control%l_cnv_cloud_extinction)) THEN
        DO i=0, n_layer+1
          DO l=1, n_profile
            radout%flux_direct_clear_sph_band(l, i, i_band) &
              = radout%flux_direct_clear_sph_band(l, i, i_band) &
              + weight_incr*sph%clear%flux_direct(l, i)
          END DO
        END DO
      END IF
      IF (control%l_flux_down_clear_band) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_down_clear_band(l, i, i_band) &
              = radout%flux_down_clear_band(l, i, i_band) &
              + weight_incr*flux_total_incr_clear(l, 2*i+2)
          END DO
        END DO
      END IF
      IF (control%l_flux_up_clear_band .OR. &
          control%l_cloud_absorptivity .OR. &
          control%l_ls_cloud_absorptivity .OR. &
          control%l_cnv_cloud_absorptivity) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_up_clear_band(l, i, i_band) &
              = radout%flux_up_clear_band(l, i, i_band) &
              + weight_incr*flux_total_incr_clear(l, 2*i+1)
          END DO
        END DO
      END IF
      IF (control%l_flux_div_clear_band .AND. .NOT.control%l_map_sub_bands) THEN
        IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear_band(l, i, i_band) &
                = radout%flux_div_clear_band(l, i, i_band) + weight_incr * &
                ( flux_total_incr_clear(l, 2*i+1) &
                - flux_total_incr_clear(l, 2*i-1) &
                + flux_total_incr_clear(l, 2*i) &
                - flux_total_incr_clear(l, 2*i+2) &
                + sph%clear%flux_direct_div(l, i) )
            END DO
          END DO
        ELSE
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear_band(l, i, i_band) &
                = radout%flux_div_clear_band(l, i, i_band) + weight_incr * &
                ( flux_total_incr_clear(l, 2*i+1) &
                - flux_total_incr_clear(l, 2*i-1) &
                + flux_total_incr_clear(l, 2*i) &
                - flux_total_incr_clear(l, 2*i+2) )
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar .AND. control%l_orog .AND. &
          .NOT. control%l_spherical_solar) THEN
          DO l=1, n_profile
            radout%flux_div_clear_band(l, n_layer, i_band) &
              = radout%flux_div_clear_band(l, n_layer, i_band) &
              + weight_incr*flux_direct_incr_clear(l, n_layer) &
              * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
          END DO
        END IF
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear_band(l, i, i_band) &
                = radout%flux_div_clear_band(l, i, i_band) &
                - photolysis_div_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF
      IF (control%l_actinic_flux_clear_band) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%actinic_flux_clear_band(l, i, i_band) &
              = radout%actinic_flux_clear_band(l, i, i_band) &
              + weight_incr*actinic_flux_incr_clear(l, i)
          END DO
        END DO
      END IF
    END IF

    IF (control%l_contrib_func_band) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          radout%contrib_funci_band(l, i, i_band) &
            = radout%contrib_funci_band(l, i, i_band) &
            + weight_incr*contrib_funci_incr(l, i)
          radout%contrib_funcf_band(l, i, i_band) &
            = radout%contrib_funcf_band(l, i, i_band) &
            + weight_incr*contrib_funcf_incr(l, i)
        END DO
      END DO
    END IF

  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE augment_radiance
END MODULE augment_radiance_mod

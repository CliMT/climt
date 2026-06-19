! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Increment radiances or fluxes for each channel.  
!
! Method:
!   The arrays holding the channel fluxes or radiances are incremented
!   by a weighted sum of the supplied increments. This routine may be
!   called multiple times from augment_radiance when sub-bands are
!   mapped to channels (where each k-term can contribute to multiple
!   sub-bands).
!
!------------------------------------------------------------------------------
MODULE augment_channel_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'AUGMENT_CHANNEL_MOD'

CONTAINS

SUBROUTINE augment_channel( &
! Derived types for the radiance_calc interface
  control, sp, bound, radout, &
! Controlling logicals
  l_initial_channel, l_clear, &
! Channel to be incremented and array sizes for looping
  i_channel, n_profile, n_layer, n_viewing_level, n_direction, &
! Weight for summing increments
  weight_channel_incr, &
! Increments to radiances and fluxes
  flux_direct_incr, flux_total_incr, actinic_flux_incr, &
  photolysis_div_incr, photolysis_rate_incr, radiance_incr, photolysis_incr, &
  flux_direct_incr_clear, flux_total_incr_clear, actinic_flux_incr_clear, &
! Increments to pseudo-spherical fluxes and contribution functions
  sph, contrib_funci_incr, contrib_funcf_incr, &
! Dimensions of arrays
  nd_flux_profile, nd_radiance_profile, nd_j_profile, &
  nd_layer, nd_viewing_level, nd_direction, nd_channel)

  USE realtype_rd, ONLY: RealK
  USE def_bound, ONLY: StrBound
  USE def_control, ONLY: StrCtrl
  USE def_out, ONLY: StrOut
  USE def_spectrum, ONLY: StrSpecData
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_solar, ip_spherical_harmonic, ip_two_stream,    &
                     ip_ir_gauss, ip_sph_mode_flux, ip_sph_mode_rad,    &
                     ip_sph_mode_j

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Control options:
  TYPE(StrCtrl), INTENT(IN)     :: control

! Output fields:
  TYPE(StrOut), INTENT(INOUT)   :: radout

! Spherical geometry fields
  TYPE(StrSphGeo), INTENT(IN)   :: sph

! Boundary conditions:
  TYPE(StrBound), INTENT(IN)    :: bound

! Spectral data:
  TYPE(StrSpecData), INTENT(IN) :: sp

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_flux_profile                                                   &
!       Size allocated for points where fluxes are calculated
    , nd_j_profile                                                      &
!       Size allocated for points where photolysis is calculated
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_radiance_profile                                               &
!       Size allocated for points where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_channel
!       Size allocated for output channels

!                 Increments to Fluxes
  REAL (RealK), INTENT(IN) ::                                           &
      flux_direct_incr(nd_flux_profile, 0: nd_layer)                    &
!       Increment to direct flux
    , photolysis_div_incr(nd_flux_profile, nd_layer, sp%dim%nd_pathway) &
!       Flux divergence for photolysis increment for the sub-band
    , actinic_flux_incr(nd_flux_profile, nd_layer)                      &
!       Increment to actinic flux
    , flux_direct_incr_clear(nd_flux_profile, 0: nd_layer)              &
!       Increment to clear direct flux
    , flux_total_incr_clear(nd_flux_profile, 2*nd_layer+2)              &
!       Increment to clear total flux
    , actinic_flux_incr_clear(nd_flux_profile, nd_layer)
!       Increment to clear actinic flux

!       Flux divergence for photolysis increment for the sub-band
  REAL (RealK), INTENT(IN) ::                                           &
      photolysis_rate_incr(nd_flux_profile, nd_layer,                   &
                                       sp%dim%nd_pathway)               &
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

  REAL (RealK), INTENT(IN) :: weight_channel_incr
!       Weight to apply to channel increments
  REAL (RealK), INTENT(IN) ::                                           &
      flux_total_incr(nd_flux_profile, 2*nd_layer+2)
!       Increment to total flux

  LOGICAL, INTENT(INOUT) :: l_initial_channel(nd_channel)
!       Perform initialisation rather than incrementing channel fluxes

  LOGICAL, INTENT(IN) :: l_clear
!       Clear fluxes calculated

  INTEGER, INTENT(IN) :: i_channel
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , n_direction
!       Number of viewing directions

! Local variables
  INTEGER :: i, k, l, i_path

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='AUGMENT_CHANNEL'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (l_initial_channel(i_channel)) THEN
!   Initialization of the radiance field takes place here.
    IF ( (control%i_angular_integration == ip_two_stream).OR. &
         (control%i_angular_integration == ip_ir_gauss).OR. &
       ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
          (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

      IF (control%l_flux_div) THEN
        IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div(l, i, i_channel) &
                = weight_channel_incr * &
                ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
                + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) &
                + sph%allsky%flux_direct_div(l, i) )
            END DO
          END DO
        ELSE
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div(l, i, i_channel) &
                = weight_channel_incr * &
                ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
                + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) )
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar .AND. control%l_orog .AND. &
          .NOT. control%l_spherical_solar) THEN
          DO l=1, n_profile
            radout%flux_div(l, n_layer, i_channel) &
              = radout%flux_div(l, n_layer, i_channel) &
              + weight_channel_incr*flux_direct_incr(l, n_layer) &
              * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
          END DO
        END IF
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div(l, i, i_channel) &
                = radout%flux_div(l, i, i_channel) &
                - photolysis_div_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF
      IF (control%isolir == ip_solar) THEN
        IF (control%l_spherical_solar) THEN
          DO i=0, n_layer+1
            DO l=1, n_profile
              radout%flux_direct_sph(l, i, i_channel) &
                = weight_channel_incr*sph%allsky%flux_direct(l, i)
            END DO
          END DO
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_direct_div(l, i, i_channel) &
                = weight_channel_incr*sph%allsky%flux_direct_div(l, i)
            END DO
          END DO
        ELSE
          DO i=0, n_layer
            DO l=1, n_profile
              radout%flux_direct(l, i, i_channel) &
                = weight_channel_incr*flux_direct_incr(l, i)
            END DO
          END DO
        END IF
      END IF
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_up(l, i, i_channel) &
            = weight_channel_incr*flux_total_incr(l, 2*i+1)
          radout%flux_down(l, i, i_channel) &
            = weight_channel_incr*flux_total_incr(l, 2*i+2)
        END DO
      END DO
      IF (control%l_actinic_flux) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%actinic_flux(l, i, i_channel) &
              = weight_channel_incr*actinic_flux_incr(l, i)
          END DO
        END DO
      END IF
      IF (control%l_photolysis_rate) THEN
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%photolysis_rate(l, i, i_path, i_channel) &
                = photolysis_rate_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF
      IF (control%l_photolysis_div) THEN
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%photolysis_div(l, i, i_path, i_channel) &
                = photolysis_div_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF

      IF (l_clear) THEN
        IF (control%l_flux_div) THEN
          IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_div_clear(l, i, i_channel) &
                  = weight_channel_incr * &
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
                radout%flux_div_clear(l, i, i_channel) &
                  = weight_channel_incr * &
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
              radout%flux_div_clear(l, n_layer, i_channel) &
                = radout%flux_div_clear(l, n_layer, i_channel) &
                + weight_channel_incr*flux_direct_incr_clear(l, n_layer) &
                * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
            END DO
          END IF
          DO i_path=1, sp%photol%n_pathway
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_div_clear(l, i, i_channel) &
                  = radout%flux_div_clear(l, i, i_channel) &
                  - photolysis_div_incr(l, i, i_path)
              END DO
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar) THEN
          IF (control%l_spherical_solar) THEN
            DO i=0, n_layer+1
              DO l=1, n_profile
                radout%flux_direct_clear_sph(l,i, i_channel) &
                  = weight_channel_incr*sph%clear%flux_direct(l, i)
              END DO
            END DO
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_direct_clear_div(l,i, i_channel) &
                  = weight_channel_incr*sph%clear%flux_direct_div(l, i)
              END DO
            END DO
          ELSE
            DO i=0, n_layer
              DO l=1, n_profile
                radout%flux_direct_clear(l, i, i_channel) &
                  = weight_channel_incr*flux_direct_incr_clear(l, i)
              END DO
            END DO
          END IF
        END IF
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_up_clear(l, i, i_channel) &
              = weight_channel_incr*flux_total_incr_clear(l, 2*i+1)
            radout%flux_down_clear(l, i, i_channel) &
              = weight_channel_incr*flux_total_incr_clear(l, 2*i+2)
          END DO
        END DO
        IF (control%l_actinic_flux_clear) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%actinic_flux_clear(l, i, i_channel) &
                = weight_channel_incr*actinic_flux_incr_clear(l, i)
            END DO
          END DO
        END IF
      ELSE
        IF (control%l_flux_div) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div_clear(l, i, i_channel) = 0.0_RealK
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar) THEN
          IF (control%l_spherical_solar) THEN
            DO i=0, n_layer+1
              DO l=1, n_profile
                radout%flux_direct_clear_sph(l,i, i_channel) = 0.0_RealK
              END DO
            END DO
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_direct_clear_div(l,i, i_channel) = 0.0_RealK
              END DO
            END DO
          ELSE
            DO i=0, n_layer
              DO l=1, n_profile
                radout%flux_direct_clear(l, i, i_channel) = 0.0_RealK
              END DO
            END DO
          END IF
        END IF
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_up_clear(l, i, i_channel) = 0.0_RealK
            radout%flux_down_clear(l, i, i_channel) = 0.0_RealK
          END DO
        END DO
        IF (control%l_actinic_flux_clear) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%actinic_flux_clear(l, i, i_channel) = 0.0_RealK
            END DO
          END DO
        END IF
      END IF

      IF (control%l_contrib_func) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%contrib_funci(l, i, i_channel) &
              = weight_channel_incr*contrib_funci_incr(l, i)
            radout%contrib_funcf(l, i, i_channel) &
              = weight_channel_incr*contrib_funcf_incr(l, i)
          END DO
        END DO
      END IF

    ELSE IF ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
              (control%i_sph_mode == ip_sph_mode_rad) ) THEN

      DO k=1, n_direction
        DO i=1, n_viewing_level
          DO l=1, n_profile
            radout%radiance(l, i, k, i_channel) &
              = weight_channel_incr*radiance_incr(l, i, k)
          END DO
        END DO
      END DO

    ELSE IF ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
              (control%i_sph_mode == ip_sph_mode_j) ) THEN

      DO i=1, n_viewing_level
        DO l=1, n_profile
          radout%photolysis(l, i, i_channel) &
            = weight_channel_incr*photolysis_incr(l, i)
        END DO
      END DO

    END IF

    l_initial_channel(i_channel) = .FALSE.

  ELSE
!   Increment the channel diagnostics.
    IF ( (control%i_angular_integration == ip_two_stream).OR. &
         (control%i_angular_integration == ip_ir_gauss).OR. &
       ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
         (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

      IF (control%l_flux_div) THEN
        IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div(l, i, i_channel) &
                = radout%flux_div(l, i, i_channel) + weight_channel_incr * &
                ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
                + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) &
                + sph%allsky%flux_direct_div(l, i) )
            END DO
          END DO
        ELSE
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div(l, i, i_channel) &
                = radout%flux_div(l, i, i_channel) + weight_channel_incr * &
                ( flux_total_incr(l, 2*i+1) - flux_total_incr(l, 2*i-1) &
                + flux_total_incr(l, 2*i)   - flux_total_incr(l, 2*i+2) )
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar .AND. control%l_orog .AND. &
          .NOT. control%l_spherical_solar) THEN
          DO l=1, n_profile
            radout%flux_div(l, n_layer, i_channel) &
              = radout%flux_div(l, n_layer, i_channel) &
              + weight_channel_incr*flux_direct_incr(l, n_layer) &
              * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
          END DO
        END IF
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_div(l, i, i_channel) &
                = radout%flux_div(l, i, i_channel) &
                - photolysis_div_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF
      IF (control%isolir == ip_solar) THEN
        IF (control%l_spherical_solar) THEN
          DO i=0, n_layer+1
            DO l=1, n_profile
              radout%flux_direct_sph(l, i, i_channel) &
                = radout%flux_direct_sph(l, i, i_channel) &
                + weight_channel_incr*sph%allsky%flux_direct(l, i)
            END DO
          END DO
          DO i=1, n_layer
            DO l=1, n_profile
              radout%flux_direct_div(l, i, i_channel) &
                = radout%flux_direct_div(l, i, i_channel) &
                + weight_channel_incr*sph%allsky%flux_direct_div(l, i)
            END DO
          END DO
        ELSE
          DO i=0, n_layer
            DO l=1, n_profile
              radout%flux_direct(l, i, i_channel) &
                = radout%flux_direct(l, i, i_channel) &
                + weight_channel_incr*flux_direct_incr(l, i)
            END DO
          END DO
        END IF
      END IF
      DO i=0, n_layer
        DO l=1, n_profile
          radout%flux_up(l, i, i_channel) &
            = radout%flux_up(l, i, i_channel) &
            + weight_channel_incr*flux_total_incr(l, 2*i+1)
          radout%flux_down(l, i, i_channel) &
            = radout%flux_down(l, i, i_channel) &
            + weight_channel_incr*flux_total_incr(l, 2*i+2)
        END DO
      END DO
      IF (control%l_actinic_flux) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%actinic_flux(l, i, i_channel) &
              = radout%actinic_flux(l, i, i_channel) &
              + weight_channel_incr*actinic_flux_incr(l, i)
          END DO
        END DO
      END IF
      IF (control%l_photolysis_rate) THEN
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%photolysis_rate(l, i, i_path, i_channel) &
                = radout%photolysis_rate(l, i, i_path, i_channel) &
                + photolysis_rate_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF
      IF (control%l_photolysis_div) THEN
        DO i_path=1, sp%photol%n_pathway
          DO i=1, n_layer
            DO l=1, n_profile
              radout%photolysis_div(l, i, i_path, i_channel) &
                = radout%photolysis_div(l, i, i_path, i_channel) &
                + photolysis_div_incr(l, i, i_path)
            END DO
          END DO
        END DO
      END IF

      IF (l_clear) THEN
        IF (control%l_flux_div) THEN
          IF (control%isolir == ip_solar .AND. control%l_spherical_solar) THEN
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_div_clear(l, i, i_channel) &
                  = radout%flux_div_clear(l, i, i_channel) &
                  + weight_channel_incr * &
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
                radout%flux_div_clear(l, i, i_channel) &
                  = radout%flux_div_clear(l, i, i_channel) &
                  + weight_channel_incr * &
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
              radout%flux_div_clear(l, n_layer, i_channel) &
                = radout%flux_div_clear(l, n_layer, i_channel) &
                + weight_channel_incr*flux_direct_incr_clear(l, n_layer) &
                * (bound%orog_corr(l) - 1.0_RealK)/bound%orog_corr(l)
            END DO
          END IF
          DO i_path=1, sp%photol%n_pathway
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_div_clear(l, i, i_channel) &
                  = radout%flux_div_clear(l, i, i_channel) &
                  - photolysis_div_incr(l, i, i_path)
              END DO
            END DO
          END DO
        END IF
        IF (control%isolir == ip_solar) THEN
          IF (control%l_spherical_solar) THEN
            DO i=0, n_layer+1
              DO l=1, n_profile
                radout%flux_direct_clear_sph(l, i, i_channel) &
                 = radout%flux_direct_clear_sph(l, i, i_channel) &
                 + weight_channel_incr*sph%clear%flux_direct(l, i)
              END DO
            END DO
            DO i=1, n_layer
              DO l=1, n_profile
                radout%flux_direct_clear_div(l, i, i_channel) &
                 = radout%flux_direct_clear_div(l, i, i_channel) &
                 + weight_channel_incr*sph%clear%flux_direct_div(l, i)
              END DO
            END DO
          ELSE
            DO i=0, n_layer
              DO l=1, n_profile
                radout%flux_direct_clear(l, i, i_channel) &
                  = radout%flux_direct_clear(l, i, i_channel) &
                  + weight_channel_incr*flux_direct_incr_clear(l, i)
              END DO
            END DO
          END IF
        END IF
        DO i=0, n_layer
          DO l=1, n_profile
            radout%flux_up_clear(l, i, i_channel) &
              = radout%flux_up_clear(l, i, i_channel) &
              + weight_channel_incr*flux_total_incr_clear(l, 2*i+1)
            radout%flux_down_clear(l, i, i_channel) &
              = radout%flux_down_clear(l, i, i_channel) &
              + weight_channel_incr*flux_total_incr_clear(l, 2*i+2)
          END DO
        END DO
        IF (control%l_actinic_flux_clear) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              radout%actinic_flux_clear(l, i, i_channel) &
                = radout%actinic_flux_clear(l, i, i_channel) &
                + weight_channel_incr*actinic_flux_incr_clear(l, i)
            END DO
          END DO
        END IF
      END IF

      IF (control%l_contrib_func) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            radout%contrib_funci(l, i, i_channel) &
              = radout%contrib_funci(l, i, i_channel) &
              + weight_channel_incr*contrib_funci_incr(l, i)
            radout%contrib_funcf(l, i, i_channel) &
              = radout%contrib_funcf(l, i, i_channel) &
              + weight_channel_incr*contrib_funcf_incr(l, i)
          END DO
        END DO
      END IF

    ELSE IF ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
              (control%i_sph_mode == ip_sph_mode_rad) ) THEN

      DO k=1, n_direction
        DO i=1, n_viewing_level
          DO l=1, n_profile
            radout%radiance(l, i, k, i_channel) &
              = radout%radiance(l, i, k, i_channel) &
              + weight_channel_incr*radiance_incr(l, i, k)
          END DO
        END DO
      END DO

    ELSE IF ( (control%i_angular_integration == ip_spherical_harmonic).AND. &
              (control%i_sph_mode == ip_sph_mode_j) ) THEN

      DO i=1, n_viewing_level
        DO l=1, n_profile
          radout%photolysis(l, i, i_channel) &
            = radout%photolysis(l, i, i_channel) &
            + weight_channel_incr*photolysis_incr(l, i)
        END DO
      END DO

    END IF
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE augment_channel
END MODULE augment_channel_mod

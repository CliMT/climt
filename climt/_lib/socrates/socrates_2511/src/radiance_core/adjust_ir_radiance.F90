! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convert differential IR radiances or fluxes to actual ones.
!
!------------------------------------------------------------------------------
MODULE adjust_ir_radiance_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ADJUST_IR_RADIANCE_MOD'
CONTAINS
SUBROUTINE adjust_ir_radiance(control, spectrum, atm, radout, &
    planck, i_band, l_clear)

  USE realtype_rd,  ONLY: RealK
  USE def_control,  ONLY: StrCtrl
  USE def_spectrum, ONLY: StrSpecData
  USE def_atm,      ONLY: StrAtm
  USE def_out,      ONLY: StrOut
  USE def_planck,   ONLY: StrPlanck
  USE rad_pcf, ONLY: ip_ir_gauss, ip_sph_mode_flux, ip_spherical_harmonic,    &
                     ip_two_stream, ip_sph_mode_rad
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),     INTENT(IN)    :: control

! Spectral data:
  TYPE(StrSpecData), INTENT(IN)    :: spectrum

! Atmospheric properties:
  TYPE(StrAtm),      INTENT(IN)    :: atm

! Output fields:
  TYPE(StrOut),      INTENT(INOUT) :: radout


! Planckian emission fields
  TYPE(StrPlanck),   INTENT(IN)    :: planck

  INTEGER, INTENT(IN) :: i_band
!     Band being considered

  LOGICAL, INTENT(IN) :: l_clear
!     Calculate clear-sky fluxes


! Local arguments
  INTEGER :: i, id, l, iex, i_sub_k, i_sub, i_channel, i_major
!       Loop variables
  REAL (RealK) :: weight_channel_incr
!       Weight to apply to channel increments

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='ADJUST_IR_RADIANCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (control%l_map_sub_bands) THEN
    ! Adjust the fluxes with sub-bands mapping to channels
    IF (spectrum%gas%n_band_absorb(i_band) > 0) THEN
      ! Sub-bands correspond to the major gas
      i_major = spectrum%gas%index_absorb(1, i_band)
      DO iex=1, spectrum%gas%i_band_k(i_band, i_major)
        DO i_sub_k=1, spectrum%map%n_sub_band_k(iex, i_band)
          i_sub = spectrum%map%list_sub_band_k(i_sub_k, iex, i_band)
          i_channel = control%map_channel(i_sub)
          weight_channel_incr = control%weight_band(i_band) &
            * spectrum%map%weight_sub_band_k(i_sub_k, iex, i_band)
          CALL adjust_ir_channel()
        END DO
      END DO
    ELSE
      ! No k-terms so there will be a single sub-band for the band
      i_sub = spectrum%map%list_sub_band_k(1, 1, i_band)
      i_channel = control%map_channel(i_sub)
      weight_channel_incr = control%weight_band(i_band)
      CALL adjust_ir_channel()
    END IF
  ELSE
    ! Adjust the fluxes with bands mapping to channels
    i_channel = control%map_channel(i_band)
    weight_channel_incr = control%weight_band(i_band)
    CALL adjust_ir_channel()
  END IF

  ! Adjust the band-by-band fluxes
  IF ( (control%i_angular_integration == ip_two_stream).OR. &
       (control%i_angular_integration == ip_ir_gauss) ) THEN

    IF (control%l_flux_up_band) THEN
      DO i=0, atm%n_layer
        DO l=1, atm%n_profile
          radout%flux_up_band(l, i, i_band) = &
            radout%flux_up_band(l, i, i_band) &
            + planck%flux(l, i)*control%weight_band(i_band)
        END DO
      END DO
    END IF
    IF (control%l_flux_down_band) THEN
      DO i=0, atm%n_layer
        DO l=1, atm%n_profile
          radout%flux_down_band(l, i, i_band) = &
            radout%flux_down_band(l, i, i_band) &
            + planck%flux(l, i)*control%weight_band(i_band)
        END DO
      END DO
    END IF
    IF (l_clear) THEN
      IF (control%l_flux_up_clear_band) THEN
        DO i=0, atm%n_layer
          DO l=1, atm%n_profile
            radout%flux_up_clear_band(l, i, i_band) = &
              radout%flux_up_clear_band(l, i, i_band) &
              + planck%flux(l, i)*control%weight_band(i_band)
          END DO
        END DO
      END IF
      IF (control%l_flux_down_clear_band) THEN
        DO i=0, atm%n_layer
          DO l=1, atm%n_profile
            radout%flux_down_clear_band(l, i, i_band) = &
              radout%flux_down_clear_band(l, i, i_band) &
              + planck%flux(l, i)*control%weight_band(i_band)
          END DO
        END DO
      END IF
    END IF

  ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

!   Planckian radiances are always used with spherical harmonics,
!   even when calculating fluxes. The number of levels should
!   be set appropriately above.
    IF (control%i_sph_mode == ip_sph_mode_flux) THEN
      IF (control%l_flux_up_band) THEN
        DO i=0, atm%n_layer
          DO l=1, atm%n_profile
            radout%flux_up_band(l, i, i_band) = &
              radout%flux_up_band(l, i, i_band) &
              + pi*planck%radiance(l, i+1)*control%weight_band(i_band)
          END DO
        END DO
      END IF
      IF (control%l_flux_down_band) THEN
        DO i=0, atm%n_layer
          DO l=1, atm%n_profile
            radout%flux_down_band(l, i, i_band) = &
              radout%flux_down_band(l, i, i_band) &
              + pi*planck%radiance(l, i+1)*control%weight_band(i_band)
          END DO
        END DO
      END IF
    END IF

  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


CONTAINS


! Adjust the fluxes or radiances for each channel.
! An internal subroutine is used to allow multiple channels to be updated
! where sub-bands map to channels.
SUBROUTINE adjust_ir_channel()
  IMPLICIT NONE

  IF ( (control%i_angular_integration == ip_two_stream) .OR. &
       (control%i_angular_integration == ip_ir_gauss) ) THEN

    DO i=0, atm%n_layer
      DO l=1, atm%n_profile
        radout%flux_up(l, i, i_channel) = &
          radout%flux_up(l, i, i_channel) &
          + planck%flux(l, i)*weight_channel_incr
        radout%flux_down(l, i, i_channel) = &
          radout%flux_down(l, i, i_channel) &
          + planck%flux(l, i)*weight_channel_incr
      END DO
    END DO
    IF (l_clear) THEN
      DO i=0, atm%n_layer
        DO l=1, atm%n_profile
          radout%flux_up_clear(l, i, i_channel) = &
            radout%flux_up_clear(l, i, i_channel) &
            + planck%flux(l, i)*weight_channel_incr
          radout%flux_down_clear(l, i, i_channel) = &
            radout%flux_down_clear(l, i, i_channel) &
            + planck%flux(l, i)*weight_channel_incr
        END DO
      END DO
    END IF

  ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

!   Planckian radiances are always used with spherical harmonics,
!   even when calculating fluxes. The number of levels should
!   be set appropriately above.
    IF (control%i_sph_mode == ip_sph_mode_flux) THEN
      DO i=0, atm%n_layer
        DO l=1, atm%n_profile
          radout%flux_up(l, i, i_channel) = &
            radout%flux_up(l, i, i_channel) &
            + pi*planck%radiance(l, i+1)*weight_channel_incr
          radout%flux_down(l, i, i_channel) = &
            radout%flux_down(l, i, i_channel) &
            + pi*planck%radiance(l, i+1)*weight_channel_incr
        END DO
      END DO
    ELSE IF (control%i_sph_mode == ip_sph_mode_rad) THEN
      DO id=1, atm%n_direction
        DO i=1, atm%n_viewing_level
          DO l=1, atm%n_profile
            radout%radiance(l, i, id, i_channel) = &
              radout%radiance(l, i, id, i_channel) &
              + planck%radiance(l, i)*weight_channel_incr
          END DO
        END DO
      END DO
    END IF

  END IF

END SUBROUTINE adjust_ir_channel

END SUBROUTINE adjust_ir_radiance
END MODULE adjust_ir_radiance_mod

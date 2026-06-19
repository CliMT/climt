! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculations for the photolysis increments that are required within
! the loop over sub-bands.
!
!------------------------------------------------------------------------------
MODULE calc_photolysis_incr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_PHOTOLYSIS_INCR_MOD'

CONTAINS

SUBROUTINE calc_photolysis_incr(                                      &
  sp, photolysis_rate_incr, nd_profile, nd_flux_profile, nd_layer,    &
  nd_esft_term, nd_abs, weight_channel_incr, i_band, l_path,          &
  iex_minor, i_sub, n_layer, n_profile, k_abs_layer, photol )

  USE realtype_rd, ONLY: RealK
  USE def_spectrum, ONLY: StrSpecData
  USE def_qy, ONLY: StrQy

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Spectral data:
  TYPE(StrSpecData), INTENT(IN) :: sp

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_flux_profile                                                   &
!       Size allocated for points where fluxes are calculated
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_abs                                                            &
!       Size allocated for absorbers
    , nd_esft_term
!       Size allocated for ESFT terms

  REAL (RealK), INTENT(INOUT) ::                                        &
                  photolysis_rate_incr(nd_flux_profile, nd_layer,       &
                                       sp%dim%nd_pathway)
!       Photolysis rate increment for the sub-band

  REAL (RealK), INTENT(IN) :: weight_channel_incr
!       Weight to apply to channel increments
  INTEGER, INTENT(IN) :: i_band
!       Band being considered
  LOGICAL, INTENT(INOUT) :: l_path(sp%dim%nd_pathway)
!       Flag to calculate pathway increments for these sub-bands
  INTEGER, INTENT(IN) ::                                                &
      iex_minor(nd_abs)
!       k-term of minor gases (only used with exact_major overlap)

  INTEGER, INTENT(IN) :: i_sub
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) :: &
    k_abs_layer(nd_profile, nd_layer, nd_esft_term, nd_abs)
!       Scaled absorption terms

  TYPE (StrQy), INTENT(IN) :: photol(sp%photol%n_pathway)
!   Photolysis quantum yields interpolated to model grid temperatures

! Local variables
  INTEGER :: i_wl, i_gas_last, i_path, i_gas, i_abs, i_k_sub
  INTEGER :: i, l, i_k
  REAL (RealK) :: sub_band_work, photol_work(n_profile, n_layer)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_PHOTOLYSIS_INCR'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Photolysis calculations that depend on sub-band
  i_gas_last=0
  sub_band_work = weight_channel_incr * sp%var%wavelength_sub_band(0, i_sub)
  DO i_path=1, sp%photol%n_pathway
    i_wl = sp%photol%qy_sub(i_sub, i_path)
    i_gas = sp%photol%pathway_absorber(i_path)
    IF (i_wl > 0 .AND. sp%map%n_k_sub_band(i_gas, i_sub) > 0) THEN
      l_path(i_path) = .TRUE.        
      IF (i_gas /= i_gas_last) THEN
        ! Do calculations that only depend on the absorber for this pathway
        DO i_abs=1, sp%gas%n_band_absorb(i_band)
          ! Find the absorber index for the k_abs_layer array
          IF (sp%gas%index_absorb(i_abs, i_band) == i_gas) EXIT
        END DO
        IF (iex_minor(i_abs) > 0) THEN
          ! If this is a minor gas for the sub-band and we are using the
          ! random overlap method then we can directly use the absorption
          ! for the minor gas k-term.
          DO i_k_sub=1, sp%map%n_k_sub_band(i_gas, i_sub)
            i_k = sp%map%list_k_sub_band(i_k_sub, i_gas, i_sub)
            IF (i_k == iex_minor(i_abs)) THEN
              DO i=1, n_layer
                DO l=1, n_profile
                  photol_work(l, i) &
                    = k_abs_layer(l, i, i_k, i_abs) * sub_band_work
                END DO
              END DO
              EXIT
            ELSE IF (i_k_sub == sp%map%n_k_sub_band(i_gas, i_sub)) THEN
              photol_work = 0.0_RealK
            END IF
          END DO            
        ELSE IF (sp%map%n_k_sub_band(i_gas, i_sub) > 1) THEN
          ! If there is more than one minor gas k-term in the sub-band
          ! and we are not using random overlap then calculate a
          ! weighted mean absorption for the sub-band
          DO i=1, n_layer
            DO l=1, n_profile
              i_k = sp%map%list_k_sub_band(1, i_gas, i_sub)
              photol_work(l, i) &
                = sp%map%weight_k_sub_band(1, i_gas, i_sub) &
                * k_abs_layer(l, i, i_k, i_abs)
              DO i_k_sub=2, sp%map%n_k_sub_band(i_gas, i_sub)
                i_k = sp%map%list_k_sub_band(i_k_sub, i_gas, i_sub)
                photol_work(l, i) = photol_work(l, i) &
                  + sp%map%weight_k_sub_band(i_k_sub, i_gas, i_sub) &
                  * k_abs_layer(l, i, i_k, i_abs)
              END DO
              photol_work(l, i) = photol_work(l, i) * sub_band_work
            END DO
          END DO
        ELSE
          ! Only one term in the sub-band so use the absorption for that term.
          i_k = sp%map%list_k_sub_band(1, i_gas, i_sub)
          DO i=1, n_layer
            DO l=1, n_profile
              photol_work(l, i) = k_abs_layer(l, i, i_k, i_abs) * sub_band_work
            END DO
          END DO
        END IF
      END IF
      ! Finally do calculations that depend on this pathway
      IF (sp%photol%n_t_lookup_photol(i_path) > 1) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            photolysis_rate_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) + photol_work(l, i) &
              * photol(i_path)%qy(l, i, i_wl)
          END DO
        END DO
      ELSE
        DO i=1, n_layer
          DO l=1, n_profile
            photolysis_rate_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) + photol_work(l, i) &
              * sp%photol%quantum_yield(1, i_wl, i_path)
          END DO
        END DO
      END IF
      i_gas_last=i_gas
    END IF
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_photolysis_incr
END MODULE calc_photolysis_incr_mod

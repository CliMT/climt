! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculations for the photolysis increments that can be done after
! the loop over sub-bands.
!
!------------------------------------------------------------------------------
MODULE finalise_photol_incr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FINALISE_PHOTOL_INCR_MOD'

CONTAINS

SUBROUTINE finalise_photol_incr( &
  control, sp, atm, l_path, nd_flux_profile, nd_layer, n_profile, n_layer, &
  photolysis_div_incr, photolysis_rate_incr, actinic_flux_incr)

  USE realtype_rd, ONLY: RealK
  USE gas_list_pcf, ONLY: molar_weight
  USE rad_ccf, ONLY: h_planck, c_light, n_avogadro
  USE def_control, ONLY: StrCtrl
  USE def_spectrum, ONLY: StrSpecData
  USE def_atm, ONLY: StrAtm

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Control options:
  TYPE(StrCtrl), INTENT(IN)     :: control

! Spectral data:
  TYPE(StrSpecData), INTENT(IN) :: sp

! Atmospheric properties:
  TYPE(StrAtm), INTENT(IN)      :: atm

  LOGICAL, INTENT(IN) :: l_path(sp%dim%nd_pathway)
!       Flag to calculate pathway increments for these sub-bands

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_flux_profile                                                   &
!       Size allocated for points where fluxes are calculated
    , nd_layer
!       Size allocated for layers

  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer   
!       Number of layers

  REAL (RealK), INTENT(INOUT) ::                                        &
                  photolysis_div_incr(nd_flux_profile, nd_layer,        &
                                      sp%dim%nd_pathway)
!       Flux divergence for photolysis increment for the sub-band
  REAL (RealK), INTENT(INOUT) ::                                        &
                  photolysis_rate_incr(nd_flux_profile, nd_layer,       &
                                       sp%dim%nd_pathway)
!       Photolysis rate increment for the sub-band
  REAL (RealK), INTENT(IN) ::                                           &
                  actinic_flux_incr(nd_flux_profile, nd_layer)
!       Increment to actinic flux

! Local variables
  INTEGER :: i_gas_last, i_path, i_gas, i, l
  REAL (RealK) :: photol_work, threshold_wavenumber
  REAL (RealK) :: photol_rate_work(n_profile, n_layer)
  REAL (RealK) :: tol = SQRT(TINY(1.0_RealK))

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='FINALISE_PHOTOL_INCR'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Photolysis calculations that do not depend on sub-band
  i_gas_last=0
  DO i_path=1, sp%photol%n_pathway
    IF (l_path(i_path)) THEN
      i_gas = sp%photol%pathway_absorber(i_path)
      IF (i_gas /= i_gas_last) THEN
        photol_work = molar_weight(sp%gas%type_absorb(i_gas)) &
             * 1.0E-03_RealK / ( n_avogadro * h_planck * c_light )
        IF (.NOT.control%l_photol_only(sp%gas%type_absorb(i_gas))) THEN
          WHERE (atm%gas_mix_ratio(1:n_profile, 1:n_layer, i_gas) > tol)
            photol_rate_work &
              = photol_work / atm%gas_mix_ratio(1:n_profile, 1:n_layer, i_gas)
          ELSEWHERE
            photol_rate_work = 0.0_RealK
          END WHERE
        END IF
        i_gas_last=i_gas
      END IF
      IF (control%l_photol_only(sp%gas%type_absorb(i_gas))) THEN
        ! Absorption of flux is neglected for this gas. Only the photolysis
        ! rates are calculated so the gas mixing ratio is not required.
        DO i=1, n_layer
          DO l=1, n_profile
            ! Photolysis reactions per molecule per second
            photolysis_rate_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) * actinic_flux_incr(l, i) &
              * photol_work
          END DO
        END DO
      ELSE IF (sp%photol%l_thermalise(i_path)) THEN
        ! Energy used for photolysis is considered to be thermalised
        ! immediately and will be included in the radiative heating.
        ! In this case the photolysis_div_incr remains at zero.
        DO i=1, n_layer
          DO l=1, n_profile
            ! Photolysis reactions per molecule per second
            photolysis_rate_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) * actinic_flux_incr(l, i) &
              * photol_rate_work(l, i)
          END DO
        END DO
      ELSE
        threshold_wavenumber = 1.0_RealK &
                             / sp%photol%threshold_wavelength(i_path)
        DO i=1, n_layer
          DO l=1, n_profile
            photolysis_rate_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) * actinic_flux_incr(l, i)
            ! Energy used for photolysis (Wm-2)
            photolysis_div_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) * atm%mass(l, i) &
              * threshold_wavenumber
            ! Photolysis reactions per molecule per second
            photolysis_rate_incr(l, i, i_path) &
              = photolysis_rate_incr(l, i, i_path) * photol_rate_work(l, i)
          END DO
        END DO
      END IF
    END IF
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE finalise_photol_incr
END MODULE finalise_photol_incr_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to check the number of terms in the phase function.
!
! Purpose:
!   This subroutine checks the prescription of the phase function
!   against the specified options to ensure that information is
!   present to define all required moments.
!
!- ---------------------------------------------------------------------
MODULE check_phf_term_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CHECK_PHF_TERM_MOD'
CONTAINS
SUBROUTINE check_phf_term(ierr                                          &
    , l_aerosol, n_aerosol, i_aerosol_parametrization                   &
    , n_aerosol_phf_term                                                &
    , n_phase_term_aerosol_prsc                                         &
    , l_cloud, n_condensed, i_condensed_param, i_phase_cmp              &
    , condensed_n_phf                                                   &
    , n_phase_term_drop_prsc, n_phase_term_ice_prsc                     &
    , n_order_phase, l_henyey_greenstein_pf                             &
    , l_rescale, n_order_forward, l_solar_phf, n_order_phase_solar      &
    , nd_aerosol_species, nd_condensed                                  &
    )


  USE rad_pcf, ONLY: i_err_fatal, i_normal, ip_ackerman_stephens,       &
                     ip_aerosol_param_dry, ip_aerosol_param_moist,      &
                     ip_aerosol_param_phf_dry, ip_drop_pade_2,          &
                     ip_aerosol_param_phf_moist, ip_ice_adt,            &
                     ip_aerosol_unparametrized, ip_ice_adt_10,          &
                     ip_drop_unparametrized, ip_ice_fu_ir,              &
                     ip_ice_fu_phf, ip_ice_fu_solar, ip_phase_ice,      &
                     ip_ice_unparametrized, ip_phase_water,             &
                     ip_ps_size_phf, ip_slingo_schr_ice_phf,            &
                     ip_slingo_schr_phf, ip_slingo_schrecker,           &
                     ip_slingo_schrecker_ice, ip_ice_pade_2_phf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Dummy arguments:
  INTEGER                                                               &
      ierr
!       Error flag

! Dimensions of arrays:
  INTEGER, INTENT(IN) ::                                                &
      nd_aerosol_species                                                &
!       Size allocated for species of aerosols
    , nd_condensed
!       Size allocated for condensed components


! Generic variables:
  LOGICAL, INTENT(IN) ::                                                &
      l_henyey_greenstein_pf                                            &
!       Flag for Henyey-Greenstein phase functions
    , l_rescale                                                         &
!       Flag for rescaling of the phase functions
    , l_solar_phf
!       Flag to use a separate treatment of the solar beam
  INTEGER, INTENT(IN) ::                                                &
      n_order_phase                                                     &
!       Order of terms required in the phase function
    , n_order_forward                                                   &
!       Order of the term in the phase function used for rescaling
    , n_order_phase_solar
!       Order of the phase function used in solar calculations

! Aerosol Fields
  LOGICAL, INTENT(IN) ::                                                &
      l_aerosol
!       Flag to use aerosols
  INTEGER, INTENT(IN) ::                                                &
      n_aerosol                                                         &
!       Number of aerosols
    , i_aerosol_parametrization(nd_aerosol_species)                     &
!       Parametrizations adopted for aerosols
    , n_aerosol_phf_term(nd_aerosol_species)                            &
!       Number of terms in the phase function
    , n_phase_term_aerosol_prsc(nd_aerosol_species)
!       Number of terms in the prescribed phase functions
!       for each species

! Cloudy Fields
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud
!       Flag to include clouds
  INTEGER, INTENT(IN) ::                                                &
      n_condensed                                                       &
!       Number of condensed components
    , i_condensed_param(nd_condensed)                                   &
!       Parametrizations adopted for condensed components
    , i_phase_cmp(nd_condensed)                                         &
!       Phases of the condensed components
    , condensed_n_phf(nd_condensed)                                     &
!       Number of terms in the phase function
    , n_phase_term_drop_prsc                                            &
!       Number of terms in the prescribed phase functions
!       for each droplets
    , n_phase_term_ice_prsc
!       Number of terms in the prescribed phase functions
!       for each ice crystals


! Local variables:
  LOGICAL                                                               &
      l_inadequate
!       Flag for inadequate information
  INTEGER                                                               &
      j                                                                 &
!       Loop variable
    , n_order_required
!       Order of terms which are required in the phase function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CHECK_PHF_TERM'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Determine the order of terms for which information in the
! phase function will be required.
  IF (l_henyey_greenstein_pf) THEN
    n_order_required=1
  ELSE
    IF (l_rescale) THEN
      n_order_required=MAX(n_order_phase, n_order_forward)
    ELSE
      n_order_required=n_order_phase
    END IF
!   If the solar beam is to be treated separately more terms
!   may be required.
    IF (l_solar_phf) THEN
      n_order_required = MAX(n_order_phase, n_order_phase_solar)
      IF (l_rescale) n_order_required=n_order_required+1
    END IF
  END IF

! If aerosols are included carry out the required checks.
  IF (l_aerosol) THEN
    l_inadequate=.FALSE.
    DO j=1, n_aerosol
      IF ( (i_aerosol_parametrization(j) ==                             &
            ip_aerosol_param_dry).OR.                                   &
           (i_aerosol_parametrization(j) ==                             &
            ip_aerosol_param_moist) ) THEN
!       In this case information will be extended as a
!       Henyey-Greenstein phase function; and the available
!       information will include the asymmetry.
        CONTINUE
      ELSE IF ( (i_aerosol_parametrization(j) ==                        &
                 ip_aerosol_param_phf_dry).OR.                          &
                (i_aerosol_parametrization(j) ==                        &
                 ip_aerosol_param_phf_moist) ) THEN
        l_inadequate=(n_order_required >  n_aerosol_phf_term(j))
      ELSE IF (  i_aerosol_parametrization(j) ==                        &
                 ip_aerosol_unparametrized) THEN
        l_inadequate=(n_order_required > n_phase_term_aerosol_prsc(j))
      END IF

      IF (l_inadequate) THEN
        WRITE(cmessage, '(/A, /A, I3, A)')                              &
          '*** Error: There is not enough information to define'        &
          , 'the phase function for aerosol ', j                        &
          , ' to the desired order.'
        ierr=i_err_fatal
        GO TO 9999
      END IF
    END DO
  END IF

  IF (l_cloud) THEN
    l_inadequate=.FALSE.
    DO j=1, n_condensed
      IF (i_phase_cmp(j) == ip_phase_water) THEN
        IF ( (i_condensed_param(j) == ip_slingo_schrecker).OR.          &
             (i_condensed_param(j) == ip_ackerman_stephens).OR.         &
             (i_condensed_param(j) == ip_drop_pade_2) ) THEN
!         The phase function will be extended as a
!         Henyey-Greenstein phase function from information
!         already present.
          CONTINUE
        ELSE IF ( (i_condensed_param(j) == ip_slingo_schr_phf).OR.      &
                  (i_condensed_param(j) == ip_ps_size_phf) ) THEN
          l_inadequate=(n_order_required >  condensed_n_phf(j))
        ELSE IF (  i_condensed_param(j) == ip_drop_unparametrized) THEN
          l_inadequate=(n_order_required > n_phase_term_drop_prsc)
        END IF
        IF (l_inadequate) THEN
          WRITE(cmessage, '(/A, /A, I3, A, /A)')                        &
            '*** Error: There is not enough information to define'      &
            , 'the phase function for condensed species ', j            &
            , ' (water droplets) ', 'to the desired order.'
          ierr=i_err_fatal
          GO TO 9999
        END IF
      ELSE IF (i_phase_cmp(j) == ip_phase_ice) THEN
        IF ( (i_condensed_param(j) == ip_slingo_schrecker_ice).OR.      &
             (i_condensed_param(j) == ip_ice_adt).OR.                   &
             (i_condensed_param(j) == ip_ice_fu_solar).OR.              &
             (i_condensed_param(j) == ip_ice_fu_ir).OR.                 &
             (i_condensed_param(j) == ip_ice_adt_10) ) THEN
!         The phase function will be extended as a
!         Henyey-Greenstein phase function from information
!         already present.
          CONTINUE
        ELSE IF ( (i_condensed_param(j) == ip_slingo_schr_ice_phf).OR.  &
                  (i_condensed_param(j) == ip_ice_fu_phf).OR.           &
                  (i_condensed_param(j) == ip_ice_pade_2_phf) ) THEN
          l_inadequate=(n_order_required >  condensed_n_phf(j))
        ELSE IF (i_condensed_param(j) == ip_ice_unparametrized) THEN
          l_inadequate=(n_order_required > n_phase_term_ice_prsc)
        END IF
        IF (l_inadequate) THEN
          WRITE(cmessage, '(/A, /A, I3, A)')                            &
            '*** Error: There is not enough information to define'      &
            , 'the phase function for condensed species ', j            &
            , ' (ice crystals) to the desired order.'
          ierr=i_err_fatal
          GO TO 9999
        END IF
      END IF

    END DO
  END IF


! Check error condition
 9999 IF (ierr /= i_normal) THEN
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_phf_term
END MODULE check_phf_term_mod

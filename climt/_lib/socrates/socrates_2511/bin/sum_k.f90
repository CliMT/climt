! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to sum two k-distributions.
!
! Method:
!   Calculates the sum of two k-distributions. The k-distributions
!   need to have the same weights, otherwise an error is raised.
!
!------------------------------------------------------------------------------
MODULE sum_k_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SUM_K_MOD'
CONTAINS
SUBROUTINE sum_k(ierr, n_profile, n_layer                                      &
     , n_abs_esft, k_abs_layer, w_abs_esft                                     &
     , n_add_esft, k_add_layer, w_add_esft                                     &
     , nd_profile, nd_layer, nd_k_term, nd_k_term_add)

  USE realtype_rd, ONLY: RealK
  USE rad_pcf,     ONLY: i_err_fatal
  USE yomhook,     ONLY: lhook, dr_hook
  USE parkind1,    ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                       &
      nd_profile                                                               &
!       Size allocated for atmospheric profiles
    , nd_layer                                                                 &
!       Size allocated for atmospheric layers
    , nd_k_term                                                                &
!       Size allocated for combined ESFT terms
    , nd_k_term_add
!       Size allocated for continuum ESFT terms

! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                                    &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                       &
      n_profile                                                                &
!       Number of profiles
    , n_layer
!       Number of layers
  INTEGER, INTENT(IN) ::                                                       &
      n_abs_esft                                                               &
!       Number of ESFT terms for absorbers
    , n_add_esft
!       Number of ESFT terms for k-distribution to be added
  REAL(RealK), INTENT(INOUT) ::                                                &
      k_abs_layer(nd_profile, nd_layer, nd_k_term)
!       Absorber ESFT terms
  REAL(RealK), INTENT(IN) ::                                                   &
      w_abs_esft(nd_k_term)                                                    &
!       Absorber ESFT weights
    , k_add_layer(nd_profile, nd_layer, nd_k_term_add)                         &
!       ESFT terms to be added
    , w_add_esft(nd_k_term_add)
!       ESFT weights to be added

! Local variables.
  REAL(RealK) ::                                                               &
      eps = 1.0E-08_RealK
!       Tolerance used to take into account loss of precision due to
!       number of decimal points on weights in spectral files
  LOGICAL ::                                                                   &
      l_weights_match
!       Flag indicating if the k-term weights match

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SUM_K'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  l_weights_match = .FALSE.
  IF (n_abs_esft == n_add_esft) THEN
    IF (ALL(ABS(w_abs_esft(1:n_abs_esft) - w_add_esft(1:n_add_esft)) <         &
            eps*w_abs_esft(1:n_abs_esft))) THEN
      k_abs_layer(1:n_profile, 1:n_layer, 1:n_abs_esft) =                      &
        k_abs_layer(1:n_profile, 1:n_layer, 1:n_abs_esft) +                    &
        k_add_layer(1:n_profile, 1:n_layer, 1:n_add_esft)
      l_weights_match = .TRUE.
    END IF
  END IF

  IF (.NOT.l_weights_match) THEN
    cmessage = '*** Error: k-distributions do not have identical weights.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sum_k
END MODULE sum_k_mod

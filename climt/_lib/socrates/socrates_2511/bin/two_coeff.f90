! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!   The basic two-stream coefficients in the differential equations
!   are calculated. These are then used to determine the
!   transmission and reflection coefficients. Coefficients for
!   determining the solar or infra-red source terms are calculated.
!
!- ---------------------------------------------------------------------
MODULE two_coeff_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TWO_COEFF_MOD'
CONTAINS
SUBROUTINE two_coeff(ierr, control                                      &
     , n_profile, i_layer_first, i_layer_last                           &
     , i_2stream                                                        &
     , asymmetry, omega, tau_dir, tau                                   &
     , isolir, sec_0, path_div                                          &
     , trans, reflect, trans_0_dir, trans_0                             &
     , source_coeff                                                     &
     , nd_profile                                                       &
     , id_op_lt, id_op_lb, id_trs_lt, id_trs_lb                         &
     , nd_source_coeff                                                  &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_solar
  USE vectlib_mod, ONLY : sqrt_v
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE def_control, ONLY: StrCtrl
  USE solar_coefficient_basic_mod, ONLY: solar_coefficient_basic
  USE trans_source_coeff_mod, ONLY: trans_source_coeff
  USE two_coeff_basic_mod, ONLY: two_coeff_basic

  IMPLICIT NONE

! Control options:
  TYPE(StrCtrl), INTENT(IN) :: control

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , id_op_lt                                                          &
!       Topmost declared layer for optical properties
    , id_op_lb                                                          &
!       Bottom declared layer for optical properties
    , id_trs_lt                                                         &
!       Topmost declared layer for transmission coefficients
    , id_trs_lb                                                         &
!       Bottom declared layer for transmission coefficients
    , nd_source_coeff
!       Size allocated for source coefficients


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to consider
    , i_layer_last                                                      &
!       Last layer to consider
    , isolir                                                            &
!       Spectral region
    , i_2stream
!       Two stream scheme

! Optical properties of layer:
  REAL (RealK), INTENT(IN) ::                                           &
      asymmetry(nd_profile, id_op_lt: id_op_lb)                         &
!       Asymmetry factor
    , omega(nd_profile, id_op_lt: id_op_lb)                             &
!       Albedo of single scattering
    , tau(nd_profile, id_op_lt: id_op_lb)                               &
!       Optical depth
    , tau_dir(nd_profile, id_op_lt: id_op_lb)
!       Unscaled optical depth 

! Solar beam
  REAL (RealK), INTENT(IN) ::                                           &
      sec_0(nd_profile)                                                 &
!       Secant of zenith angle
    , path_div(nd_profile, id_trs_lt: id_trs_lb)
!       Path scaling for spherical geometry

! Coefficients in the two-stream equations:
  REAL (RealK), INTENT(OUT) ::                                          &
      trans(nd_profile, id_trs_lt: id_trs_lb)                           &
!       Diffuse transmission coefficient
    , reflect(nd_profile, id_trs_lt: id_trs_lb)                         &
!       Diffuse reflection coefficient
    , trans_0(nd_profile, id_trs_lt: id_trs_lb)                         &
!       Direct transmission coefficient
    , trans_0_dir(nd_profile, id_trs_lt: id_trs_lb)                     &
!       Direct transmission coefficient without scaling
    , source_coeff(nd_profile, id_trs_lt: id_trs_lb                     &
        , nd_source_coeff)
!       Source coefficients in two-stream equations


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

! Coefficients in the two-stream equations:
  REAL (RealK) ::                                                       &
      lambda(nd_profile, id_op_lt: id_op_lb)                            &
!       Coefficients in two-stream equations
    , sum(nd_profile, id_op_lt: id_op_lb)                               &
!       Sum of alpha_1 and alpha_2
    , diff(nd_profile, id_op_lt: id_op_lb)                              &
!       Difference of alpha_1 and alpha_2
    , gamma_up(nd_profile, id_op_lt: id_op_lb)                          &
!       Basic solar coefficient for upward radiation
    , gamma_down(nd_profile, id_op_lt: id_op_lb)
!       Basic solar coefficient for downward radiation

  REAL (RealK) ::                                                       &
       temp(nd_profile)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='TWO_COEFF'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the basic two-stream coefficients. (The single
! scattering albedo has already been perturbed away from 1 in
! SINGLE_SCATTERING.)
  CALL two_coeff_basic(ierr                                             &
    , n_profile, i_layer_first, i_layer_last                            &
    , i_2stream                                                         &
    , asymmetry, omega                                                  &
    , sum, diff                                                         &
    , nd_profile, id_op_lt, id_op_lb                                    &
    )

! LAMBDA is now calculated.
  DO i=i_layer_first, i_layer_last
     DO l=1, n_profile
        temp(l)=sum(l,i)*diff(l,i)
     END DO
     CALL sqrt_v(n_profile,temp,lambda(1,i))
  END DO


! Calculate the basic coefficients for the solar source terms.
  IF (isolir == ip_solar) THEN
!   LAMBDA may be perturbed by this routine to avoid
!   ill-conditioning for the singular zenith angle.
    CALL solar_coefficient_basic(control                                &
      , n_profile, i_layer_first, i_layer_last                          &
      , omega, asymmetry, sec_0, path_div                               &
      , sum, diff, lambda                                               &
      , gamma_up, gamma_down                                            &
      , nd_profile, id_op_lt, id_op_lb, id_trs_lt, id_trs_lb            &
      )
  END IF


! Determine the transmission and reflection coefficients.
  CALL trans_source_coeff(control                                       &
    , n_profile, i_layer_first, i_layer_last                            &
    , tau_dir, tau, sum, diff, lambda, sec_0, path_div                  &
    , gamma_up, gamma_down                                              &
    , trans, reflect, trans_0_dir, trans_0, source_coeff                &
    , nd_profile                                                        &
    , id_op_lt, id_op_lb, id_trs_lt, id_trs_lb                          &
    , nd_source_coeff                                                   &
    )


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE two_coeff
END MODULE two_coeff_mod

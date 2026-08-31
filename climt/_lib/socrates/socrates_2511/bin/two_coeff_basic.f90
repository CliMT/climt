! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate basic coefficients in two-stream equations.
!
! Method:
!   Depending on the two-stream equations employed, the
!   appropriate coefficients for the fluxes are calculated.
!
!- ---------------------------------------------------------------------
MODULE two_coeff_basic_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TWO_COEFF_BASIC_MOD'
CONTAINS
SUBROUTINE two_coeff_basic(ierr                                         &
     , n_profile, i_layer_first, i_layer_last                           &
     , i_2stream                                                        &
     , asymmetry, omega                                                 &
     , sum, diff                                                        &
     , nd_profile, id_lt, id_lb                                         &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_2s_test, ip_coakley_chylek_1,      &
                     ip_coakley_chylek_2, ip_discrete_ord, ip_eddington,&
                     ip_elsasser, ip_hemi_mean, ip_ifm, ip_krschg_flux, &
                     ip_meador_weaver, ip_pifm80, ip_pifm85, ip_zdk_flux
  USE diffusivity_factor, ONLY: elsasser_factor
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , id_lt                                                             &
!       Topmost declared layer
    , id_lb
!       Bottom declared layer

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
    , i_2stream
!       Two stream scheme

! Optical properties of layer:
  REAL (RealK), INTENT(IN) ::                                           &
      asymmetry(nd_profile, id_lt: id_lb)                               &
!       Asymmetry factor
    , omega(nd_profile, id_lt: id_lb)
!       Albedo of single scattering


! coefficients in the two-stream equations:
  REAL (RealK), INTENT(OUT) ::                                          &
      sum(nd_profile, id_lt: id_lb)                                     &
!       Sum of alpha_1 and alpha_2
    , diff(nd_profile, id_lt: id_lb)
!       Difference of alpha_1 and alpha_2

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK), PARAMETER ::                                            &
      root_3 = 1.7320508075688772e+00_RealK
!       Square root of 3

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'TWO_COEFF_BASIC'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (i_2stream == ip_eddington) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=1.5e+00_RealK*(1.0e+00_RealK                          &
          -omega(l, i)*asymmetry(l, i))
        diff(l, i)=2.0e+00_RealK*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_elsasser) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=elsasser_factor                                       &
          -1.5e+00_RealK*omega(l, i)*asymmetry(l, i)
        diff(l, i)=elsasser_factor*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_discrete_ord) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=root_3*(1.0e+00_RealK                                 &
          -omega(l, i)*asymmetry(l, i))
        diff(l, i)=root_3*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_pifm85) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=2.0e+00_RealK                                         &
          -1.5e+00_RealK*omega(l, i)*asymmetry(l, i)
        diff(l, i)=2.0e+00_RealK*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_2s_test) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=1.5e+00_RealK                                         &
          -1.5e+00_RealK*omega(l, i)*asymmetry(l, i)
        diff(l, i)=1.5e+00_RealK*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_hemi_mean) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=2.0e+00_RealK                                         &
          *(1.0e+00_RealK-omega(l, i)*asymmetry(l, i))
        diff(l, i)=2.0e+00_RealK*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_pifm80) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum(l, i)=2.0e+00_RealK                                         &
          -1.5e+00_RealK*omega(l, i)*asymmetry(l, i)                    &
          -0.5e+00_RealK*omega(l, i)
        diff(l, i)=2.0e+00_RealK*(1.0e+00_RealK-omega(l, i))
      END DO
    END DO

  ELSE IF (i_2stream == ip_ifm) THEN
    cmessage = '*** Error: The improved flux mathod has '               &
      //'not been implemented.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE IF (i_2stream == ip_zdk_flux) THEN
    cmessage = '*** Error: Zdunkowski''s flux method has '              &
      //'not been implemented.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE IF (i_2stream == ip_krschg_flux) THEN
    cmessage = '*** Error: Kerschgen''s flux method has '               &
      //'not been implemented.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE IF (i_2stream == ip_coakley_chylek_1) THEN
    cmessage = '*** Error: Coakley and Chylek''s first method has '     &
      //'not been implemented.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE IF (i_2stream == ip_coakley_chylek_2) THEN
    cmessage = '*** Error: Coakley and Chylek''s second method has '    &
      //'not been implemented.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE IF (i_2stream == ip_meador_weaver) THEN
    cmessage = '*** Error: Meador and Weaver''s method has '            &
      //'not been implemented.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE
    cmessage = '*** Error: An unrecognized value has been specified '   &
      //'to define the two-stream scheme.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE two_coeff_basic
END MODULE two_coeff_basic_mod

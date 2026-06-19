! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to gather data for spline fitting.
!
! Method:
! Splined fits to the given data at the corresponding pressure
! levels are carried out. Optical properties at the required
! pressure levels are calculated from the splined fits.
! This routine is not intended to run in vector mode,
! and has therefore not been optimized for such use.
!
!- ---------------------------------------------------------------------
MODULE prsc_gather_spline_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE prsc_gather_spline(ierr                                      &
    , n_profile, i_first_layer, i_last_layer, p_eval                    &
    , n_opt_level_prsc, prsc_pressure, prsc_opt_property                &
    , opt_property                                                      &
    , nd_profile, nd_layer, id_1, id_2                                  &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    )

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_range, i_normal
  USE spline_evaluate_mod, ONLY: spline_evaluate
  USE spline_fit_mod, ONLY: spline_fit

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , id_1                                                              &
!       Topmost declared layer
    , id_2                                                              &
!       Bottom declared layer for optical properties
    , nd_profile_prsc                                                   &
!       Size allowed for profiles of prescribed properties
    , nd_opt_level_prsc
!       Size allowed for levels of prescribed properties

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_first_layer                                                     &
!       First layer in which the field is required
    , i_last_layer
!       Last layer in which the field is required

  INTEGER, INTENT(IN) ::                                                &
      n_opt_level_prsc
!       Number of levels of prescribed optical data
  REAL (RealK), INTENT(IN) ::                                           &
      prsc_pressure(nd_profile_prsc, nd_opt_level_prsc)                 &
!       Pressures of specified levels
    , prsc_opt_property(nd_profile_prsc, nd_opt_level_prsc)             &
!       Prescribed optical properties
    , p_eval(nd_profile, nd_layer)
!       Pressures where the property is to be evaluated

  REAL (RealK), INTENT(OUT) ::                                          &
      opt_property(nd_profile, id_1: id_2)
!       Calculated optical property

! Local variables.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , i
!       Loop variable
  REAL (RealK) ::                                                       &
      prsc_pressure_g(nd_opt_level_prsc)                                &
!       Pressures of specified levels
!       gathered to the current profile
    , prsc_opt_property_g(nd_opt_level_prsc)                            &
!       Prescribed optical property
!       gathered to the current profile
    , d2(nd_opt_level_prsc)
!       Second derivatives for spline fits



! Calculate the second derivatives for the spline fits.
  DO l=1, n_profile

!   Because of checking for data which are out of range the
!   splining routines do not work in vector mode, so points
!   are gathered to a single profile.

    DO i=1, n_opt_level_prsc
      prsc_pressure_g(i)=prsc_pressure(l, i)
      prsc_opt_property_g(i)=prsc_opt_property(l, i)
    END DO

!   Calculate second derivatives for fitting.
    CALL spline_fit(n_opt_level_prsc, prsc_pressure_g                   &
      , prsc_opt_property_g, d2)

    DO i=i_first_layer, i_last_layer
      CALL spline_evaluate(ierr, n_opt_level_prsc                       &
        , prsc_pressure_g, prsc_opt_property_g                          &
        , d2, p_eval(l, i)                                              &
        , opt_property(l, i)                                            &
        )
!     Here, values which are out of range are silently set to 0.
      IF (ierr /= i_normal) THEN
        IF (ierr == i_err_range) THEN
          opt_property(l, i)=0.0e+00_RealK
!         Recover from this error.
          ierr=i_normal
        ELSE
          RETURN
        END IF
      END IF

    END DO

  END DO

END SUBROUTINE prsc_gather_spline
END MODULE prsc_gather_spline_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the layers in which radiances are required.
!
! Purpose:
!   This determines the layers of the atmosphere where the analytic
!   expression for the radiance must be intercepted to give values
!   on the correct levels.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
MODULE set_rad_layer_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_RAD_LAYER_MOD'
CONTAINS
SUBROUTINE set_rad_layer(ierr                                           &
    , n_layer, n_viewing_level, viewing_level                           &
    , i_rad_layer, frac_rad_layer                                       &
    , nd_viewing_level                                                  &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


  INTEGER, INTENT(IN) ::                                                &
      nd_viewing_level
!       Size allocated for levels where radiances are calculated

! Dummy arguments
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

  INTEGER, INTENT(IN) ::                                                &
      n_viewing_level                                                   &
!       Number of levels on which to calculate the radiance
    , n_layer
!       Number of atmospheric layers
  REAL (RealK), INTENT(IN) ::                                           &
      viewing_level(nd_viewing_level)
!       Levels where radiances are calculated
  INTEGER, INTENT(OUT) ::                                               &
      i_rad_layer(nd_viewing_level)
!       Layers in which to intercept radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers where radiances
!       are calculated


! Local Variables
  INTEGER                                                               &
      i
!       Loop variable
  REAL (RealK) ::                                                       &
      tol_bdy
!       The tolerance detecting the closeness of boundaries

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SET_RAD_LAYER'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerance for detecting boundaries.
  tol_bdy=1.6e+01_RealK*EPSILON(tol_bdy)

  DO i=1, n_viewing_level

!   Check that a level is not above the top of the atmosphere.
    IF (viewing_level(i) <  0.0e+00_RealK) THEN
      cmessage = '*** Error: A viewing level is above the TOA.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    i_rad_layer(i)=INT(viewing_level(i))+1
    frac_rad_layer(i)=1.0e+00_RealK+viewing_level(i)                    &
      -REAL(i_rad_layer(i), RealK)

!   At the bottom of the atmosphere this will give a value greater
!   than N_LAYER, so we reset, but check that an unreasonable
!   attempt to get radiances below the column has not been made:
!   this will give a fatal error.
    IF (i_rad_layer(i) >  n_layer) THEN
      IF (frac_rad_layer(i) <  tol_bdy) THEN
        i_rad_layer(i)=i_rad_layer(i)-1
        frac_rad_layer(i)=1.0e+00_RealK
      ELSE
        cmessage = '*** Error: A viewing level is below the surface.'
        ierr=i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    END IF

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_rad_layer
END MODULE set_rad_layer_mod

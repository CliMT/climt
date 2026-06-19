! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate solar scattering angles.
!
! Purpose:
!   This routine returns the cosines of the angles of scattering
!   from the solar beam for each viewing direction.
!
! Method:
!   A scalar product of the solar and viewing directions is
!   evaluated. This routine is called only when radiances are to
!   be calculated, so ND_PROFILE can be used for all horizontal
!   dimensions.
!
!- ---------------------------------------------------------------------
MODULE sol_scat_cos_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOL_SCAT_COS_MOD'
CONTAINS
SUBROUTINE sol_scat_cos(n_profile, n_direction                          &
    , mu_0, direction, cos_sol_view                                     &
    , nd_profile, nd_direction)


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_direction
!       Size allocated for viewing directions

  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric profiles
    , n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      mu_0(nd_profile)                                                  &
!       Cosines of solar zenith angles
    , direction(nd_profile, nd_direction, 2)
!       Viewing directions stored as the cosine of the polar
!       viewing angle and the azimuthal viewing angle itself
!       realative to the solar direction
  REAL (RealK), INTENT(OUT) ::                                          &
      cos_sol_view(nd_profile, nd_direction)
!       Cosines of the angles between the solar beam and the
!       viewing directions


! Local variables
  INTEGER                                                               &
      id                                                                &
!       Loop variable
    , l
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOL_SCAT_COS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO id=1, n_direction
    DO l=1, n_profile
      cos_sol_view(l, id)=-mu_0(l)*direction(l, id, 1)                  &
        +SQRT((1.0e+00_RealK-mu_0(l)*mu_0(l))                           &
        *(1.0e+00_RealK-direction(l, id, 1)*direction(l, id, 1)))       &
        *COS(direction(l, id, 2))
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sol_scat_cos
END MODULE sol_scat_cos_mod

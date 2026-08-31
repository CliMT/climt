! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to increment radiances for radiation from the top.
!
! Purpose:
!   The contribution to the solution of radiances transmitted from
!   the top boundary is evaluated. In the IR where differential
!   radiances are used the radiance at the top will be the Planckian
!   radiance at that temperature. In idealized tests an incident
!   flux may be prescribed.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
MODULE calc_top_rad_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_TOP_RAD_MOD'
CONTAINS
SUBROUTINE calc_top_rad(n_profile, tau                                  &
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
    , n_direction, mu_v                                                 &
    , flux_inc_down                                                     &
    , radiance                                                          &
    , nd_profile, nd_layer, nd_direction, nd_viewing_level              &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_direction
!       Size allocated for viewing directions


! The atmosphere:
  INTEGER, INTENT(IN) ::                                                &
      n_profile
!       Number of atmospheric profiles
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)
!       Optical depths

! Viewing geometry:
  INTEGER, INTENT(IN) ::                                                &
      n_direction                                                       &
!       Number of directions
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , i_rad_layer(nd_viewing_level)
!       Indices of layers containing viewing levels
  REAL (RealK), INTENT(IN) ::                                           &
      frac_rad_layer(nd_viewing_level)                                  &
!       Fraction optical depth into its layer of the
!       viewing level
    , mu_v(nd_profile, nd_direction)
!       Cosines of polar viewing angles
  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_down(nd_profile)
!       Isotropic incident flux

  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_profile, nd_viewing_level, nd_direction)
!       Radiances (to be incremented by the contribution of
!       the particular integral)

! Local variables
  INTEGER                                                               &
      l                                                                 &
!       Loop variable (points)
    , iv                                                                &
!       Loop variable (viewing levels)
    , id                                                                &
!       Loop variable (directions)
    , i                                                                 &
!       Loop variable
    , ll
!       Loop variable
  INTEGER                                                               &
      n_list_down                                                       &
!       Number of points where the viewing direction is upward
    , list_down(nd_profile)
!       List of points where the viewing direction is upward
  REAL (RealK) ::                                                       &
      tau_c(nd_profile, nd_viewing_level)
!       Cumulative optical depths to the viewing level

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_TOP_RAD'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the cumulative optical depth from the
! top of the atmosphere to each viewing level.
  DO iv=1, n_viewing_level

    DO l=1, n_profile
      tau_c(l, iv)                                                      &
        =frac_rad_layer(iv)*tau(l, i_rad_layer(iv))
    END DO
    DO i=i_rad_layer(iv)-1, 1, -1
      DO l=1, n_profile
        tau_c(l, iv)=tau_c(l, iv)+tau(l, i)
      END DO
    END DO

  END DO


  DO id=1, n_direction

!   Collect downward directions.
    n_list_down=0
    DO l=1, n_profile
      IF (mu_v(l, id) <  0.0e+00_RealK) THEN
        n_list_down=n_list_down+1
        list_down(n_list_down)=l
      END IF
    END DO

    DO iv=1, n_viewing_level
      DO ll=1, n_list_down
        l=list_down(ll)
        radiance(l, iv, id)=radiance(l, iv, id)                         &
          +(flux_inc_down(l)/pi)*EXP(tau_c(l, iv)/mu_v(l, id))
      END DO
    END DO

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_top_rad
END MODULE calc_top_rad_mod

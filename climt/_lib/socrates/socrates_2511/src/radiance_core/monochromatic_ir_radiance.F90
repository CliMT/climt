! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the infra-red radiance ignoring scattering.
!
! Method:
!   Using the secant of the ray transmission coefficients for
!   each layer may be defined and source terms may be calculated.
!   The upward and downward radiances are integrated along
!   their paths.
!
!- ---------------------------------------------------------------------
MODULE monochromatic_ir_radiance_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MONOCHROMATIC_IR_RADIANCE_MOD'
CONTAINS
SUBROUTINE monochromatic_ir_radiance(n_profile, n_layer                 &
     , tau                                                              &
     , rad_inc_down                                                     &
     , diff_planck, source_ground, albedo_surface_diff                  &
     , secant_ray                                                       &
     , radiance                                                         &
     , nd_profile, nd_layer                                             &
     )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer
!       Maximum number of layers

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)                                         &
!       Optical depths of layers
    , rad_inc_down(nd_profile)                                          &
!       Incident downward radiance
    , source_ground(nd_profile)                                         &
!       Source function of ground
    , albedo_surface_diff(nd_profile)                                   &
!       Diffuse albedo
    , diff_planck(nd_profile, nd_layer)                                 &
!       Difference in Planckian function
    , secant_ray
!       Secant of angle with vertical
  REAL (RealK), INTENT(OUT) ::                                          &
      radiance(nd_profile, 2*nd_layer+2)
!       Diffuse radiance

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      trans(nd_profile, nd_layer)                                       &
!       Transmissivities
    , source_up(nd_profile, nd_layer)                                   &
!       Upward source function
    , source_down(nd_profile, nd_layer)
!       Downward source function

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r                                                             &
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision
    , sq_eps_r
!       The square root of the above

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MONOCHROMATIC_IR_RADIANCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(tau(1, 1))
  sq_eps_r=SQRT(eps_r)

  DO i=1, n_layer
    DO l=1, n_profile
      trans(l, i)=EXP(-secant_ray*tau(l, i))
    END DO
  END DO

  DO i=1, n_layer
    DO l=1, n_profile
      source_up(l, i)=(1.0e+00_RealK-trans(l, i)+sq_eps_r)              &
        *diff_planck(l, i)                                              &
        /(secant_ray*tau(l, i)+sq_eps_r)
      source_down(l, i)=-source_up(l, i)
    END DO
  END DO

! Downward radiance.
  DO l=1, n_profile
    radiance(l, 2)=rad_inc_down(l)
  END DO
  DO i=1, n_layer
    DO l=1, n_profile
      radiance(l, 2*i+2)=trans(l, i)*radiance(l, 2*i)                   &
        +source_down(l, i)
    END DO
  END DO

! Upward radiance.
  DO l=1, n_profile
    radiance(l, 2*n_layer+1)=source_ground(l)                           &
      +albedo_surface_diff(l)*radiance(l, 2*n_layer+2)
  END DO
  DO i=n_layer, 1, -1
    DO l=1, n_profile
      radiance(l, 2*i-1)=trans(l, i)*radiance(l, 2*i+1)                 &
        +source_up(l, i)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE monochromatic_ir_radiance
END MODULE monochromatic_ir_radiance_mod

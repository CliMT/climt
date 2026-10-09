! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE spherical_solar_source_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SPHERICAL_SOLAR_SOURCE_MOD'
CONTAINS

! Subroutine to calculate the solar flux and source terms
! using spherical geometry.
!
!------------------------------------------------------------------------------
SUBROUTINE spherical_solar_source(bound, sph_common, sph_comp,                 &
    n_profile, n_layer,                                                        &
    trans_0, source_coeff, l_scale_solar, adjust_solar_ke, s_down, s_up,       &
    nd_profile, nd_layer, nd_source_coeff)


  USE realtype_rd,            ONLY: RealK
  USE def_bound,              ONLY: StrBound
  USE def_spherical_geometry, ONLY: StrSphCommon, StrSphComp
  USE rad_pcf,                ONLY: ip_scf_solar_up, ip_scf_solar_down
  USE parkind1,               ONLY: jprb, jpim
  USE yomhook,                ONLY: lhook, dr_hook

  IMPLICIT NONE

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Spherical geometry fields common to all-sky / clear-sky
  TYPE(StrSphCommon), INTENT(IN)    :: sph_common

! Spherical geometry fields for component (either all-sky or clear-sky)
  TYPE(StrSphComp),   INTENT(INOUT) :: sph_comp

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_source_coeff
!       Size allocated for coefficients in the source terms


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers

  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Scaling applied to solar beam

  REAL (RealK), INTENT(IN) ::                                           &
      trans_0(nd_profile, nd_layer)                                     &
!       Direct transmission coefficient
    , source_coeff(nd_profile, nd_layer, nd_source_coeff)               &
!       Reflection coefficient
    , adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment to solar flux


  REAL (RealK), INTENT(OUT) ::                                          &
      s_down(nd_profile, nd_layer)                                      &
!       Downward source function
    , s_up(nd_profile, nd_layer)
!       Upward source function


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) :: f_forward, f_backward
!       Working variables for calculating up and down scattered flux
  
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SPHERICAL_SOLAR_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (l_scale_solar) THEN
!   The solar flux may be multiplied by a scaling factor if an
!   equivalent extinction is used.
    DO i=0, n_layer+1
      DO l=1, n_profile
        sph_comp%flux_direct(l, i) = sph_common%flux_inc_direct(l,i) &
          *sph_comp%trans_0(l,i)*sph_common%adjust_solar_ke(l, i)
      END DO
    END DO
    DO i=1, n_layer
      DO l=1, n_profile
        ! The flux scattered vertically up and down the column is
        ! calculated assuming the scattered flux in the forward and
        ! backwards directions along the slant path are isotropic 
        ! across the respective hemispheres.
        f_forward  = 0.5_RealK * &
          ( 1.0_RealK + bound%cos_zen(l, i)*sph_common%path_div(l, i) )
        f_backward = 0.5_RealK * &
          ( 1.0_RealK - bound%cos_zen(l, i)*sph_common%path_div(l, i) )
        s_up(l, i) = &
          ( f_forward  * source_coeff(l, i, ip_scf_solar_up) &
          + f_backward * source_coeff(l, i, ip_scf_solar_down) ) &
          * sph_comp%flux_direct(l, i)
        s_down(l, i) = &
          ( f_forward  * source_coeff(l, i, ip_scf_solar_down) &
          + f_backward * source_coeff(l, i, ip_scf_solar_up) ) &
          *sph_comp%flux_direct(l, i)
        sph_comp%flux_direct_div(l, i)=sph_comp%flux_direct(l, i) &
          *(1.0_RealK-trans_0(l, i)*adjust_solar_ke(l, i))
        ! Convert flux to Watts per square metre normal to the beam
        ! for diagnostic output
        sph_comp%flux_direct(l, i) &
          = sph_comp%flux_direct(l, i) * sph_common%path_div(l, i)
      END DO
    END DO
  ELSE
    DO i=0, n_layer+1
      DO l=1, n_profile
        sph_comp%flux_direct(l, i) = sph_common%flux_inc_direct(l,i) &
          *sph_comp%trans_0(l,i)
      END DO
    END DO
    DO i=1, n_layer
      DO l=1, n_profile
        f_forward  = 0.5_RealK * &
          ( 1.0_RealK + bound%cos_zen(l, i)*sph_common%path_div(l, i) )
        f_backward = 0.5_RealK * &
          ( 1.0_RealK - bound%cos_zen(l, i)*sph_common%path_div(l, i) )
        s_up(l, i) = &
          ( f_forward  * source_coeff(l, i, ip_scf_solar_up) &
          + f_backward * source_coeff(l, i, ip_scf_solar_down) ) &
          * sph_comp%flux_direct(l, i)
        s_down(l, i) = &
          ( f_forward  * source_coeff(l, i, ip_scf_solar_down) &
          + f_backward * source_coeff(l, i, ip_scf_solar_up) ) &
          *sph_comp%flux_direct(l, i)
        sph_comp%flux_direct_div(l, i)=sph_comp%flux_direct(l, i) &
          *(1.0_RealK-trans_0(l, i))
        sph_comp%flux_direct(l, i) &
          = sph_comp%flux_direct(l, i) * sph_common%path_div(l, i)
      END DO
    END DO
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE spherical_solar_source
END MODULE spherical_solar_source_mod

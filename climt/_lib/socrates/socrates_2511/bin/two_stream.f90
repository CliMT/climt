! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve the two-stream equations in a column.
!
! Method:
!   The coefficients of the two-stream equations are calculated.
!   From these we obtain the transmission and reflection
!   coefficients and the source terms. Depending on the solver
!   selected, an appropriate set of matrix equations is formulated
!   and solved to give the fluxes.
!
!- ---------------------------------------------------------------------
MODULE two_stream_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TWO_STREAM_MOD'
CONTAINS
SUBROUTINE two_stream(ierr                                              &
    , control, bound                                                    &
!                 Atmospheric Properties
    , n_profile, n_layer                                                &
!                 Two-stream Scheme
    , i_2stream                                                         &
!                 Options for Solver
    , i_solver, i_scatter_method                                        &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , diff_planck                                                       &
    , l_ir_source_quad, diff_planck_2                                   &
!                 Conditions at TOA
    , flux_inc_down, flux_inc_direct, sec_0                             &
!                 Surface Conditions
    , diffuse_albedo, direct_albedo, d_planck_flux_surface              &
!                 Spherical geometry
    , sph                                                               &
!                 Single Scattering Properties
    , tau_dir, tau, omega, asymmetry                                    &
!                 Fluxes Calculated
    , flux_direct, flux_total                                           &
!                 Dimensions
    , nd_profile, nd_layer, nd_source_coeff                             &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_infra_red, ip_no_scatter_abs, ip_no_scatter_ext, &
                     ip_scatter_approx, ip_scatter_full
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE column_solver_mod, ONLY: column_solver
  USE ir_source_mod, ONLY: ir_source
  USE two_coeff_fast_lw_mod, ONLY: two_coeff_fast_lw
  USE two_coeff_mod, ONLY: two_coeff

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_source_coeff
!       Size allocated for source coefficients

! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , isolir                                                            &
!       Spectral region
    , i_solver                                                          &
!       Solver employed
    , i_2stream                                                         &
!       Two-stream scheme
    , i_scatter_method
!       Method of treating scattering
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar                                                     &
!       Scaling applied to solar flux
    , l_ir_source_quad
!       Use quadratic source term
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)                                         &
!       Optical depth
    , tau_dir(nd_profile, nd_layer)                                     &
!       Optical depth for direct flux
    , omega(nd_profile, nd_layer)                                       &
!       Albedo of single scattering
    , asymmetry(nd_profile, nd_layer)                                   &
!       Asymmetry
    , sec_0(nd_profile)                                                 &
!       Secants of solar zenith angles
    , diffuse_albedo(nd_profile)                                        &
!       Diffuse albedo
    , direct_albedo(nd_profile)                                         &
!       Direct albedo
    , flux_inc_down(nd_profile)                                         &
!       Incident total flux
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , diff_planck(nd_profile, nd_layer)                                 &
!       Difference in Planckian fluxes across layers
    , d_planck_flux_surface(nd_profile)                                 &
!       Ground source function
    , adjust_solar_ke(nd_profile, nd_layer)                             &
!       Adjustment of solar beam with equivalent extinction
    , diff_planck_2(nd_profile, nd_layer)
!       2x2nd differences of Planckian
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_total(nd_profile, 2*nd_layer+2)
!       Total fluxes
  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!       Spherical geometry fields


! Local variables.
  REAL (RealK) ::                                                       &
      trans(nd_profile, nd_layer)                                       &
!       Transmission of layer
    , reflect(nd_profile, nd_layer)                                     &
!       Reflectance of layer
    , trans_0(nd_profile, nd_layer)                                     &
!       Direct transmittance
    , trans_0_dir(nd_profile, nd_layer)                                 &
!       Direct transmittance without scaling
    , source_coeff(nd_profile, nd_layer, nd_source_coeff)               &
!       Source coefficients
    , s_down(nd_profile, nd_layer)                                      &
!       Downward source
    , s_up(nd_profile, nd_layer)
!       Upward source

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'TWO_STREAM'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the two-stream coefficients.
  IF ( (i_scatter_method == ip_scatter_full) .OR.                       &
       (i_scatter_method == ip_scatter_approx) ) THEN
    CALL two_coeff(ierr, control                                        &
      , n_profile, 1, n_layer                                           &
      , i_2stream                                                       &
      , asymmetry, omega, tau_dir, tau                                  &
      , isolir, sec_0, sph%common%path_div                              &
      , trans, reflect, trans_0_dir, trans_0                            &
      , source_coeff                                                    &
      , nd_profile, 1, nd_layer, 1, nd_layer, nd_source_coeff           &
      )
  ELSE IF ( (i_scatter_method == ip_no_scatter_abs) .OR.                &
            (i_scatter_method == ip_no_scatter_ext) ) THEN
    CALL two_coeff_fast_lw(n_profile, 1, n_layer                        &
      , l_ir_source_quad, tau                                           &
      , trans, source_coeff                                             &
      , nd_profile, nd_layer, 1, nd_layer, nd_source_coeff)
    reflect(1:n_profile, 1:n_layer)=0.0e+00_RealK
  END IF

! Calculate the appropriate source terms.
  IF (isolir == ip_infra_red) THEN
    CALL ir_source(n_profile, 1, n_layer                                &
      , source_coeff, diff_planck                                       &
      , l_ir_source_quad, diff_planck_2                                 &
      , s_down, s_up                                                    &
      , nd_profile, nd_layer, nd_source_coeff                           &
      )
  END IF

  CALL column_solver(ierr, control, bound, sph%common, sph%allsky       &
    , n_profile, n_layer                                                &
    , i_scatter_method, i_solver                                        &
    , trans, reflect, trans_0_dir, trans_0, source_coeff                &
    , isolir, flux_inc_direct, flux_inc_down                            &
    , s_down, s_up                                                      &
    , diffuse_albedo, direct_albedo                                     &
    , d_planck_flux_surface                                             &
    , l_scale_solar, adjust_solar_ke                                    &
    , flux_direct, flux_total                                           &
    , nd_profile, nd_layer, nd_source_coeff)


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE two_stream
END MODULE two_stream_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate fluxes in a homogeneous column.
!
! Method:
!   The optical properties of the column are passed to a solver
!   appropriate for the treatment of scattering. Only homogeneous
!   columns are treated here so this may be used to calculate
!   clear-sky fluxes or overcast cloudy fluxes (such as sub-columns
!   with MCICA).
!
!- ---------------------------------------------------------------------
MODULE column_solver_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'COLUMN_SOLVER_MOD'
CONTAINS
SUBROUTINE column_solver(ierr, control, bound, sph_common, sph_comp     &
    , n_profile, n_layer                                                &
    , i_scatter_method, i_solver_clear                                  &
    , trans, reflect, trans_0_dir, trans_0                              &
    , source_coeff                                                      &
    , isolir, flux_inc_direct, flux_inc_down                            &
    , s_down, s_up                                                      &
    , albedo_surface_diff, albedo_surface_dir                           &
    , d_planck_flux_surface                                             &
    , l_scale_solar, adjust_solar_ke                                    &
    , flux_direct, flux_total                                           &
    , nd_profile, nd_layer, nd_source_coeff                             &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE def_spherical_geometry, ONLY: StrSphCommon, StrSphComp
  USE rad_pcf, ONLY: i_err_fatal, ip_solar, ip_no_scatter_abs,          &
                     ip_no_scatter_ext, ip_solver_pentadiagonal,        &
                     ip_solver_homogen_direct
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE band_solver_mod, ONLY: band_solver
  USE spherical_solar_source_mod, ONLY: spherical_solar_source
  USE set_matrix_pentadiagonal_mod, ONLY: set_matrix_pentadiagonal
  USE solar_source_mod, ONLY: solar_source
  USE solver_homogen_direct_mod, ONLY: solver_homogen_direct
  USE solver_no_scat_mod, ONLY: solver_no_scat

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

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
!       Size allocated for layers

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , isolir                                                            &
!       Spectral region
    , i_solver_clear                                                    &
!       Solver for clear fluxes
    , i_scatter_method
!       Method of treating scattering
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Scaling applied to solar beam
  REAL (RealK), INTENT(IN) ::                                           &
      trans(nd_profile, nd_layer)                                       &
!       Transmission coefficients
    , reflect(nd_profile, nd_layer)                                     &
!       Reflection coefficients
    , trans_0(nd_profile, nd_layer)                                     &
!       Direct transmission coefficients
    , trans_0_dir(nd_profile, nd_layer)                                 &
!       Direct transmission with direct tau
    , source_coeff(nd_profile, nd_layer, nd_source_coeff)               &
!       Coefficients in source terms
    , albedo_surface_diff(nd_profile)                                   &
!       Diffuse albedo
    , albedo_surface_dir(nd_profile)                                    &
!       Direct albedo
    , flux_inc_down(nd_profile)                                         &
!       Incident total flux
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , d_planck_flux_surface(nd_profile)                                 &
!       Difference between the Planckian flux at the surface
!       temperature and that of the overlaying air
    , adjust_solar_ke(nd_profile, nd_layer)
!       Scaling of solar beam

  REAL (RealK), INTENT(INOUT) ::                                        &
      s_down(nd_profile, nd_layer)                                      &
!       Downward source
    , s_up(nd_profile, nd_layer)
!       Upward source

  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_total(nd_profile, 2*nd_layer+2)
!       Total fluxes (diffuse only when using spherical geometry)


! Local variabales.
  INTEGER ::                                                            &
      l                                                                 &
!       Loop variable
    , n_equation
!       Number of equations
  REAL (RealK) ::                                                       &
      source_ground(nd_profile)                                         &
!       Source from ground
    , a5(nd_profile, 5, 2*nd_layer+2)                                   &
!       Pentadiagonal matrix
    , b(nd_profile, 2*nd_layer+2)                                       &
!       Rhs of matrix equation
    , work_1(nd_profile, 2*nd_layer+2)
!       Working array for solver

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'COLUMN_SOLVER'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The source functions only need to be recalculated in the visible.
  IF (isolir == ip_solar) THEN
    IF (control%l_spherical_solar) THEN
      CALL spherical_solar_source(bound, sph_common, sph_comp,          &
        n_profile, n_layer,                                             &
        trans_0, source_coeff,                                          &
        l_scale_solar, adjust_solar_ke,                                 &
        s_down, s_up,                                                   &
        nd_profile, nd_layer, nd_source_coeff)
      DO l=1, n_profile
        source_ground(l)=albedo_surface_dir(l)                          &
          *sph_comp%flux_direct(l, n_layer+1)
      END DO
    ELSE
      CALL solar_source(control, bound, n_profile, n_layer,             &
        flux_inc_direct, trans_0_dir, trans_0, source_coeff,            &
        l_scale_solar, adjust_solar_ke,                                 &
        flux_direct, s_down, s_up,                                      &
        nd_profile, nd_layer, nd_source_coeff)
      DO l=1, n_profile
        source_ground(l)=(albedo_surface_dir(l)-albedo_surface_diff(l)) &
          *flux_direct(l, n_layer)
      END DO
    END IF
  ELSE
    DO l=1, n_profile
      source_ground(l) = (1.0_RealK-albedo_surface_diff(l))             &
        *d_planck_flux_surface(l)
    END DO    
  END IF


! Select an appropriate solver for the equations of transfer.
  IF ( (i_scatter_method == ip_no_scatter_abs) .OR.                     &
       (i_scatter_method == ip_no_scatter_ext) ) THEN

!   Solve for the fluxes ignoring scattering
    CALL solver_no_scat(n_profile, n_layer                              &
      , trans, s_down, s_up                                             &
      , albedo_surface_diff, flux_inc_down, d_planck_flux_surface       &
      , flux_total                                                      &
      , nd_profile, nd_layer)

  ELSE IF (i_solver_clear == ip_solver_pentadiagonal) THEN

!   Calculate the elements of the matrix equations.
    CALL set_matrix_pentadiagonal(n_profile, n_layer                    &
      , trans, reflect                                                  &
      , s_down, s_up                                                    &
      , albedo_surface_diff                                             &
      , flux_inc_down                                                   &
      , source_ground                                                   &
      , a5, b                                                           &
      , nd_profile, nd_layer                                            &
      )
    n_equation=2*n_layer+2

    CALL band_solver(n_profile, n_equation                              &
      , 2, 2                                                            &
      , a5, b                                                           &
      , flux_total                                                      &
      , work_1                                                          &
      , nd_profile, 5, 2*nd_layer+2                                     &
      )

  ELSE IF (i_solver_clear == ip_solver_homogen_direct) THEN

!   Solve for the fluxes in the column directly.
    CALL solver_homogen_direct(n_profile, n_layer                       &
      , trans, reflect                                                  &
      , s_down, s_up                                                    &
      , albedo_surface_diff                                             &
      , flux_inc_down                                                   &
      , source_ground                                                   &
      , flux_total                                                      &
      , nd_profile, nd_layer                                            &
      )

  ELSE

     cmessage = '*** Error: The specified solver is not valid.'
     ierr=i_err_fatal
     CALL ereport(RoutineName, ierr, cmessage)

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE column_solver
END MODULE column_solver_mod

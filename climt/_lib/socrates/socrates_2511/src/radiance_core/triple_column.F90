! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve the two-stream equations in a triple column.
!
! Method:
!   The atmospheric column is divided into three regions
!   in each layer and the two-stream coefficients are determined
!   for each region. The equations are then solved using
!   appropriate coupling of the fluxes at the boundaries
!   of layers.
!
!- ---------------------------------------------------------------------
MODULE triple_column_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TRIPLE_COLUMN_MOD'
CONTAINS
SUBROUTINE triple_column(ierr                                           &
     , control, bound                                                   &
!                 Atmospheric properties
     , n_profile, n_layer                                               &
!                 Two-stream scheme
     , i_2stream                                                        &
!                 Options for solver
     , i_solver, i_scatter_method                                       &
!                 Options for equivalent extinction
     , l_scale_solar, adjust_solar_ke                                   &
!                 Spectral region
     , isolir                                                           &
!                 Infra-red properties
     , diff_planck                                                      &
     , l_ir_source_quad, diff_planck_2                                  &
!                 Conditions at TOA
     , flux_inc_down, flux_inc_direct, sec_0                            &
!                 Conditions at surface
     , diffuse_albedo, direct_albedo, d_planck_flux_surface             &
!                 Spherical geometry
     , sph                                                              &
!                 Optical Properties
     , ss_prop                                                          &
!                 Cloud geometry
     , n_cloud_top                                                      &
     , n_cloud_type, frac_cloud                                         &
     , n_region, i_region_cloud, frac_region                            &
     , w_free, w_cloud                                                  &
     , cloud_overlap                                                    &
!                 fluxes calculated
     , flux_direct, flux_total                                          &
!                 flags for clear-sky calculations
     , l_clear, i_solver_clear                                          &
!                 Clear-sky fluxes calculated
     , flux_direct_clear, flux_total_clear                              &
!                 Dimensions of arrays
     , nd_profile, nd_layer, nd_layer_clr, id_ct                        &
     , nd_max_order, nd_source_coeff                                    &
     , nd_cloud_type, nd_region, nd_overlap_coeff                       &
     )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: i_err_fatal, ip_infra_red, ip_solar,               &
                     ip_scatter_full, ip_scatter_approx,                &
                     ip_no_scatter_abs, ip_no_scatter_ext,              &
                     ip_region_clear, ip_region_strat, ip_region_conv,  &
                     ip_solver_triple, ip_solver_triple_hogan,          &
                     ip_solver_triple_app_scat
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE column_solver_mod, ONLY: column_solver
  USE set_n_source_coeff_mod, ONLY: set_n_source_coeff
  USE ir_source_mod, ONLY: ir_source
  USE solver_triple_app_scat_mod, ONLY: solver_triple_app_scat
  USE solver_triple_mod, ONLY: solver_triple
  USE solver_triple_hogan_mod, ONLY: solver_triple_hogan
  USE triple_solar_source_mod, ONLY: triple_solar_source
  USE two_coeff_region_mod, ONLY: two_coeff_region
  USE two_coeff_region_fast_lw_mod, ONLY: two_coeff_region_fast_lw

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
    , nd_layer_clr                                                      &
!       Size allocated for completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_cloud_type                                                     &
!       Maximum number of types of cloud
    , nd_region                                                         &
!       Maximum number of cloudy regions
    , nd_overlap_coeff
!       Maximum number of overlap coefficients

! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Top cloudy layer
    , n_cloud_type                                                      &
!       Number of types of clouds
    , isolir                                                            &
!       Spectral region
    , i_2stream                                                         &
!       Two-stream scheme
    , i_solver                                                          &
!       Solver used
    , i_solver_clear                                                    &
!       Solver for clear-sky fluxes
    , i_scatter_method
!       Method of treating scattering
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  LOGICAL, INTENT(IN) ::                                                &
      l_clear                                                           &
!       Calculate clear-sky fluxes
    , l_scale_solar                                                     &
!       Flag to scale solar
    , l_ir_source_quad
!       Use quadratic source term

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!       Spherical geometry fields

! Optical properties:
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

! Cloud geometry:
  INTEGER, INTENT(IN) ::                                                &
      n_region                                                          &
!       Number of cloudy regions
    , i_region_cloud(nd_cloud_type)
!       Regions in which types of clouds fall
  REAL (RealK), INTENT(IN) ::                                           &
      w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Cloudy fractions in each layer
    , w_free(nd_profile, id_ct: nd_layer)                               &
!       Clear sky fractions in each layer
    , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)            &
!       Fractions of different types of cloud
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)                                             &
!       Energy transfer coefficients
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region
  REAL (RealK), INTENT(IN) ::                                           &
      sec_0(nd_profile)                                                 &
!       Secant of solar zenith angle
    , diffuse_albedo(nd_profile)                                        &
!       Diffuse albedo
    , direct_albedo(nd_profile)                                         &
!       Direct albedo
    , flux_inc_down(nd_profile)                                         &
!       Incident total flux
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , diff_planck(nd_profile, nd_layer)                                 &
!       Change in Planckian function
    , d_planck_flux_surface(nd_profile)                                 &
!       Difference in Planckian fluxes between the surface
!       and the air above
    , adjust_solar_ke(nd_profile, nd_layer)                             &
!       Adjustment of solar beam with equivalent extinction
    , diff_planck_2(nd_profile, nd_layer)
!         2x2nd difference of Planckian

! Fluxes calculated
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_total(nd_profile, 2*nd_layer+2)                              &
!       Long flux vector
    , flux_direct_clear(nd_profile, 0: nd_layer)                        &
!       Clear direct flux
    , flux_total_clear(nd_profile, 2*nd_layer+2)
!       Clear total flux



! Local variabales.
  INTEGER                                                               &
      n_source_coeff                                                    &
!       Number of source coefficients
    , i                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , n_top
!       Top-most layer for calculation


! Clear-sky coefficients:
  REAL (RealK) ::                                                       &
      trans(nd_profile, nd_layer, nd_region)                            &
!       Transmission coefficients
    , reflect(nd_profile, nd_layer, nd_region)                          &
!       Reflection coefficients
    , trans_0(nd_profile, nd_layer, nd_region)                          &
!       Direct transmission coefficients
    , trans_0_dir(nd_profile, nd_layer, nd_region)                      &
!       Unscaled Direct transmission coefficients
    , source_coeff(nd_profile, nd_layer, nd_source_coeff, nd_region)    &
!       Source coefficients
    , s_down(nd_profile, nd_layer, nd_region)                           &
!       Free downward source
    , s_up(nd_profile, nd_layer, nd_region)                             &
!       Free upward source
    , s_down_clear(nd_profile, nd_layer)                                &
!       Clear downward source
    , s_up_clear(nd_profile, nd_layer)
!       Clear upward source

! Source functions at the ground
  REAL (RealK) ::                                                       &
      source_flux_ground(nd_profile, nd_region)                         &
!       Source of flux from ground
    , flux_direct_ground(nd_profile, nd_region)
!       Direct flux at ground in each region


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'TRIPLE_COLUMN'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the number of source coefficients for the approximation
  n_source_coeff=set_n_source_coeff(isolir, l_ir_source_quad)


  IF ( (i_scatter_method == ip_scatter_full) .OR.                       &
       (i_scatter_method == ip_scatter_approx) ) THEN
    CALL two_coeff_region(ierr, control                                 &
      , n_profile, n_layer, n_cloud_top                                 &
      , i_2stream, n_source_coeff                                       &
      , n_cloud_type, frac_cloud                                        &
      , n_region, i_region_cloud, frac_region                           &
      , ss_prop%phase_fnc_clr, ss_prop%omega_clr                        &
      , ss_prop%tau_clr_dir, ss_prop%tau_clr                            &
      , ss_prop%phase_fnc, ss_prop%omega                                &
      , ss_prop%tau_dir, ss_prop%tau                                    &
      , isolir, sec_0, sph                                              &
      , trans, reflect, trans_0_dir, trans_0, source_coeff              &
      , nd_profile, nd_layer, nd_layer_clr, id_ct                       &
      , nd_max_order, nd_source_coeff                                   &
      , nd_cloud_type, nd_region                                        &
      )
  ELSE IF ( (i_scatter_method == ip_no_scatter_abs) .OR.                &
            (i_scatter_method == ip_no_scatter_ext) ) THEN
    CALL two_coeff_region_fast_lw(ierr                                  &
      , n_profile, n_layer, n_cloud_top                                 &
      , l_ir_source_quad, n_source_coeff                                &
      , n_cloud_type, frac_cloud                                        &
      , n_region, i_region_cloud, frac_region                           &
      , ss_prop%tau_clr, ss_prop%tau                                    &
      , isolir                                                          &
      , trans, reflect, source_coeff                                    &
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_source_coeff      &
      , nd_cloud_type, nd_region                                        &
      )
  END IF


  IF (isolir == ip_infra_red) THEN

    DO k=1, n_region
      IF (k == ip_region_clear) THEN
        n_top=1
      ELSE
        n_top=n_cloud_top
      END IF

      CALL ir_source(n_profile, n_top, n_layer                          &
        , source_coeff(1, 1, 1, k), diff_planck                         &
        , l_ir_source_quad, diff_planck_2                               &
        , s_down(1, 1, k), s_up(1, 1, k)                                &
        , nd_profile, nd_layer, nd_source_coeff                         &
        )
    END DO


!   Weight the source functions by the area fractions, but
!   save the clear-sky fractions for diagnostic use if
!   required.
    IF (l_clear) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          s_down_clear(l, i)=s_down(l, i, ip_region_clear)
          s_up_clear(l, i)=s_up(l, i, ip_region_clear)
        END DO
      END DO
    END IF
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        s_down(l, i, ip_region_clear)                                   &
          =w_free(l, i)*s_down(l, i, ip_region_clear)
        s_up(l, i, ip_region_clear)                                     &
          =w_free(l, i)*s_up(l, i, ip_region_clear)
        s_down(l, i, ip_region_strat)                                   &
          =w_cloud(l, i)                                                &
          *frac_region(l, i, ip_region_strat)                           &
          *s_down(l, i, ip_region_strat)
        s_up(l, i, ip_region_strat)                                     &
          =w_cloud(l, i)                                                &
          *frac_region(l, i, ip_region_strat)                           &
          *s_up(l, i, ip_region_strat)
        s_down(l, i, ip_region_conv)                                    &
          =w_cloud(l, i)                                                &
          *frac_region(l, i, ip_region_conv)                            &
          *s_down(l, i, ip_region_conv)
        s_up(l, i, ip_region_conv)                                      &
          =w_cloud(l, i)                                                &
          *frac_region(l, i, ip_region_conv)                            &
          *s_up(l, i, ip_region_conv)
      END DO
    END DO

  END IF


! Calculate the appropriate source terms for the solar: cloudy
! and clear properties are both needed here.

  IF (isolir == ip_solar) THEN

    CALL triple_solar_source(control, bound                             &
      , n_profile, n_layer, n_cloud_top                                 &
      , n_region, flux_inc_direct                                       &
      , l_scale_solar, adjust_solar_ke                                  &
      , trans_0_dir, trans_0, source_coeff                              &
      , cloud_overlap(1, id_ct-1, 1)                                    &
      , cloud_overlap(1, id_ct-1, 2)                                    &
      , cloud_overlap(1, id_ct-1, 3)                                    &
      , cloud_overlap(1, id_ct-1, 4)                                    &
      , cloud_overlap(1, id_ct-1, 5)                                    &
      , cloud_overlap(1, id_ct-1, 6)                                    &
      , cloud_overlap(1, id_ct-1, 7)                                    &
      , cloud_overlap(1, id_ct-1, 8)                                    &
      , cloud_overlap(1, id_ct-1, 9)                                    &
      , flux_direct, flux_direct_ground                                 &
      , s_up, s_down                                                    &
      , nd_profile, nd_layer, id_ct, nd_source_coeff, nd_region         &
      )
  END IF

!   Set the partitioned source functions at the ground.
    IF (isolir == ip_solar) THEN
      DO l=1, n_profile
        source_flux_ground(l, ip_region_clear)                          &
          =(direct_albedo(l)-diffuse_albedo(l))                         &
          *flux_direct_ground(l, ip_region_clear)
        source_flux_ground(l, ip_region_strat)                          &
          =(direct_albedo(l)-diffuse_albedo(l))                         &
          *flux_direct_ground(l, ip_region_strat)
        source_flux_ground(l, ip_region_conv)                           &
          =(direct_albedo(l)-diffuse_albedo(l))                         &
          *flux_direct_ground(l, ip_region_conv)
      END DO
    ELSE
      DO l=1, n_profile
        source_flux_ground(l, ip_region_clear)                          &
          =cloud_overlap(l, n_layer, 10)                                &
          *(1.0e+00_RealK-diffuse_albedo(l))                            &
          *d_planck_flux_surface(l)
        source_flux_ground(l, ip_region_strat)                          &
          =cloud_overlap(l, n_layer, 13)                                &
          *(1.0e+00_RealK-diffuse_albedo(l))                            &
          *d_planck_flux_surface(l)
        source_flux_ground(l, ip_region_conv)                           &
          =cloud_overlap(l, n_layer, 16)                                &
          *(1.0e+00_RealK-diffuse_albedo(l))                            &
          *d_planck_flux_surface(l)
      END DO
    END IF



  SELECT CASE (i_solver)

  CASE (ip_solver_triple)

    CALL solver_triple(n_profile, n_layer, n_cloud_top                  &
      , trans(1, 1, ip_region_clear)                                    &
      , reflect(1, 1, ip_region_clear)                                  &
      , s_down(1, 1, ip_region_clear)                                   &
      , s_up(1, 1, ip_region_clear)                                     &
      , trans(1, 1, ip_region_strat)                                    &
      , reflect(1, 1, ip_region_strat)                                  &
      , s_down(1, 1, ip_region_strat)                                   &
      , s_up(1, 1, ip_region_strat)                                     &
      , trans(1, 1, ip_region_conv)                                     &
      , reflect(1, 1, ip_region_conv)                                   &
      , s_down(1, 1, ip_region_conv)                                    &
      , s_up(1, 1, ip_region_conv)                                      &
      , cloud_overlap(1, id_ct-1, 1)                                    &
      , cloud_overlap(1, id_ct-1, 2)                                    &
      , cloud_overlap(1, id_ct-1, 3)                                    &
      , cloud_overlap(1, id_ct-1, 4)                                    &
      , cloud_overlap(1, id_ct-1, 5)                                    &
      , cloud_overlap(1, id_ct-1, 6)                                    &
      , cloud_overlap(1, id_ct-1, 7)                                    &
      , cloud_overlap(1, id_ct-1, 8)                                    &
      , cloud_overlap(1, id_ct-1, 9)                                    &
      , cloud_overlap(1, id_ct-1, 10)                                   &
      , cloud_overlap(1, id_ct-1, 11)                                   &
      , cloud_overlap(1, id_ct-1, 12)                                   &
      , cloud_overlap(1, id_ct-1, 13)                                   &
      , cloud_overlap(1, id_ct-1, 14)                                   &
      , cloud_overlap(1, id_ct-1, 15)                                   &
      , cloud_overlap(1, id_ct-1, 16)                                   &
      , cloud_overlap(1, id_ct-1, 17)                                   &
      , cloud_overlap(1, id_ct-1, 18)                                   &
      , flux_inc_down                                                   &
      , source_flux_ground(1, ip_region_clear)                          &
      , source_flux_ground(1, ip_region_strat)                          &
      , source_flux_ground(1, ip_region_conv)                           &
      , diffuse_albedo                                                  &
      , flux_total                                                      &
      , nd_profile, nd_layer, id_ct                                     &
      )

  CASE (ip_solver_triple_hogan)

    CALL solver_triple_hogan(n_profile, n_layer, n_cloud_top            &
      , trans(1, 1, ip_region_clear)                                    &
      , reflect(1, 1, ip_region_clear)                                  &
      , s_down(1, 1, ip_region_clear)                                   &
      , s_up(1, 1, ip_region_clear)                                     &
      , trans(1, 1, ip_region_strat)                                    &
      , reflect(1, 1, ip_region_strat)                                  &
      , s_down(1, 1, ip_region_strat)                                   &
      , s_up(1, 1, ip_region_strat)                                     &
      , trans(1, 1, ip_region_conv)                                     &
      , reflect(1, 1, ip_region_conv)                                   &
      , s_down(1, 1, ip_region_conv)                                    &
      , s_up(1, 1, ip_region_conv)                                      &
      , cloud_overlap(1, id_ct-1, 1)                                    &
      , cloud_overlap(1, id_ct-1, 2)                                    &
      , cloud_overlap(1, id_ct-1, 3)                                    &
      , cloud_overlap(1, id_ct-1, 4)                                    &
      , cloud_overlap(1, id_ct-1, 5)                                    &
      , cloud_overlap(1, id_ct-1, 6)                                    &
      , cloud_overlap(1, id_ct-1, 7)                                    &
      , cloud_overlap(1, id_ct-1, 8)                                    &
      , cloud_overlap(1, id_ct-1, 9)                                    &
      , cloud_overlap(1, id_ct-1, 10)                                   &
      , cloud_overlap(1, id_ct-1, 11)                                   &
      , cloud_overlap(1, id_ct-1, 12)                                   &
      , cloud_overlap(1, id_ct-1, 13)                                   &
      , cloud_overlap(1, id_ct-1, 14)                                   &
      , cloud_overlap(1, id_ct-1, 15)                                   &
      , cloud_overlap(1, id_ct-1, 16)                                   &
      , cloud_overlap(1, id_ct-1, 17)                                   &
      , cloud_overlap(1, id_ct-1, 18)                                   &
      , flux_inc_down                                                   &
      , source_flux_ground(1, ip_region_clear)                          &
      , source_flux_ground(1, ip_region_strat)                          &
      , source_flux_ground(1, ip_region_conv)                           &
      , diffuse_albedo                                                  &
      , flux_total                                                      &
      , nd_profile, nd_layer, id_ct                                     &
      )

  CASE (ip_solver_triple_app_scat)

    CALL solver_triple_app_scat(n_profile, n_layer, n_cloud_top         &
      , trans(1, 1, ip_region_clear)                                    &
      , reflect(1, 1, ip_region_clear)                                  &
      , s_down(1, 1, ip_region_clear)                                   &
      , s_up(1, 1, ip_region_clear)                                     &
      , trans(1, 1, ip_region_strat)                                    &
      , reflect(1, 1, ip_region_strat)                                  &
      , s_down(1, 1, ip_region_strat)                                   &
      , s_up(1, 1, ip_region_strat)                                     &
      , trans(1, 1, ip_region_conv)                                     &
      , reflect(1, 1, ip_region_conv)                                   &
      , s_down(1, 1, ip_region_conv)                                    &
      , s_up(1, 1, ip_region_conv)                                      &
      , cloud_overlap(1, id_ct-1, 1)                                    &
      , cloud_overlap(1, id_ct-1, 2)                                    &
      , cloud_overlap(1, id_ct-1, 3)                                    &
      , cloud_overlap(1, id_ct-1, 4)                                    &
      , cloud_overlap(1, id_ct-1, 5)                                    &
      , cloud_overlap(1, id_ct-1, 6)                                    &
      , cloud_overlap(1, id_ct-1, 7)                                    &
      , cloud_overlap(1, id_ct-1, 8)                                    &
      , cloud_overlap(1, id_ct-1, 9)                                    &
      , cloud_overlap(1, id_ct-1, 10)                                   &
      , cloud_overlap(1, id_ct-1, 11)                                   &
      , cloud_overlap(1, id_ct-1, 12)                                   &
      , cloud_overlap(1, id_ct-1, 13)                                   &
      , cloud_overlap(1, id_ct-1, 14)                                   &
      , cloud_overlap(1, id_ct-1, 15)                                   &
      , cloud_overlap(1, id_ct-1, 16)                                   &
      , cloud_overlap(1, id_ct-1, 17)                                   &
      , cloud_overlap(1, id_ct-1, 18)                                   &
      , flux_inc_down                                                   &
      , source_flux_ground(1, ip_region_clear)                          &
      , source_flux_ground(1, ip_region_strat)                          &
      , source_flux_ground(1, ip_region_conv)                           &
      , diffuse_albedo                                                  &
      , flux_total                                                      &
      , nd_profile, nd_layer, id_ct                                     &
      )

  CASE DEFAULT

    cmessage = '***Error: The solver specified is not valid here.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END SELECT



  IF (l_clear) THEN

    CALL column_solver(ierr, control, bound, sph%common, sph%clear      &
      , n_profile, n_layer                                              &
      , i_scatter_method, i_solver_clear                                &
      , trans(1, 1, ip_region_clear)                                    &
      , reflect(1, 1, ip_region_clear)                                  &
      , trans_0_dir(1, 1, ip_region_clear)                              &
      , trans_0(1, 1, ip_region_clear)                                  &
      , source_coeff(1, 1, 1, ip_region_clear)                          &
      , isolir, flux_inc_direct, flux_inc_down                          &
      , s_down_clear, s_up_clear                                        &
      , diffuse_albedo, direct_albedo                                   &
      , d_planck_flux_surface                                           &
      , l_scale_solar, adjust_solar_ke                                  &
      , flux_direct_clear, flux_total_clear                             &
      , nd_profile, nd_layer, nd_source_coeff                           &
      )
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triple_column
END MODULE triple_column_mod

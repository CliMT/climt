! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve the two-stream equations in a mixed column.
!
! Method:
!   The two-stream coefficients are calculated in clear regions
!   and in stratiform and convective clouds. From these
!   coefficients transmission and reflection coefficients are
!   determined. The coefficients for convective and stratiform
!   clouds are appropriately mixed to form single cloudy values
!   and an appropriate solver is called.
!
!- ---------------------------------------------------------------------
MODULE mix_column_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MIX_COLUMN_MOD'
CONTAINS
SUBROUTINE mix_column(ierr                                              &
     , control, bound                                                   &
!                 Atmospheric properties
     , n_profile, n_layer, k_clr                                        &
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
     , w_free, w_cloud                                                  &
     , cloud_overlap                                                    &
!                 Calculated fluxes
     , flux_direct, flux_total                                          &
!                 Flags for clear-sky calculations
     , l_clear, i_solver_clear                                          &
!                 Calculated clear-sky fluxes
     , flux_direct_clear, flux_total_clear                              &
!                 Dimensions of arrays
     , nd_profile, nd_layer, nd_layer_clr, id_ct                        &
     , nd_max_order, nd_source_coeff                                    &
     , nd_cloud_type, nd_overlap_coeff                                  &
     )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_solver_mix_app_scat,ip_solver_mix_direct,       &
                     ip_solver_mix_direct_hogan, ip_solar, ip_infra_red,&
                     i_err_fatal
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE column_solver_mod, ONLY: column_solver
  USE set_n_source_coeff_mod, ONLY: set_n_source_coeff
  USE ir_source_mod, ONLY: ir_source
  USE mix_app_scat_mod, ONLY: mix_app_scat
  USE mixed_solar_source_mod, ONLY: mixed_solar_source
  USE solver_mix_direct_mod, ONLY: solver_mix_direct
  USE solver_mix_direct_hogan_mod, ONLY: solver_mix_direct_hogan
  USE two_coeff_mod, ONLY: two_coeff
  USE two_coeff_cloud_mod, ONLY: two_coeff_cloud

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_layer_clr                                                      &
!       Size allocated for completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_source_coeff                                                   &
!       Size allocated for coefficients in the source function
    , nd_cloud_type                                                     &
!       Size allocated for types of clouds
    , nd_overlap_coeff
!       Size allocated for overlpa coefficients

! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Top cloudy layer
    , k_clr                                                             &
!       Index of the clear-sky region
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

! Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

! Cloud geometry:
  REAL (RealK), INTENT(IN) ::                                           &
      w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Cloudy fractions in each layer
    , w_free(nd_profile, id_ct: nd_layer)                               &
!       Clear sky fractions in each layer
    , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)            &
!       Fractions of different types of cloud
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)
!       Energy transfer coefficients
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
!       Flux from surface
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
    , l
!       Loop variable

! Pointers to sections of the array of overlap coefficients:
! Here F denotes the clear-sky region and C the cloudy region;
! the ordering of the final suffix is such that the suffix for
! the area being left appears last: this is convenient to agree
! with the documented notation, but is not directly consistent
! with older versions of the code.
  INTEGER                                                               &
      i_ovp_dn_ff                                                       &
!       Pointer to section of the array of overlaps for downward
!       transmission from clear-sky to clear-sky
    , i_ovp_dn_fc                                                       &
!       Pointer to section of the array of overlaps for downward
!       transmission from cloud to clear-sky
    , i_ovp_dn_cf                                                       &
!       Pointer to section of the array of overlaps for downward
!       transmission from clear-sky to cloud
    , i_ovp_dn_cc                                                       &
!       Pointer to section of the array of overlaps for downward
!       transmission from cloud to cloud
    , i_ovp_up_ff                                                       &
!       Pointer to section of the array of overlaps for upward
!       transmission from clear-sky to clear-sky
    , i_ovp_up_fc                                                       &
!       Pointer to section of the array of overlaps for upward
!       transmission from cloud to clear-sky
    , i_ovp_up_cf                                                       &
!       Pointer to section of the array of overlaps for upward
!       transmission from clear-sky to cloud
    , i_ovp_up_cc
!       Pointer to section of the array of overlaps for upward
!       transmission from cloud to cloud


! Clear-sky coefficients:
  REAL (RealK) ::                                                       &
      trans_free(nd_profile, nd_layer)                                  &
!       Free transmission of layer
    , reflect_free(nd_profile, nd_layer)                                &
!       Free reflectance of layer
    , trans_0_free(nd_profile, nd_layer)                                &
!       Free direct transmission of layer
    , trans_0_free_dir(nd_profile, nd_layer)                            &
!       Free direct transmission of layer without scaling
    , source_coeff_free(nd_profile, nd_layer, nd_source_coeff)          &
!       Free source coefficients
    , s_down_free(nd_profile, nd_layer)                                 &
!       Free downward source
    , s_up_free(nd_profile, nd_layer)                                   &
!       Free upward source
    , s_down_clear(nd_profile, nd_layer)                                &
!       Clear downward source
    , s_up_clear(nd_profile, nd_layer)
!       Clear upward source

! Cloudy coefficients:
  REAL (RealK) ::                                                       &
      trans_cloud(nd_profile, nd_layer)                                 &
!       Cloudy transmission of layer
    , reflect_cloud(nd_profile, nd_layer)                               &
!       Cloudy reflectance of layer
    , trans_0_cloud(nd_profile, nd_layer)                               &
!       Cloudy direct transmission of layer
    , trans_0_cloud_dir(nd_profile, nd_layer)                           &
!       Cloudy direct transmission of layer without scaling
    , source_coeff_cloud(nd_profile, nd_layer, nd_source_coeff)         &
!       Cloudy source coefficients
    , s_down_cloud(nd_profile, nd_layer)                                &
!       Cloudy downward source
    , s_up_cloud(nd_profile, nd_layer)
!       Cloudy upward source

! Source functions at the surface
  REAL (RealK) ::                                                       &
      source_ground_free(nd_profile)                                    &
!       Source from ground under clear skies
    , source_ground_cloud(nd_profile)                                   &
!       Source from ground under cloudy skies
    , flux_direct_ground_cloud(nd_profile)
!       Direct flux at ground under cloudy skies

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'MIX_COLUMN'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the pointers to the various types of transition.
  i_ovp_dn_ff=3*k_clr-2
  i_ovp_dn_fc=k_clr+1
  i_ovp_dn_cf=4-k_clr
  i_ovp_dn_cc=7-3*k_clr
  i_ovp_up_ff=4+i_ovp_dn_ff
  i_ovp_up_fc=4+i_ovp_dn_fc
  i_ovp_up_cf=4+i_ovp_dn_cf
  i_ovp_up_cc=4+i_ovp_dn_cc

! Calculate the transmission and reflection coefficients and
! source terms for the clear and cloudy parts of the column

! Set the number of source coefficients for the approximation
  n_source_coeff=set_n_source_coeff(isolir, l_ir_source_quad)

  CALL two_coeff(ierr, control                                          &
    , n_profile, 1, n_cloud_top-1                                       &
    , i_2stream                                                         &
    , ss_prop%phase_fnc_clr                                             &
    , ss_prop%omega_clr, ss_prop%tau_clr_dir, ss_prop%tau_clr           &
    , isolir, sec_0, sph%common%path_div                                &
    , trans_free, reflect_free, trans_0_free_dir, trans_0_free          &
    , source_coeff_free                                                 &
    , nd_profile, 1, nd_layer_clr, 1, nd_layer, nd_source_coeff         &
    )
  CALL two_coeff(ierr, control                                          &
    , n_profile, n_cloud_top, n_layer                                   &
    , i_2stream                                                         &
    , ss_prop%phase_fnc(:, :, :, 0)                                     &
    , ss_prop%omega(:, :, 0)                                            &
    , ss_prop%tau_dir(:, :, 0), ss_prop%tau(:, :, 0)                    &
    , isolir, sec_0, sph%common%path_div                                &
    , trans_free, reflect_free, trans_0_free_dir, trans_0_free          &
    , source_coeff_free                                                 &
    , nd_profile, id_ct, nd_layer, 1, nd_layer, nd_source_coeff         &
    )


! Infra-red source terms depend only on the layer and may be
! calculated now. Solar terms depend on conditions in cloud
! in overlying layers and must be calculated later.

  IF (isolir == ip_infra_red) THEN

    CALL ir_source(n_profile, 1, n_layer                                &
      , source_coeff_free, diff_planck                                  &
      , l_ir_source_quad, diff_planck_2                                 &
      , s_down_free, s_up_free                                          &
      , nd_profile, nd_layer, nd_source_coeff                           &
      )

!   If a clear-sky calculation is required these source terms must
!   be stored.
    IF (l_clear) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          s_down_clear(l, i)=s_down_free(l, i)
          s_up_clear(l, i)=s_up_free(l, i)
        END DO
      END DO
    END IF

!   Scale the sources by the clear-sky fractions in the cloudy
!   layers. In higher layers the clear-sky fraction is 1.
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        s_down_free(l, i)=w_free(l, i)*s_down_free(l, i)
        s_up_free(l, i)=w_free(l, i)*s_up_free(l, i)
      END DO
    END DO

  END IF



! Repeat the calculation for cloudy regions.

! Clouds are indexed beginning with index 1 in the last
! dimension of arrays of optical properties.


  CALL two_coeff_cloud(ierr, control                                    &
    , n_profile, n_cloud_top, n_layer                                   &
    , i_2stream, n_source_coeff                                         &
    , n_cloud_type, frac_cloud                                          &
    , ss_prop%phase_fnc(:, :, :, 1:)                                    &
    , ss_prop%omega(:, :, 1:)                                           &
    , ss_prop%tau_dir(:, :, 1:), ss_prop%tau(:, :, 1:)                  &
    , isolir, sec_0, sph                                                &
    , trans_cloud, reflect_cloud, trans_0_cloud_dir, trans_0_cloud      &
    , source_coeff_cloud                                                &
    , nd_profile, nd_layer, id_ct, nd_max_order                         &
    , nd_source_coeff, nd_cloud_type                                    &
    )


  IF (isolir == ip_infra_red) THEN

    CALL ir_source(n_profile, n_cloud_top, n_layer                      &
      , source_coeff_cloud, diff_planck                                 &
      , l_ir_source_quad, diff_planck_2                                 &
      , s_down_cloud, s_up_cloud                                        &
      , nd_profile, nd_layer, nd_source_coeff                           &
      )

    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        s_down_cloud(l, i)=w_cloud(l, i)*s_down_cloud(l, i)
        s_up_cloud(l, i)=w_cloud(l, i)*s_up_cloud(l, i)
      END DO
    END DO

  END IF


! Calculate the appropriate source terms for the solar: cloudy
! and clear properties are both needed here.

  IF (isolir == ip_solar) THEN

    CALL mixed_solar_source(control, bound                              &
      , n_profile, n_layer, n_cloud_top                                 &
      , flux_inc_direct                                                 &
      , l_scale_solar, adjust_solar_ke                                  &
      , trans_0_free_dir, trans_0_free, source_coeff_free               &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cc)                          &
      , trans_0_cloud_dir, trans_0_cloud, source_coeff_cloud            &
      , flux_direct                                                     &
      , flux_direct_ground_cloud                                        &
      , s_up_free, s_down_free                                          &
      , s_up_cloud, s_down_cloud                                        &
      , nd_profile, nd_layer, id_ct, nd_source_coeff                    &
    )
  END IF



! Formulate the matrix equation for the fluxes.

  SELECT CASE (i_solver)

  CASE (ip_solver_mix_app_scat)

    CALL mix_app_scat(n_profile, n_layer, n_cloud_top                   &
      , trans_free, reflect_free, s_down_free, s_up_free                &
      , trans_cloud, reflect_cloud                                      &
      , s_down_cloud, s_up_cloud                                        &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_cc)                          &
      , flux_inc_down                                                   &
      , d_planck_flux_surface, diffuse_albedo                           &
      , flux_total                                                      &
      , nd_profile, nd_layer, id_ct                                     &
      )

  CASE (ip_solver_mix_direct, ip_solver_mix_direct_hogan)

!   Set the partitioned source functions at the ground.
    IF (isolir == ip_solar) THEN
      DO l=1, n_profile
        source_ground_free(l)=(direct_albedo(l)                         &
          -diffuse_albedo(l))                                           &
          *(flux_direct(l, n_layer)                                     &
          -flux_direct_ground_cloud(l))
        source_ground_cloud(l)=(direct_albedo(l)                        &
          -diffuse_albedo(l))                                           &
          *flux_direct_ground_cloud(l)
      END DO
    ELSE
      DO l=1, n_profile
        source_ground_free(l)                                           &
          =cloud_overlap(l, n_layer, i_ovp_up_ff)                       &
          *(1.0_RealK-diffuse_albedo(l))*d_planck_flux_surface(l)
        source_ground_cloud(l)                                          &
          =cloud_overlap(l, n_layer, i_ovp_up_cf)                       &
          *(1.0_RealK-diffuse_albedo(l))*d_planck_flux_surface(l)
      END DO
    END IF

    IF (i_solver == ip_solver_mix_direct) THEN
      CALL solver_mix_direct(n_profile, n_layer, n_cloud_top            &
      , trans_free, reflect_free, s_down_free, s_up_free                &
      , trans_cloud, reflect_cloud                                      &
      , s_down_cloud, s_up_cloud                                        &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_cc)                          &
      , flux_inc_down                                                   &
      , source_ground_free, source_ground_cloud                         &
      , diffuse_albedo                                                  &
      , flux_total                                                      &
      , nd_profile, nd_layer, id_ct                                     &
      )

    ELSE IF (i_solver == ip_solver_mix_direct_hogan) THEN
      CALL solver_mix_direct_hogan(n_profile, n_layer, n_cloud_top      &
      , trans_free, reflect_free, s_down_free, s_up_free                &
      , trans_cloud, reflect_cloud                                      &
      , s_down_cloud, s_up_cloud                                        &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_dn_cc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_ff)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_fc)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_cf)                          &
      , cloud_overlap(1, id_ct-1, i_ovp_up_cc)                          &
      , flux_inc_down                                                   &
      , source_ground_free, source_ground_cloud                         &
      , diffuse_albedo                                                  &
      , flux_total                                                      &
      , nd_profile, nd_layer, id_ct                                     &
      )

    END IF

  CASE DEFAULT

    cmessage = '***Error: The solver specified is not valid here.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END SELECT



  IF (l_clear) THEN

    CALL column_solver(ierr, control, bound, sph%common, sph%clear      &
      , n_profile, n_layer                                              &
      , i_scatter_method, i_solver_clear                                &
      , trans_free, reflect_free, trans_0_free_dir                      &
      , trans_0_free, source_coeff_free                                 &
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

END SUBROUTINE mix_column
END MODULE mix_column_mod

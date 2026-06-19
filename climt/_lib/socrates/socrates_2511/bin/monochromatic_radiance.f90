! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve for the monochromatic radiances.
!
! Method:
!   The final single scattering properties are calculated
!   and rescaled. An appropriate subroutine is called to
!   calculate the radiances depending on the treatment of
!   cloudiness.
!
!- ---------------------------------------------------------------------
MODULE monochromatic_radiance_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MONOCHROMATIC_RADIANCE_MOD'
CONTAINS
SUBROUTINE monochromatic_radiance(ierr                                  &
    , control, atm, cld, bound                                          &
!                 Atmospheric Propertries
    , n_profile, n_layer, d_mass                                        &
!                 Angular Integration
    , i_angular_integration, i_2stream                                  &
    , l_rescale, n_order_gauss                                          &
    , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc       &
    , accuracy_adaptive, euler_factor, i_sph_algorithm                  &
    , i_sph_mode                                                        &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Treatment of Scattering
    , i_scatter_method                                                  &
!                 Options for Solver
    , i_solver                                                          &
!                 Gaseous Properties
    , n_k_term_inner, k_gas_abs                                         &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , planck                                                            &
!                 Conditions at TOA
    , zen_0, flux_inc_direct, flux_inc_down                             &
    , i_direct                                                          &
!                 Surface Properties
    , d_planck_flux_surface                                             &
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
!                 Spherical geometry
    , sph                                                               &
!                 Optical Properties
    , ss_prop                                                           &
!                 Cloudy Properties
    , l_cloud, i_cloud                                                  &
!                 Cloud Geometry
    , n_cloud_top, index_subcol                                         &
    , n_region, k_clr, i_region_cloud, frac_region                      &
    , w_free, cloud_overlap                                             &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                 Viewing geometry
    , n_direction, direction                                            &
!                 Calculated Fluxes
    , flux_direct, flux_total, l_actinic, actinic_flux                  &
!                 Calculated Radiances
    , radiance                                                          &
!                 Calculated mean radiances
    , j_radiance                                                        &
!                 Flags for Clear-sky Calculation
    , l_clear, i_solver_clear                                           &
!                 Clear-sky Fluxes Calculated
    , flux_direct_clear, flux_total_clear, actinic_flux_clear           &
!                 Contribution function
    , contrib_funci_part, contrib_funcf_part                            &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_cloud_type, nd_region, nd_overlap_coeff                        &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                &
    , nd_direction, nd_source_coeff                                     &
    , nd_k_term_inner                                                   &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_atm,     ONLY: StrAtm
  USE def_cld,     ONLY: StrCld
  USE def_bound,   ONLY: StrBound
  USE def_planck,  ONLY: StrPlanck
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_direct_csr_scaling, ip_direct_noscaling,        &
                     ip_ir_gauss, ip_spherical_harmonic, ip_two_stream, &
                     ip_surf_alb_diff
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE calc_contrib_func_mod, ONLY: calc_contrib_func 
  USE spherical_trans_coeff_mod, ONLY: spherical_trans_coeff
  USE gauss_angle_mod, ONLY: gauss_angle
  USE monochromatic_radiance_sph_mod, ONLY: monochromatic_radiance_sph
  USE monochromatic_radiance_tseq_mod, ONLY: monochromatic_radiance_tseq
  USE rescale_tau_csr_mod, ONLY: rescale_tau_csr
  USE rescale_tau_omega_mod, ONLY: rescale_tau_omega
  USE single_scattering_all_mod, ONLY: single_scattering_all
  USE single_scattering_mod, ONLY: single_scattering

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Atmospheric properties:
  TYPE(StrAtm),       INTENT(IN)    :: atm

! Cloud properties:
  TYPE(StrCld),       INTENT(IN)    :: cld

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_layer_clr                                                      &
!       Maximum number of completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_flux_profile                                                   &
!       Maximum number of profiles in arrays of fluxes
    , nd_radiance_profile                                               &
!       Maximum number of profiles in arrays of radiances
    , nd_j_profile                                                      &
!       Maximum number of profiles in arrays of mean radiances
    , nd_column                                                         &
!       Number of columns per point
    , nd_cloud_type                                                     &
!       Maximum number of types of cloud
    , nd_region                                                         &
!       Maximum number of cloudy regions
    , nd_overlap_coeff                                                  &
!       Maximum number of overlap coefficients
    , nd_max_order                                                      &
!       Maximum order of spherical harmonics used
    , nd_sph_coeff                                                      &
!       Allocated size for spherical coefficients
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for orders of BRDFs
    , nd_viewing_level                                                  &
!       Allocated size for levels where radiances are calculated
    , nd_direction                                                      &
!       Allocated size for viewing directions
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_k_term_inner
!       Size allocated for k-terms in inner loop

! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

!                 Atmospheric properties
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                           &
      d_mass(nd_profile, nd_layer)
!       Mass thickness of each layer

!                 Angular integration
  INTEGER, INTENT(IN) ::                                                &
      i_angular_integration                                             &
!       Angular integration scheme
    , i_2stream                                                         &
!       Two-stream scheme
    , n_order_gauss                                                     &
!       Order of Gaussian integration
    , n_order_phase                                                     &
!       Highest order retained in the phase function
    , i_truncation                                                      &
!       Type of spherical truncation adopted
    , ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                   &
!       Orders of truncation at each azimuthal order
    , i_sph_mode                                                        &
!       Mode in which teh spherical harmonic solver is being used
    , i_sph_algorithm
!       Algorithm used for spherical harmonic calculation
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale
!       Rescale optical properties
  REAL (RealK) ::                                                       &
      cg_coeff(nd_sph_coeff)                                            &
!       Clebsch-Gordan coefficients
    , uplm_zero(nd_sph_coeff)                                           &
!       Values of spherical harmonics at polar angles pi/2
    , uplm_sol(nd_radiance_profile, nd_sph_coeff)
!       Values of spherical harmonics in the solar direction
  REAL (RealK), INTENT(IN) ::                                           &
      accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series

!                 Treatment of scattering
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method

!                 Options for solver
  INTEGER, INTENT(IN) ::                                                &
      i_solver
!       Solver used

!                 Gaseous properties
  INTEGER, INTENT(IN) ::                                                &
      n_k_term_inner
!       Number of k-terms in inner loop
  REAL (RealK), INTENT(IN) ::                                           &
      k_gas_abs(nd_profile, nd_layer, nd_k_term_inner)
!       Gaseous absorptive extinctions

!                 Variables for equivalent extinction
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Apply scaling to solar flux
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment of solar beam with equivalent extinction

!                 Spectral region
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Visible or IR

!                 Conditions at TOA
  REAL (RealK), INTENT(IN) ::                                           &
      zen_0(nd_profile)                                                 &
!       Secants (two-stream) or cosines (spherical harmonics)
!       of the solar zenith angles
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)
!       Incident downward flux
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_radiance_profile, 0: nd_layer, nd_k_term_inner)
!       Direct radiance (the first row contains the incident
!       solar radiance: the other rows are calculated)

!                 Surface properties
  REAL (RealK), INTENT(IN) ::                                           &
      d_planck_flux_surface(nd_profile)
!       Differential Planckian flux from the surface
  INTEGER, INTENT(IN) ::                                                &
      ls_brdf_trunc                                                     &
!       Order of trunation of BRDFs
    , n_brdf_basis_fnc
!       Number of BRDF basis functions
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb(nd_profile, nd_brdf_basis_fnc)                            &
!       Weights of the basis functions
    , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                         &
!       Array of BRDF basis terms
    , brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)             &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction

  TYPE(StrPlanck), INTENT(INOUT) :: planck
!       Planckian emission fields

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!       Spherical geometry fields
  
!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

!                 Cloudy properties
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud
!       Clouds required
  INTEGER, INTENT(IN) ::                                                &
      i_cloud
!       Cloud scheme used

!                 Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top                                                       &
!       Topmost cloudy layer
    , index_subcol                                                      &
!       Index of current sub-grid cloud column
    , n_region                                                          &
!       Number of cloudy regions
    , k_clr                                                             &
!       Index of clear-sky region
    , i_region_cloud(nd_cloud_type)
!       Regions in which types of clouds fall

! Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_column_slv(nd_profile)                                          &
!       Number of columns to be solved in each profile
    , list_column_slv(nd_profile, nd_column)                            &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(nd_profile, nd_column)                              &
!       Layer in the current column to change
    , i_clm_cld_typ(nd_profile, nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK), INTENT(IN) ::                                           &
      w_free(nd_profile, id_ct: nd_layer)                               &
!       Clear-sky fraction
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)                                             &
!       Coefficients for energy transfer at interfaces
    , area_column(nd_profile, nd_column)                                &
!       Areas of columns
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region


!                 Levels where radiance are calculated
  INTEGER, INTENT(IN) ::                                                &
      n_viewing_level                                                   &
!       Number of levels where radiances are calculated
    , i_rad_layer(nd_viewing_level)
!       Layers in which radiances are calculated
  REAL (RealK), INTENT(IN) ::                                           &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers

!                   Viewing Geometry
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)
!       Viewing directions

!                 Calculated Fluxes
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer, nd_k_term_inner)        &
!       Direct flux
    , flux_total(nd_flux_profile, 2*nd_layer+2, nd_k_term_inner)        &
!       Total flux
    , actinic_flux(nd_flux_profile, nd_layer, nd_k_term_inner)
!       Actinic flux
  LOGICAL, INTENT(IN) :: l_actinic
!       Flag for calculation of actinic flux

!                 Calculated radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      radiance(nd_radiance_profile, nd_viewing_level, nd_direction,     &
               nd_k_term_inner)
!       Radiances
!                 Calculated mean radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      j_radiance(nd_j_profile, nd_viewing_level, nd_k_term_inner)
!       Mean radiances

!                 Flags for clear-sky calculations
  LOGICAL, INTENT(IN) ::                                                &
      l_clear
!       Calculate clear-sky properties
  INTEGER, INTENT(IN) ::                                                &
      i_solver_clear
!       Clear solver used

!                 Clear-sky fluxes calculated
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct_clear(nd_flux_profile, 0: nd_layer, nd_k_term_inner)  &
!       Clear-sky direct flux
    , flux_total_clear(nd_flux_profile, 2*nd_layer+2, nd_k_term_inner)  &
!       Clear-sky total flux
    , actinic_flux_clear(nd_flux_profile, nd_layer, nd_k_term_inner)
!       Clear-sky actinic flux

  REAL (RealK), INTENT(INOUT) ::                                        &
      contrib_funci_part(nd_flux_profile, nd_layer, nd_k_term_inner)
!       Contribution function (intensity)
  REAL (RealK), INTENT(INOUT) ::                                        &
      contrib_funcf_part(nd_flux_profile, nd_layer, nd_k_term_inner)
!       Contribution function (flux)

! Local variables.
  INTEGER                                                               &
      k                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , i_n_term_inner

  REAL (RealK), ALLOCATABLE ::                                          &
      tau_clr_f(:, :)
!       Clear-sky optical depths for the whole column

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MONOCHROMATIC_RADIANCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO i_n_term_inner = 1 , n_k_term_inner

!   Calculate the optical depths and albedos of single scattering.
!   The phase function is not involved here as that is constant
!   across the band, whereas these parameters vary with the gaseous
!   absorption.


    CALL single_scattering_all(i_scatter_method                           &
!                 Atmospheric properties
      , n_profile, n_layer, d_mass                                        &
!                 Cloudy properties
      , l_cloud, n_cloud_top, cld%n_cloud_type                            &
!                 Optical properties
      , ss_prop, k_gas_abs(:,:,i_n_term_inner)                            &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct                         &
      )


    IF ( (i_angular_integration == ip_two_stream).OR.                     &
         (i_angular_integration == ip_spherical_harmonic) ) THEN

!     Rescale the optical depth and albedo of single scattering.
      IF (l_rescale) THEN

        !-----------------------------------------------------------------
        ! For direct solar flux, the optical depth can be scaled by three
        ! options: 
        ! ip_direct_csr_scaling - scaling optical depth with a forward
        ! scattering fraction with a FOV of pyrheliometer
        ! ip_direct_noscaling - no scaling 
        ! ip_direct_delta_scaling - Delta-Eddington scaling
        !-----------------------------------------------------------------
        IF (control%i_direct_tau == ip_direct_csr_scaling) THEN
          ! Rescale tau by CSR forward fraction
          ! Above cloud top.
          CALL rescale_tau_csr(n_profile                                  &
             , 1, n_cloud_top-1                                           &
             , ss_prop%forward_scatter_clr_csr                            &
             , ss_prop%tau_clr, ss_prop%tau_clr_dir                       &
             , ss_prop%omega_clr                                          &
             , nd_profile, nd_layer_clr, 1                                &
             )
          ! Below cloud top.
          CALL rescale_tau_csr(n_profile                                  &
             , n_cloud_top, n_layer                                       &
             , ss_prop%forward_scatter_csr(:, :, 0)                       &
             , ss_prop%tau(:, :, 0), ss_prop%tau_dir(:, :, 0)             &
             , ss_prop%omega(:, :, 0)                                     &
             , nd_profile, nd_layer, id_ct                                &
             )
        ELSE IF (control%i_direct_tau == ip_direct_noscaling) THEN
          DO i=1, n_cloud_top-1
            DO l=1, n_profile
              ss_prop%tau_clr_dir(l,i) = ss_prop%tau_clr(l,i)
            END DO
          END DO
          DO i=n_cloud_top, n_layer
            DO l=1, n_profile
              ss_prop%tau_dir(l,i,0)=ss_prop%tau(l,i,0)
            END DO
          END DO
        END IF

        CALL rescale_tau_omega(n_profile, 1, n_cloud_top-1                &
          , ss_prop%tau_clr, ss_prop%omega_clr                            &
          , ss_prop%forward_scatter_clr                                   &
          , nd_profile, nd_layer_clr, 1                                   &
          )
        CALL rescale_tau_omega(n_profile, n_cloud_top, n_layer            &
          , ss_prop%tau(:, :, 0), ss_prop%omega(:, :, 0)                  &
          , ss_prop%forward_scatter(:, :, 0)                              &
          , nd_profile, nd_layer, id_ct                                   &
          )

        IF (l_cloud) THEN

          IF (control%i_direct_tau == ip_direct_csr_scaling) THEN
            DO k=1, cld%n_cloud_type  
              CALL rescale_tau_csr(n_profile                              &
                 , n_cloud_top, n_layer                                   &
                 , ss_prop%forward_scatter_csr(:, :, k)                   &
                 , ss_prop%tau(:, :, k)                                   &
                 , ss_prop%tau_dir(:, :, k)                               &
                 , ss_prop%omega(:, :, k)                                 &
                 , nd_profile, nd_layer, id_ct                            &
                 )
            END DO
          ELSE IF (control%i_direct_tau == ip_direct_noscaling) THEN
            DO k=1, cld%n_cloud_type
              DO i=n_cloud_top, n_layer
                DO l=1, n_profile
                  ss_prop%tau_dir(l,i,k)=ss_prop%tau(l,i,k)
                END DO
              END DO
            END DO
          END IF
          DO k=1, cld%n_cloud_type
            CALL rescale_tau_omega(n_profile, n_cloud_top, n_layer        &
              , ss_prop%tau(:, :, k), ss_prop%omega(:, :, k)              &
              , ss_prop%forward_scatter(:, :, k)                          &
              , nd_profile, nd_layer, id_ct                               &
              )
          END DO

        END IF

      END IF

    END IF


    IF (control%l_spherical_solar) THEN
      CALL spherical_trans_coeff(n_profile, n_layer, n_cloud_top,         &
        ss_prop, sph)
    END IF


!   Now divide the algorithmic path depending on the option
!   for angular integration.

    IF (i_angular_integration == ip_two_stream) THEN

!     The standard two-stream approximations.
      CALL monochromatic_radiance_tseq(ierr                               &
        , control, cld, bound                                             &
!                   Atmospheric Propertries
        , n_profile, n_layer                                              &
!                   Options for Solver
        , i_2stream, i_solver, i_scatter_method                           &
!                   Optical Properties
        , l_scale_solar, adjust_solar_ke                                  &
!                   Spectral Region
        , isolir                                                          &
!                   Infra-red Properties
        , planck%diff, control%l_ir_source_quad, planck%diff_2            &
!                   Conditions at TOA
        , zen_0, flux_inc_direct, flux_inc_down                           &
!                   Surface Properties
        , d_planck_flux_surface                                           &
        , rho_alb                                                         &
!                   Spherical geometry
        , sph                                                             &
!                   Optical Properties
        , ss_prop                                                         &
!                   Cloudy Properties
        , i_cloud                                                         &
!                   Cloud Geometry
        , n_cloud_top, index_subcol                                       &
        , n_region, k_clr, i_region_cloud, frac_region                    &
        , w_free, cloud_overlap                                           &
        , n_column_slv, list_column_slv                                   &
        , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!                   Fluxes Calculated
        , flux_direct(:,:,i_n_term_inner)                                 &
        , flux_total(:,:,i_n_term_inner)                                  &
        , l_actinic                                                       &
        , actinic_flux(:,:,i_n_term_inner)                                &
!                   Flags for Clear-sky Calculation
        , l_clear, i_solver_clear                                         &
!                   Clear-sky Fluxes Calculated
        , flux_direct_clear (:,:,i_n_term_inner)                          &
        , flux_total_clear  (:,:,i_n_term_inner)                          &
        , actinic_flux_clear(:,:,i_n_term_inner)                          &
!                   Dimensions of Arrays
        , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column            &
        , nd_cloud_type, nd_region, nd_overlap_coeff                      &
        , nd_source_coeff, nd_max_order                                   &
        )


    ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN

!     The spherical harmonic option:
      CALL monochromatic_radiance_sph(ierr                                &
!                   Atmospheric Propertries
        , control, n_profile, n_layer                                     &
!                   Angular Integration
        , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc     &
        , accuracy_adaptive, euler_factor, i_sph_algorithm                &
        , i_sph_mode, l_rescale                                           &
!                   Precalculated angular arrays
        , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                        &
!                   Options for Equivalent Extinction
        , l_scale_solar, adjust_solar_ke                                  &
!                   Spectral Region
        , isolir                                                          &
!                   Infra-red Properties
        , planck%diff, control%l_ir_source_quad, planck%diff_2            &
!                   Conditions at TOA
        , zen_0, flux_inc_down                                            &
        , i_direct(:,:,i_n_term_inner)                                    &
!                   Surface Properties
        , d_planck_flux_surface                                           &
        , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                        &
        , f_brdf, brdf_sol, brdf_hemi                                     &
!                   Optical properties
        , ss_prop                                                         &
!                   Cloudy Properties
        , i_cloud                                                         &
!                   Cloud Geometry
        , n_cloud_top                                                     &
        , n_column_slv, list_column_slv                                   &
        , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!                   Levels for calculating radiances
        , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!                   Viewing geometry
        , n_direction, direction                                          &
!                   Calculated Fluxes
        , flux_direct(:,:,i_n_term_inner)                                 &
        , flux_total(:,:,i_n_term_inner)                                  &
!                   Calculated radiances
        , radiance(:,:,:,i_n_term_inner)                                  &
!                   Calculated mean radiances
        , j_radiance(:,:,i_n_term_inner)                                  &
!                   Dimensions of Arrays
        , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column            &
        , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
        , nd_max_order, nd_sph_coeff                                      &
        , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level              &
        , nd_direction                                                    &
        )

    ELSE IF (i_angular_integration == ip_ir_gauss) THEN

!     Full angular resolution using Gaussian integration.

      ALLOCATE(tau_clr_f(nd_profile, nd_layer))
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          tau_clr_f(l, i)=ss_prop%tau_clr(l, i)
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          tau_clr_f(l, i)=ss_prop%tau(l, i, 0)
        END DO
      END DO

      CALL gauss_angle(n_profile, n_layer                                 &
        , n_order_gauss                                                   &
        , tau_clr_f                                                       &
        , flux_inc_down                                                   &
        , planck%diff, d_planck_flux_surface                              &
        , rho_alb(1, ip_surf_alb_diff)                                    &
        , flux_total(:,:,i_n_term_inner)                                  &
        , nd_profile, nd_layer                                            &
        )

      DEALLOCATE(tau_clr_f)

    END IF

!   Calculate the contribution function
    IF (control%l_contrib_func .OR. control%l_contrib_func_band) THEN
      CALL calc_contrib_func(n_profile, n_layer, n_cloud_top              &
        , atm%p_level, planck%flux, ss_prop                               &
        , contrib_funci_part(:,:,i_n_term_inner)                          &
        , contrib_funcf_part(:,:,i_n_term_inner), nd_profile, nd_layer)
    END IF

  END DO  ! n_k_term_inner

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE monochromatic_radiance
END MODULE monochromatic_radiance_mod

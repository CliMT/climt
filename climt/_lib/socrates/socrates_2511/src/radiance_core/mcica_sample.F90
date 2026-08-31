! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to sample the cloud field for a single k-term
!
! Method:
!   Returns the fluxes for a single k-term. A number of cloudy
!   sub-columns are sampled (dependent on the configuration of McICA)
!   and a single clear-sky calculation is performed. The results are
!   meaned making use of the fraction of cloudy sub-columns.
!
!- ---------------------------------------------------------------------
MODULE mcica_sample_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MCICA_SAMPLE_MOD'
CONTAINS
SUBROUTINE mcica_sample(ierr                                            &
    , control, dimen, atm, cld, bound                                   &
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
    , k_gas_abs                                                         &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , planck                                                            &
!                 Conditions at TOA
    , zen_0, flux_inc_direct, flux_inc_down                             &
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
    , n_cloud_top                                                       &
    , n_region, k_clr, i_region_cloud, frac_region                      &
    , w_free, cloud_overlap                                             &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                 Additional variables required for mcica
    , l_cloud_cmp, n_cloud_profile, i_cloud_profile                     &
    , i_cloud_type, nd_cloud_component, iex, i_band                     &
    , i_cloud_representation                                            &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                 Viewing geometry
    , n_direction, direction                                            &
!                 Calculated Fluxes
    , flux_direct, flux_total, l_actinic, actinic_flux                  &
!                 Flags for Clear-sky Calculation
    , i_solver_clear                                                    &
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
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_dimen,   ONLY: StrDim
  USE def_atm,     ONLY: StrAtm
  USE def_cld,     ONLY: StrCld
  USE def_bound,   ONLY: StrBound
  USE def_planck,  ONLY: StrPlanck
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_solar, ip_cloud_homogen, ip_cloud_ice_water, &
    ip_cloud_combine_homogen, ip_cloud_combine_ice_water, &
    ip_no_scatter_abs, ip_no_scatter_ext, ip_direct_csr_scaling, &
    i_err_fatal
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE circumsolar_fraction_mod, ONLY: circumsolar_fraction
  USE monochromatic_radiance_mod, ONLY: monochromatic_radiance

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
  TYPE(StrDim),       INTENT(IN)    :: dimen

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
    , nd_source_coeff
!       Size allocated for source coefficients


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
  REAL (RealK), INTENT(IN) ::                                           &
      k_gas_abs(nd_profile, nd_layer)
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
  REAL (RealK), INTENT(INOUT) ::                                        &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct flux
    , flux_total(nd_flux_profile, 2*nd_layer+2)                         &
!       Total flux
    , actinic_flux(nd_flux_profile, nd_layer)
!       Actinic flux
  LOGICAL, INTENT(IN) :: l_actinic
!       Flag for calculation of actinic flux

!                 Flags for clear-sky calculations
  INTEGER, INTENT(IN) ::                                                &
      i_solver_clear
!       Clear solver used

!                 Clear-sky fluxes calculated
  REAL (RealK), INTENT(INOUT) ::                                        &
      flux_direct_clear(nd_flux_profile, 0: nd_layer)                   &
!       Clear-sky direct flux
    , flux_total_clear(nd_flux_profile, 2*nd_layer+2)                   &
!       Clear-sky total flux
    , actinic_flux_clear(nd_flux_profile, nd_layer)
!       Clear-sky actinic flux

! Variables required for McICA
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_profile(id_ct: nd_layer)                                  &
!       number of cloudy profiles in each layer
    , i_cloud_profile(nd_profile, id_ct: nd_layer)                      &
!       profiles containing clouds
    , nd_cloud_component                                                &
!       size allocated for components of clouds
    , i_cloud_type(nd_cloud_component)                                  &
!       types of cloud to which each component contributes
    , iex                                                               &
!       Index of ESFT term
    , i_band                                                            &
!       Band being considered
    , i_cloud_representation
!       Representation of mixing rule chosen

  LOGICAL, INTENT(IN) ::                                                &
      l_cloud_cmp(nd_cloud_component)
!       flags to activate cloudy components

  REAL (RealK), INTENT(INOUT) ::                                        &
    contrib_funci_part(nd_flux_profile, nd_layer)
!       Contribution function (intensity) increment
  REAL (RealK), INTENT(INOUT) ::                                        &
    contrib_funcf_part(nd_flux_profile, nd_layer)
!       Contribution function (flux) increment

! Local variables.
  INTEGER :: i, j, l, ll, ls, k, m
!       Loop variables

  INTEGER :: index_subcol
!       Index of current sub-grid cloud column

  INTEGER, PARAMETER ::                                                 &
      n_k_term_inner_dummy = 1                                          &
!       Number of monochromatic calculations in inner loop (dummy here)
    , nd_k_term_inner_dummy = 1
!       Maximum number of k-terms in inner loops (dummy here)

  REAL (RealK) ::                                                       &
      i_direct_subcol(nd_radiance_profile, 0: nd_layer)                 &
!       Partial solar irradiances
    , flux_direct_subcol(nd_flux_profile, 0: nd_layer)                  &
!       Partial direct flux
    , flux_total_subcol(nd_flux_profile, 2*nd_layer+2)                  &
!       Partial total flux
    , actinic_flux_subcol(nd_flux_profile, nd_layer)                    &
!       Partial actinic flux
    , flux_direct_sph(nd_flux_profile, 0: nd_layer+1)                   &
!       Direct flux through spherical geometry
    , flux_direct_div(nd_flux_profile, nd_layer)                        &
!       Direct flux divergence across layer
    , radiance_subcol(nd_radiance_profile, nd_viewing_level             &
        , nd_direction)                                                 &
!       Partial radiances
    , photolysis_subcol(nd_j_profile, nd_viewing_level)
!       Partial rate of photolysis

  REAL (RealK) ::                                                       &
      subcol_k_inv
!       1.0/cld%subcol_k

  LOGICAL ::                                                            &
      l_clear_calc
!       flag for calculating clear sky conrtribution to McICA

  REAL (RealK), PARAMETER :: eps=EPSILON(ss_prop%k_ext_scat)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER(LEN=*), PARAMETER :: RoutineName='MCICA_SAMPLE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO m=cld%first_subcol_k(i_band,iex), cld%first_subcol_k(i_band,iex+1)-1

    index_subcol=MOD(m, dimen%nd_subcol_req)
    IF (index_subcol == 0)  index_subcol = dimen%nd_subcol_req
    index_subcol=cld%subcol_reorder(index_subcol)
    subcol_k_inv=1.0e+00_RealK/cld%subcol_k(i_band,iex)


    IF (m==cld%first_subcol_k(i_band,iex))THEN
      l_clear_calc=.TRUE.
    ELSE
      l_clear_calc=.FALSE.
    END IF


!   Loop over the condensed components, calculating their optical
!   properties by simply scaling the values previously calculated
!   for the average cloud condensate amounts.
    SELECT CASE (i_cloud_representation)
    CASE (ip_cloud_ice_water, ip_cloud_combine_ice_water)
      IF (control%l_avg_phase_fnc .OR.                                  &
            (i_scatter_method == ip_no_scatter_abs) .OR.                &
            (i_scatter_method == ip_no_scatter_ext) ) THEN
        DO k=1, cld%n_condensed
          IF (l_cloud_cmp(k)) THEN
            DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                l=i_cloud_profile(ll, i)
                ss_prop%k_grey_tot(l, i, i_cloud_type(k))               &
                  =ss_prop%k_grey_tot(l, i, 0)                          &
                  +ss_prop%k_ext_tot_cloud_comp(l, i, k)                &
                  *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                ss_prop%k_ext_scat(l, i, i_cloud_type(k))               &
                  =ss_prop%k_ext_scat(l, i, 0)                          &
                  +ss_prop%k_ext_scat_cloud_comp(l, i, k)               &
                  *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
              END DO
            END DO
          END IF
        END DO
      ELSE
        IF (l_rescale) THEN
          DO k=1, cld%n_condensed
            IF (l_cloud_cmp(k)) THEN
              DO i=n_cloud_top, n_layer
!CDIR NODEP
                DO ll=1, n_cloud_profile(i)
                  l=i_cloud_profile(ll, i)
                  ss_prop%k_grey_tot(l, i, i_cloud_type(k))             &
                    =ss_prop%k_grey_tot(l, i, 0)                        &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  ss_prop%k_ext_scat(l, i, i_cloud_type(k))             &
                    =ss_prop%k_ext_scat(l, i, 0)                        &
                    +ss_prop%k_ext_scat_cloud_comp(l, i, k)             &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  ss_prop%forward_scatter(l, i, i_cloud_type(k))        &
                    =ss_prop%forward_scatter_no_cloud(l, i)             &
                    +ss_prop%forward_scatter_cloud_comp(l, i, k)        &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  DO ls=1, n_order_phase
                    ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))        &
                      =(ss_prop%phase_fnc_no_cloud(l, i, ls)            &
                      +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)        &
                      *cld%c_sub(l, i, index_subcol, i_cloud_type(k))   &
                      -ss_prop%forward_scatter(l, i, i_cloud_type(k)))  &
                      /MAX(ss_prop%k_ext_scat(l, i, i_cloud_type(k))    &
                      -ss_prop%forward_scatter(l, i, i_cloud_type(k))   &
                      ,eps)
                  END DO
                  ss_prop%forward_scatter(l, i, i_cloud_type(k))        &
                    =ss_prop%forward_scatter(l, i, i_cloud_type(k))     &
                    /MAX(ss_prop%k_ext_scat(l, i, i_cloud_type(k))      &
                    ,eps)
                END DO
              END DO
            END IF
          END DO
        ELSE
          DO k=1, cld%n_condensed
            IF (l_cloud_cmp(k)) THEN
              DO i=n_cloud_top, n_layer
!CDIR NODEP
                DO ll=1, n_cloud_profile(i)
                  l=i_cloud_profile(ll, i)
                  ss_prop%k_grey_tot(l, i, i_cloud_type(k))             &
                    =ss_prop%k_grey_tot(l, i, 0)                        &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  ss_prop%k_ext_scat(l, i, i_cloud_type(k))             &
                    =ss_prop%k_ext_scat(l, i, 0)                        &
                    +ss_prop%k_ext_scat_cloud_comp(l, i, k)             &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  DO ls=1, n_order_phase
                    ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))        &
                      =(ss_prop%phase_fnc_no_cloud(l, i, ls)            &
                      +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)        &
                      *cld%c_sub(l, i, index_subcol, i_cloud_type(k)))  &
                      /MAX(ss_prop%k_ext_scat(l, i, i_cloud_type(k))    &
                      ,eps)
                  END DO
                END DO
              END DO
            END IF
          END DO
        END IF
      END IF

    CASE (ip_cloud_homogen, ip_cloud_combine_homogen)
! Since both ice and liquid are added to same i_cloud_type, the clear-sky
! extinction, cloudy extinction and scaling of phase function must each be
! done in separate loops.

! Addition of clear-sky extinction
      DO k=1, cld%n_cloud_type
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%k_grey_tot(l, i, k)=ss_prop%k_grey_tot(l, i, 0)
            ss_prop%k_ext_scat(l, i, k)=ss_prop%k_ext_scat(l, i, 0)
          END DO
        END DO
      END DO
      IF (.NOT. control%l_avg_phase_fnc) THEN
        DO k=1, cld%n_cloud_type
          DO i=n_cloud_top, n_layer
            DO ls=1, n_order_phase
              DO l=1, n_profile
                ss_prop%phase_fnc(l, i, ls, k)                          &
                  =ss_prop%phase_fnc_no_cloud(l, i, ls)
              END DO
            END DO
          END DO
        END DO
        IF (l_rescale) THEN
          DO k=1, cld%n_cloud_type
            DO i=n_cloud_top, n_layer
              DO l=1, n_profile
                ss_prop%forward_scatter(l, i, k)                        &
                  =ss_prop%forward_scatter_no_cloud(l, i)
              END DO
            END DO
          END DO
        END IF
      END IF

! Add on the terms representing processes within clouds.
      IF (control%l_avg_phase_fnc) THEN
        DO k=1, cld%n_condensed
          IF (l_cloud_cmp(k)) THEN
            DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                l=i_cloud_profile(ll, i)
                ss_prop%k_grey_tot(l, i, i_cloud_type(k))               &
                  =ss_prop%k_grey_tot(l, i, i_cloud_type(k))            &
                  +ss_prop%k_ext_tot_cloud_comp(l, i, k)                &
                  *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                ss_prop%k_ext_scat(l, i, i_cloud_type(k))               &
                  =ss_prop%k_ext_scat(l, i, i_cloud_type(k))            &
                  +ss_prop%k_ext_scat_cloud_comp(l, i, k)               &
                  *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
              END DO
            END DO
          END IF
        END DO
      ELSE
        IF (l_rescale) THEN
          DO k=1, cld%n_condensed
            IF (l_cloud_cmp(k)) THEN
              DO i=n_cloud_top, n_layer
!CDIR NODEP
                DO ll=1, n_cloud_profile(i)
                  l=i_cloud_profile(ll, i)
                  ss_prop%k_grey_tot(l, i, i_cloud_type(k))             &
                    =ss_prop%k_grey_tot(l, i, i_cloud_type(k))          &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  ss_prop%k_ext_scat(l, i, i_cloud_type(k))             &
                    =ss_prop%k_ext_scat(l, i, i_cloud_type(k))          &
                    +ss_prop%k_ext_scat_cloud_comp(l, i, k)             &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  ss_prop%forward_scatter(l, i, i_cloud_type(k))        &
                    =ss_prop%forward_scatter(l, i, i_cloud_type(k))     &
                    +ss_prop%forward_scatter_cloud_comp(l, i, k)        &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  DO ls=1, n_order_phase
                    ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))        &
                      =ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))     &
                      +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)        &
                      *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  END DO
                END DO
              END DO
            END IF
          END DO
        ELSE
          DO k=1, cld%n_condensed
            IF (l_cloud_cmp(k)) THEN
              DO i=n_cloud_top, n_layer
!CDIR NODEP
                DO ll=1, n_cloud_profile(i)
                  l=i_cloud_profile(ll, i)
                  ss_prop%k_grey_tot(l, i, i_cloud_type(k))             &
                    =ss_prop%k_grey_tot(l, i, i_cloud_type(k))          &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  ss_prop%k_ext_scat(l, i, i_cloud_type(k))             &
                    =ss_prop%k_ext_scat(l, i, i_cloud_type(k))          &
                    +ss_prop%k_ext_scat_cloud_comp(l, i, k)             &
                    *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  DO ls=1, n_order_phase
                    ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))        &
                      =ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))     &
                      +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)        &
                      *cld%c_sub(l, i, index_subcol, i_cloud_type(k))
                  END DO
                END DO
              END DO
            END IF
          END DO
        END IF
      END IF

! Calculate the final optical properties.
! The scattering was included in the free total extinction earlier,
! but we have yet to divide the product of the phase function and
! the scattering by the mean scattering.
      IF (.NOT. control%l_avg_phase_fnc) THEN
        IF (l_rescale) THEN
          DO k=1, cld%n_cloud_type
            DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                l=i_cloud_profile(ll, i)
                DO ls=1, n_order_phase
                  ss_prop%phase_fnc(l, i, ls, k)                        &
                    =(ss_prop%phase_fnc(l, i, ls, k)                    &
                    -ss_prop%forward_scatter(l, i, k))                  &
                    /MAX(ss_prop%k_ext_scat(l, i, k)                    &
                    -ss_prop%forward_scatter(l, i, k),eps)
                END DO
                ss_prop%forward_scatter(l, i, k)                        &
                  =ss_prop%forward_scatter(l, i, k)                     &
                  /MAX(ss_prop%k_ext_scat(l, i,k)                       &
                    ,eps)
              END DO

              IF (control%i_direct_tau == ip_direct_csr_scaling ) THEN
! Calculate forward scattering fraction of direct flux within 
! the instrument FOV 
                 CALL circumsolar_fraction(n_cloud_profile(i)           &
                  , i_cloud_profile(:, i), control%half_angle           &
                  , ss_prop%phase_fnc(:, i, 1, k)                       &
                  , ss_prop%forward_scatter_csr(:, i, k)                &
                  , nd_profile                                          &
                  ) 
              END IF
            END DO
          END DO
        ELSE
          DO k=1, cld%n_cloud_type
            DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                l=i_cloud_profile(ll, i)
                DO ls=1, n_order_phase
                  ss_prop%phase_fnc(l, i, ls, k)                        &
                    =ss_prop%phase_fnc(l, i, ls, k)                     &
                    /MAX(ss_prop%k_ext_scat(l, i, k),eps)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF

    CASE DEFAULT

      cmessage = 'Cloud representation not compatible with MCICA'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)

    END SELECT


    CALL monochromatic_radiance(ierr                                    &
      , control, atm, cld, bound                                        &
!             atmospheric properties
      , n_profile, n_layer, d_mass                                      &
!             angular integration
      , i_angular_integration, i_2stream                                &
      , l_rescale, n_order_gauss                                        &
      , n_order_phase, ms_min, ms_max, i_truncation                     &
      , ls_local_trunc                                                  &
      , accuracy_adaptive, euler_factor                                 &
      , i_sph_algorithm, i_sph_mode                                     &
!             precalculated angular arrays
      , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                        &
!             treatment of scattering
      , i_scatter_method                                                &
!             options for solver
      , i_solver                                                        &
!             gaseous propreties
      , n_k_term_inner_dummy, k_gas_abs                                 &
!             options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!             spectral region
      , isolir                                                          &
!             infra-red properties
      , planck                                                          &
!             conditions at toa
      , zen_0, flux_inc_direct, flux_inc_down                           &
      , i_direct_subcol                                                 &
!             surface properties
      , d_planck_flux_surface                                           &
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                        &
      , f_brdf, brdf_sol, brdf_hemi                                     &
!             Spherical geometry
      , sph                                                             &
!             optical properties
      , ss_prop                                                         &
!             cloudy properties
      , l_cloud, i_cloud                                                &
!             cloud geometry
      , n_cloud_top, index_subcol                                       &
      , n_region, k_clr, i_region_cloud, frac_region                    &
      , w_free, cloud_overlap                                           &
      , n_column_slv, list_column_slv                                   &
      , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!               levels for the calculation of radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!               viewing geometry
      , n_direction, direction                                          &
!               calculated fluxes
      , flux_direct_subcol, flux_total_subcol                           &
      , l_actinic, actinic_flux_subcol                                  &
!               calculated radiances
      , radiance_subcol                                                 &
!               calculated rate of photolysis
      , photolysis_subcol                                               &
!             flags for clear-sky calculations
      , l_clear_calc, i_solver_clear                                    &
!             clear-sky fluxes calculated
      , flux_direct_clear, flux_total_clear, actinic_flux_clear         &
!             contribution function
      , contrib_funci_part, contrib_funcf_part                          &
!             dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column            &
      , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
      , nd_cloud_type, nd_region, nd_overlap_coeff                      &
      , nd_max_order, nd_sph_coeff                                      &
      , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level              &
      , nd_direction, nd_source_coeff, nd_k_term_inner_dummy            &
      )

    IF (m == cld%first_subcol_k(i_band,iex)) THEN

      IF (isolir == ip_solar) THEN
        DO k=0,n_layer
          DO j=1,n_profile
            flux_direct(j,k)=flux_direct_subcol(j,k)
          END DO
        END DO
        IF (control%l_spherical_solar) THEN
          DO k=0,n_layer+1
            DO j=1,n_profile
              flux_direct_sph(j,k)=sph%allsky%flux_direct(j,k)
            END DO
          END DO
          DO k=1,n_layer
            DO j=1,n_profile
              flux_direct_div(j,k)=sph%allsky%flux_direct_div(j,k)
            END DO
          END DO
        END IF
      END IF

      DO k=1,2*n_layer+2
        DO j=1,n_profile
          flux_total(j,k)=flux_total_subcol(j,k)
        END DO
      END DO

      IF (l_actinic) THEN
        DO k=1,n_layer
          DO j=1,n_profile
            actinic_flux(j,k)=actinic_flux_subcol(j,k)
          END DO
        END DO
      END IF

    ELSE

      IF (isolir == ip_solar) THEN
        DO k=0,n_layer
          DO j=1,n_profile
            flux_direct(j,k)=flux_direct(j,k)+flux_direct_subcol(j,k)
          END DO
        END DO
        IF (control%l_spherical_solar) THEN
          DO k=0,n_layer+1
            DO j=1,n_profile
              flux_direct_sph(j,k)=flux_direct_sph(j,k) &
                + sph%allsky%flux_direct(j,k)
            END DO
          END DO
          DO k=1,n_layer
            DO j=1,n_profile
              flux_direct_div(j,k)=flux_direct_div(j,k) &
                + sph%allsky%flux_direct_div(j,k)
            END DO
          END DO
        END IF
      END IF

      DO k=1,2*n_layer+2
        DO j=1,n_profile
          flux_total(j,k)=flux_total(j,k)+flux_total_subcol(j,k)
        END DO
      END DO

      IF (l_actinic) THEN
        DO k=1,n_layer
          DO j=1,n_profile
            actinic_flux(j,k)=actinic_flux(j,k)+actinic_flux_subcol(j,k)
          END DO
        END DO
      END IF

    END IF

  END DO


  IF (isolir == ip_solar) THEN
    DO k=0,n_layer
      DO j=1,n_profile
        flux_direct(j,k)=flux_direct(j,k)*subcol_k_inv
        flux_direct(j,k)=cld%frac_cloudy(j)*flux_direct(j,k) &
          + ((1.0_RealK-cld%frac_cloudy(j))*flux_direct_clear(j,k))
      END DO
    END DO
    IF (control%l_spherical_solar) THEN
      DO k=0,n_layer+1
        DO j=1,n_profile
          sph%allsky%flux_direct(j,k) &
            = cld%frac_cloudy(j)*flux_direct_sph(j,k)*subcol_k_inv &
            + ((1.0_RealK-cld%frac_cloudy(j))*sph%clear%flux_direct(j,k))
        END DO
      END DO
      DO k=1,n_layer
        DO j=1,n_profile
          sph%allsky%flux_direct_div(j,k) &
            = cld%frac_cloudy(j)*flux_direct_div(j,k)*subcol_k_inv &
            + ((1.0_RealK-cld%frac_cloudy(j))*sph%clear%flux_direct_div(j,k))
        END DO
      END DO
    END IF
  END IF

  DO k=1,2*n_layer+2
    DO j=1,n_profile
      flux_total(j,k)=flux_total(j,k)*subcol_k_inv
      flux_total(j,k)=(cld%frac_cloudy(j)*flux_total(j,k)) &
        + ((1.0_RealK-cld%frac_cloudy(j))*flux_total_clear(j,k))
    END DO
  END DO

  IF (l_actinic) THEN
    DO k=1,n_layer
      DO j=1,n_profile
        actinic_flux(j,k) = actinic_flux(j,k) * subcol_k_inv
        actinic_flux(j,k) = cld%frac_cloudy(j) * actinic_flux(j,k) &
          + ( 1.0_RealK - cld%frac_cloudy(j) ) * actinic_flux_clear(j,k)
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE mcica_sample
END MODULE mcica_sample_mod

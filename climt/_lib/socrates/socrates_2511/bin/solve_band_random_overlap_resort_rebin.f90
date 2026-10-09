! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the fluxes assuming random overlap using
! resorting and rebinning of the k-terms to reduce computational cost.
!
! Method:
!   Combinations of k-terms for two gases are calculated, reordered
!   and rebinned. This is repeated for each gas in the band. Monochromatic
!   calculations are performed for each of the resampled k-terms and the
!   results are summed.
!
!- ---------------------------------------------------------------------
MODULE solve_band_random_overlap_resort_rebin_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'SOLVE_BAND_RANDOM_OVERLAP_RESORT_REBIN_MOD'
CONTAINS
SUBROUTINE solve_band_random_overlap_resort_rebin(ierr                         &
    , control, dimen, spectrum, atm, cld, bound, radout, i_band                &
!                 Atmospheric Column
    , n_profile, n_layer, d_mass                                               &
!                 Angular Integration
    , i_angular_integration, i_2stream                                         &
    , n_order_phase, l_rescale, n_order_gauss                                  &
    , ms_min, ms_max, i_truncation, ls_local_trunc                             &
    , accuracy_adaptive, euler_factor                                          &
    , i_sph_algorithm, i_sph_mode                                              &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                                 &
!                 Treatment of Scattering
    , i_scatter_method                                                         &
!                 Options for solver
    , i_solver                                                                 &
    , n_esft_red_in, gpnt_split                                                &
!                 Gaseous Properties
    , n_abs, index_abs, n_abs_esft                                             &
    , k_abs_layer, w_abs_esft                                                  &
!                 Spectral Region
    , isolir                                                                   &
!                 Solar Properties
    , zen_0, solar_irrad, sph                                                  &
!                 Infra-red Properties
    , planck                                                                   &
!                 Surface Properties
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                                 &
    , f_brdf, brdf_sol, brdf_hemi                                              &
!                 Tiling of the surface
    , l_tile, n_point_tile, n_tile, list_tile, rho_alb_tile                    &
!                 Optical Properties
    , ss_prop, photol                                                          &
!                 Cloudy Properties
    , l_cloud, i_cloud                                                         &
!                 Cloud Geometry
    , n_cloud_top                                                              &
    , n_region, k_clr, i_region_cloud, frac_region                             &
    , w_free, cloud_overlap                                                    &
    , n_column_slv, list_column_slv                                            &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                                &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                             &
!                 Viewing Geometry
    , n_direction, direction                                                   &
!                 Calculcated radiances
    , i_direct                                                                 &
!                 Flags for initialising diagnostics
    , l_initial, l_initial_band                                                &
    , l_initial_channel, l_initial_channel_tile                                &
!                 Flags for Fluxes
    , l_actinic, l_clear, i_solver_clear                                       &
!                 Dimensions
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column                     &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                       &
    , nd_abs                                                                   &
    , nd_esft_term, nd_esft_max                                                &
    , nd_cloud_type, nd_region, nd_overlap_coeff                               &
    , nd_max_order, nd_sph_coeff                                               &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                       &
    , nd_direction, nd_source_coeff                                            &
    , nd_point_tile, nd_tile                                                   &
    )


  USE realtype_rd,  ONLY: RealK
  USE def_control,  ONLY: StrCtrl
  USE def_dimen,    ONLY: StrDim
  USE def_spectrum, ONLY: StrSpecData
  USE def_atm,      ONLY: StrAtm
  USE def_cld,      ONLY: StrCld
  USE def_bound,    ONLY: StrBound
  USE def_out,      ONLY: StrOut
  USE def_planck,   ONLY: StrPlanck
  USE def_ss_prop,  ONLY: str_ss_prop
  USE def_qy,       ONLY: StrQy
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_solar, ip_infra_red, ip_spherical_harmonic,            &
                     ip_two_stream, ip_ir_gauss
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE augment_radiance_mod, ONLY: augment_radiance
  USE augment_tiled_radiance_mod, ONLY: augment_tiled_radiance
  USE calc_gauss_weight_90_mod, ONLY: calc_gauss_weight_90
  USE monochromatic_radiance_mod, ONLY: monochromatic_radiance
  USE quicksort_mod, ONLY: quicksort
  USE rebin_esft_terms_mod, ONLY: rebin_esft_terms

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
  TYPE(StrDim),       INTENT(IN)    :: dimen

! Spectral data:
  TYPE(StrSpecData),  INTENT(IN)    :: spectrum

! Atmospheric properties:
  TYPE(StrAtm),       INTENT(IN)    :: atm

! Cloud properties:
  TYPE(StrCld),       INTENT(IN)    :: cld

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Output fields:
  TYPE(StrOut),       INTENT(INOUT) :: radout

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                       &
      nd_profile                                                               &
!       Maximum number of profiles
    , nd_layer                                                                 &
!       Maximum number of layers
    , nd_layer_clr                                                             &
!       Size allocated for totally clear layers
    , id_ct                                                                    &
!       Topmost declared cloudy layer
    , nd_flux_profile                                                          &
!       Size allocated for profiles in arrays of fluxes
    , nd_radiance_profile                                                      &
!       Size allocated for profiles in arrays of radiances
    , nd_j_profile                                                             &
!       Size allocated for profiles in arrays of mean radiances
    , nd_abs                                                                   &
!       Maximum number of absorbers
    , nd_esft_term                                                             &
!       Maximum number of ESFT terms for each gas
    , nd_esft_max                                                              &
!       Maximum number of ESFT terms needed in arrays
    , nd_column                                                                &
!       Number of columns per point
    , nd_cloud_type                                                            &
!       Size allocated for cloud types
    , nd_region                                                                &
!       Size allocated for cloudy regions
    , nd_overlap_coeff                                                         &
!       Size allocated for cloudy overlap coefficients
    , nd_max_order                                                             &
!       Size allocated for orders of spherical harmonics
    , nd_sph_coeff                                                             &
!       Size allocated for spherical harmonic coefficients
    , nd_brdf_basis_fnc                                                        &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                            &
!       Size allowed for orders of BRDFs
    , nd_viewing_level                                                         &
!       Size allocated for levels where radiances are calculated
    , nd_direction                                                             &
!       Size allocated for viewing directions
    , nd_source_coeff                                                          &
!       Size allocated for source coefficients
    , nd_point_tile                                                            &
!       Size allocated for points where the surface is tiled
    , nd_tile
!       Size allocated for surface tiles


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                                    &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                       &
      i_band
!       Band being considered

!                 Atmospheric column
  INTEGER, INTENT(IN) ::                                                       &
      n_profile                                                                &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                                  &
      d_mass(nd_profile, nd_layer)
!       Mass thickness of each layer

!                 Angular integration
  INTEGER, INTENT(IN) ::                                                       &
      i_angular_integration                                                    &
!       Angular integration scheme
    , i_2stream                                                                &
!       Two-stream scheme
    , n_order_phase                                                            &
!       Maximum order of terms in the phase function used in
!       the direct calculation of spherical harmonics
    , n_order_gauss                                                            &
!       Order of gaussian integration
    , ms_min                                                                   &
!       Lowest azimuthal order used
    , ms_max                                                                   &
!       Highest azimuthal order used
    , i_truncation                                                             &
!       Type of spherical truncation used
    , ia_sph_mm(0: nd_max_order)                                               &
!       Address of spherical coefficient of (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                          &
!       Orders of truncation at each azimuthal order
    , i_sph_mode                                                               &
!       Mode in which the spherical solver is to be used
    , i_sph_algorithm
!       Algorithm used for spherical harmonic calculation
  LOGICAL, INTENT(IN) ::                                                       &
      l_rescale
!       Rescale optical properties
  REAL (RealK) ::                                                              &
      cg_coeff(nd_sph_coeff)                                                   &
!       Clebsch-Gordan coefficients
    , uplm_zero(nd_sph_coeff)                                                  &
!       Values of spherical harmonics at polar angles pi/2
    , uplm_sol(nd_radiance_profile, nd_sph_coeff)                              &
!       Values of spherical harmonics in the solar direction
    , accuracy_adaptive                                                        &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series

!                 Treatment of scattering
  INTEGER, INTENT(IN) ::                                                       &
      i_scatter_method
!       Method of treating scattering

!                 Options for solver
  INTEGER, INTENT(IN) ::                                                       &
      i_solver                                                                 &
!       Solver used
    , n_esft_red_in
!       Number of reduced terms
  REAL (RealK), INTENT(IN) ::                                                  &
      gpnt_split
!       g-coordinate for splitting g-space

!                 Gaseous properties
  INTEGER, INTENT(IN) ::                                                       &
      n_abs                                                                    &
!       Number of gases in band
    , index_abs(nd_abs)                                                        &
!       Local indexing numbers for gases and continua
    , n_abs_esft(nd_abs)
!       Number of terms in band
  REAL (RealK), INTENT(IN) ::                                                  &
      k_abs_layer(nd_profile, nd_layer, nd_esft_term, nd_abs)                  &
!       Exponential ESFT terms at actual pressure layer
    , w_abs_esft(nd_esft_term, nd_abs)
!       Weights for ESFT

!                 Spectral region
  INTEGER, INTENT(IN) ::                                                       &
      isolir
!       Spectral region

!                 Solar properties
  REAL (RealK), INTENT(IN) ::                                                  &
      zen_0(nd_profile)                                                        &
!       Secants (two-stream) or cosines (spherical harmonics)
!       of the solar zenith angle
    , solar_irrad(nd_profile)
!       Incident solar irradiance in band

  TYPE(StrPlanck), INTENT(INOUT) :: planck
!   Planckian emission fields

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!   Spherical geometry fields

!                 Surface properties
  INTEGER, INTENT(IN) ::                                                       &
      ls_brdf_trunc                                                            &
!       Order of truncation of BRDFs
    , n_brdf_basis_fnc
!       Number of BRDF basis functions
  REAL (RealK), INTENT(IN) ::                                                  &
      rho_alb(nd_profile, nd_brdf_basis_fnc)                                   &
!       Weights of the basis functions
    , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                             &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                                &
!       Array of BRDF basis terms
    , brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)                    &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction

! Variables related to tiling of the surface
  LOGICAL, INTENT(IN) ::                                                       &
      l_tile
!       Logical to allow invoke options
  INTEGER, INTENT(IN) ::                                                       &
      n_point_tile                                                             &
!       Number of points to tile
    , n_tile                                                                   &
!       Number of tiles used
    , list_tile(nd_point_tile)
!       List of points with surface tiling
  REAL (RealK), INTENT(IN) ::                                                  &
      rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc, nd_tile)
!       Weights for the basis functions of the BRDFs
!       at the tiled points

!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

  TYPE(StrQy), INTENT(IN) :: photol(spectrum%photol%n_pathway)
!   Photolysis quantum yields interpolated to model grid temperatures

!                 Cloudy properties
  LOGICAL, INTENT(IN) ::                                                       &
      l_cloud
!       Cloud enabled
  INTEGER, INTENT(IN) ::                                                       &
      i_cloud
!       Cloud scheme used

!                 Cloud geometry
  INTEGER, INTENT(IN) ::                                                       &
      n_cloud_top                                                              &
!       Topmost cloudy layer
    , n_region                                                                 &
!       Number of cloudy regions
    , k_clr                                                                    &
!       Index of clear-sky region
    , i_region_cloud(nd_cloud_type)
!       Regions in which types of clouds fall

! Cloud geometry
  INTEGER, INTENT(IN) ::                                                       &
      n_column_slv(nd_profile)                                                 &
!       Number of columns to be solved in each profile
    , list_column_slv(nd_profile, nd_column)                                   &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(nd_profile, nd_column)                                     &
!       Layer in the current column to change
    , i_clm_cld_typ(nd_profile, nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK), INTENT(IN) ::                                                  &
      w_free(nd_profile, id_ct: nd_layer)                                      &
!       Clear-sky fraction
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                              &
        , nd_overlap_coeff)                                                    &
!       Coefficients for transfer for energy at interfaces
    , area_column(nd_profile, nd_column)                                       &
!       Areas of columns
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region



!                   Viewing Geometry
  INTEGER, INTENT(IN) ::                                                       &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                                  &
      direction(nd_radiance_profile, nd_direction, 2)
!       Viewing directions
  INTEGER, INTENT(IN) ::                                                       &
      n_viewing_level                                                          &
!       Number of levels where radiances are calculated
    , i_rad_layer(nd_viewing_level)
!       Layers in which radiances are calculated
  REAL (RealK), INTENT(IN) ::                                                  &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers

!                   Flags for initialising diagnostics
  LOGICAL, INTENT(INOUT) ::                                                    &
      l_initial                                                                &
!       Initialise rather than increment broadband diagnostics
    , l_initial_band(spectrum%dim%nd_band)                                     &
!       Initialise rather than increment band-by-band diagnostics
    , l_initial_channel(dimen%nd_channel)                                      &
!       Initialise rather than increment channel diagnostics
    , l_initial_channel_tile(dimen%nd_channel)
!       Initialise rather than increment channel diagnostics on tiles
    
!                   Flags for flux calculations
  LOGICAL, INTENT(IN) ::                                                       &
      l_actinic
!       Flag for calculation of actinic flux
  LOGICAL, INTENT(IN) ::                                                       &
      l_clear
!       Calculate clear-sky properties
  INTEGER, INTENT(IN) ::                                                       &
      i_solver_clear
!       Clear solver used

!                   Calculated radiances
  REAL (RealK), INTENT(INOUT) ::                                               &
      i_direct(nd_radiance_profile, 0: nd_layer)
!       Direct solar irradiance on levels



! Local variables.
  INTEGER ::                                                                   &
      k                                                                        &
!       Loop variable
    , l                                                                        &
!       Loop variable
    , i_profile                                                                &
!       Loop variable
    , i_layer
!       Loop variable
  INTEGER ::                                                                   &
      i_abs_band                                                               &
!       Index of active gas
    , iex                                                                      &
!       Index of ESFT term
    , iex_minor(nd_abs)                                                        &
!       Index of ESFT term for minor gas (dummy here)
    , i_band_esft_mix                                                          &
!       Number of terms in band for two absorbers combined
    , i_band_esft_mix_red                                                      &
!       Current number of terms in band for two absorbers combined after
!       resorting and rebinning
    , iex_mix                                                                  &
!       Index of term for two absorbers combined
    , iex_mix_red                                                              &
!       Index of term for two absorbers combined after resorting and
!       rebinning
    , n_esft_red
!       Number of reduced terms
  INTEGER, PARAMETER ::                                                        &
      n_k_term_inner_dummy = 1                                                 &
!       Number of monochromatic calculations in inner loop (dummy here)
    , nd_k_term_inner_dummy = 1
!       Maximum number of k-terms in inner loops (dummy here)
  REAL (RealK) ::                                                              &
      k_gas_abs(nd_profile, nd_layer)                                          &
!       Gaseous absorption
    , d_planck_flux_surface(nd_profile)                                        &
!       Difference in Planckian fluxes between the surface and
!       the air
    , flux_inc_direct(nd_profile)                                              &
!       Incident direct flux
    , flux_inc_down(nd_profile)                                                &
!       Incident downward flux
    , dummy_ke(nd_profile, nd_layer)                                           &
!       Dummy array (not used)
    , k_esft_layer_mix(nd_profile, nd_layer,                                   &
        nd_esft_max*nd_esft_max)                                               &
!       ESFT monochromatic exponents for the mixture of two absorbers
    , w_esft_mix(nd_esft_max*nd_esft_max)                                      &
!       Weights for the mixture of two absorbers
    , w_esft_mix_unsort(nd_esft_max*nd_esft_max)                               &
!       Unsorted weights for the mixture of two absorbers
    , k_esft_layer_mix_red(nd_profile, nd_layer, nd_esft_max)                  &
!       Reduced (resorted and rebinned) ESFT monochromatic exponents
!       for the mixture of two absorbers
    , w_esft_red(nd_esft_max+1)                                                &
!       Current weights for the mixture of two absorbers
    , glim_mix_red(nd_esft_max+1)                                              &
!       g-coordinate limits for the reduced terms
    , w_esft_target(nd_esft_max)                                               &
!       Target weights when performing rebinning
    , glim_target(nd_esft_max+1)                                               &
!       Target g-coordinate limits for the terms when performing
!       rebinning
    , w_esft_gauleg(nd_esft_max)                                               &
!       Weights using Gauss-Legendre quadrature
    , gpnt_gauleg(nd_esft_max+1)                                               &
!       g-coordinate points using Gauss-Legendre quadrature
    , xl
!       Variable used for transforming Gauss-Legendre quadrature points
!       and weights from [-1,1] to [a,b].

! Monochromatic incrementing radiances:
  REAL (RealK) ::                                                              &
      flux_direct_part(nd_flux_profile, 0: nd_layer)                           &
!       Partial direct flux
    , flux_direct_ground_part(nd_flux_profile)                                 &
!       Partial direct flux at the surface
    , flux_total_part(nd_flux_profile, 2*nd_layer+2)                           &
!       Partial total flux
    , actinic_flux_part(nd_flux_profile, nd_layer)                             &
!       Partial actinic flux
    , flux_direct_clear_part(nd_flux_profile, 0: nd_layer)                     &
!       Partial clear-sky direct flux
    , flux_total_clear_part(nd_flux_profile, 2*nd_layer+2)                     &
!       Partial clear-sky total flux
    , actinic_flux_clear_part(nd_flux_profile, nd_layer)
!       Clear partial actinic flux
  REAL (RealK) ::                                                              &
      i_direct_part(nd_radiance_profile, 0: nd_layer)                          &
!       Partial solar irradiances
    , radiance_part(nd_radiance_profile, nd_viewing_level                      &
        , nd_direction)
!       Partial radiances
  REAL (RealK) ::                                                              &
      photolysis_part(nd_j_profile, nd_viewing_level)
!       Partial rates of photolysis
  REAL (RealK) ::                                                              &
      weight_incr                                                              &
!       Weight applied to increments
    , weight_blue_incr                                                         &
!       Weight applied to blue increments
    , weight_sub_band_incr
!       Weight applied to sub-band increments

  INTEGER :: i_abs_major
!       Index for the major gas (first gas in band)
  
  REAL (RealK) ::                                                       &
      contrib_funci_part(nd_flux_profile, nd_layer)
!       Contribution (or weighting) function
  REAL (RealK) ::                                                       &
      contrib_funcf_part(nd_flux_profile, nd_layer)
!       Contribution (or weighting) function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  
  CHARACTER(LEN=*), PARAMETER :: RoutineName                                   &
    = 'SOLVE_BAND_RANDOM_OVERLAP_RESORT_REBIN'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
  
! Set target weights using Gaussian quadrature
  n_esft_red = n_esft_red_in
  IF (n_esft_red > 0 .AND. n_abs > 1) THEN
    IF (control%l_map_sub_bands) THEN
      ! Not compatible with sub-band mapping
      ierr = 1
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    CALL calc_gauss_weight_90(ierr, n_esft_red,                                &
      gpnt_gauleg(1:n_esft_red), w_esft_gauleg(1:n_esft_red))
    IF (gpnt_split < 1.0_RealK .AND. gpnt_split > 0.0_RealK) THEN
      xl=0.5_RealK*(gpnt_split - 0.0_RealK)
      w_esft_target(1:n_esft_red) = xl*w_esft_gauleg(1:n_esft_red)
      xl=0.5_RealK*(1.0_RealK - gpnt_split)
      w_esft_target(n_esft_red+1:2*n_esft_red) =                               &
        xl*w_esft_gauleg(1:n_esft_red)
      n_esft_red = 2*n_esft_red
    ELSE
      xl=0.5_RealK
      w_esft_target(1:n_esft_red) = xl*w_esft_gauleg(1:n_esft_red)
    END IF

! Set target weights using weights of absorber 1 (major absorber)
  ELSE
    k = 1 ! Absorber 1
    i_abs_band=index_abs(k)
    DO iex=1, n_abs_esft(i_abs_band)
      w_esft_target(iex) = w_abs_esft(iex,i_abs_band)
    END DO
    n_esft_red = n_abs_esft(i_abs_band)
  END IF

! Calculate limits in g-space for reduced k-coefficients
  glim_target(1) = 0.0_RealK
  DO iex=1,n_esft_red
    glim_target(iex+1) = glim_target(iex) + w_esft_target(iex)
  END DO

  k = 1 ! Absorber 1
  i_abs_band=index_abs(k)

! Copy ESFT terms and weights for the gas 1 into the arrays for
! the mixture
  k_esft_layer_mix_red(1:n_profile,1:n_layer,1:n_abs_esft(i_abs_band)) =       &
    k_abs_layer(1:n_profile,1:n_layer,1:n_abs_esft(i_abs_band),i_abs_band)
  
! Set current number of reduced ESFT terms
  i_band_esft_mix_red = n_abs_esft(i_abs_band)
! Set current reduced weights
  w_esft_red(1:i_band_esft_mix_red) =                                          &
    w_abs_esft(1:i_band_esft_mix_red,i_abs_band)
  
! Loop over gases and combine ESFT terms and weights
  DO k=2,n_abs
    i_abs_band=index_abs(k)
    
!   Loop over the ESFT terms of the gas to be added
    DO iex=1, n_abs_esft(i_abs_band)
!     Loop over the ESFT terms of the existing mixture and combine with
!     terms for the current gas
      DO iex_mix_red=1,i_band_esft_mix_red    
        iex_mix =                                                              &
          (iex_mix_red - 1)*n_abs_esft(i_abs_band) + iex
        k_esft_layer_mix(1:n_profile,1:n_layer,iex_mix) =                      &
          k_esft_layer_mix_red(1:n_profile,1:n_layer,iex_mix_red) +            &
          k_abs_layer(1:n_profile,1:n_layer,iex,i_abs_band)
        w_esft_mix_unsort(iex_mix) =                                           &
          w_esft_red(iex_mix_red)*w_abs_esft(iex,i_abs_band)
      END DO
    END DO
    i_band_esft_mix = i_band_esft_mix_red*n_abs_esft(i_abs_band)
    
!   Loop over all profiles and layers
    k_esft_layer_mix_red = 0.0_RealK
    DO i_profile=1,n_profile
      DO i_layer=1,n_layer
        
        ! Set weights
        w_esft_mix(1:i_band_esft_mix) =                                        &
          w_esft_mix_unsort(1:i_band_esft_mix)

!       Reorder ESFT terms and weights
        CALL quicksort(i_band_esft_mix,                                        &
          k_esft_layer_mix(i_profile,i_layer,1:i_band_esft_mix),               &
          w_esft_mix(1:i_band_esft_mix))

!       Rebin the ESFT terms
        CALL rebin_esft_terms(i_band_esft_mix, n_esft_red,                     &
          i_profile, i_layer,                                                  &
          w_esft_target, glim_target,                                          &
          k_esft_layer_mix, w_esft_mix,                                        &
          k_esft_layer_mix_red, glim_mix_red,                                  &
          nd_profile, nd_layer, nd_esft_max)
      END DO
    END DO
    
!   After inclusion of gas 2 number of reduced ESFT terms and reduced
!   weights must be reset
    i_band_esft_mix_red = n_esft_red
    w_esft_red(1:n_esft_red) = w_esft_target(1:n_esft_red)
  END DO
  
! Loop over final reduced ESFT terms
  DO iex=1,n_esft_red

!   Set the appropriate source terms for the two-stream equations.

    IF ( (i_angular_integration == ip_two_stream).OR.                          &
         (i_angular_integration == ip_ir_gauss) ) THEN

      IF (isolir == ip_solar) THEN

!       Solar region.
        IF (control%l_spherical_solar) THEN
          DO l=1, n_profile
            d_planck_flux_surface(l) = 0.0e+00_RealK
            flux_inc_down(l)         = 0.0e+00_RealK
            flux_inc_direct(l)       = 0.0e+00_RealK
          END DO
        ELSE
          DO l=1, n_profile
            d_planck_flux_surface(l)=0.0e+00_RealK
            flux_inc_down(l)=solar_irrad(l)/zen_0(l)
            flux_inc_direct(l)=solar_irrad(l)/zen_0(l)
          END DO
        END IF

      ELSE IF (isolir == ip_infra_red) THEN
!       Infra-red region.

        DO l=1, n_profile
          flux_inc_direct(l)=0.0e+00_RealK
          flux_direct_part(l, n_layer)=0.0e+00_RealK
          flux_inc_down(l)=-planck%flux(l, 0)
          d_planck_flux_surface(l)                                             &
            =planck%flux_ground(l)-planck%flux(l, n_layer)
        END DO
        IF (l_clear) THEN
          DO l=1, n_profile
            flux_direct_clear_part(l, n_layer)=0.0e+00_RealK
          END DO
        END IF

      END IF

    ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN

      IF (isolir == ip_solar) THEN
        DO l=1, n_profile
          i_direct_part(l, 0)=solar_irrad(l)
          flux_inc_down(l)=0.0e+00_RealK
        END DO
      ELSE
        DO l=1, n_profile
          flux_inc_down(l)=-planck%flux(l, 0)
          d_planck_flux_surface(l)                                             &
            =planck%flux_ground(l)-planck%flux(l, n_layer)
        END DO
      END IF

    END IF
  
!   Set the absorption for the gas mixture
    k_gas_abs(1:n_profile,1:n_layer) =                                         &
      k_esft_layer_mix_red(1:n_profile,1:n_layer,iex)
  
    CALL monochromatic_radiance(ierr                                           &
      , control, atm, cld, bound                                               &
!                 Atmospheric properties
      , n_profile, n_layer, d_mass                                             &
!                 Angular integration
      , i_angular_integration, i_2stream                                       &
      , l_rescale, n_order_gauss                                               &
      , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc            &
      , accuracy_adaptive, euler_factor                                        &
      , i_sph_algorithm, i_sph_mode                                            &
!                   Precalculated angular arrays
      , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                               &
!                 Treatment of scattering
      , i_scatter_method                                                       &
!                 Options for solver
      , i_solver                                                               &
!                 Gaseous propreties
      , n_k_term_inner_dummy, k_gas_abs                                        &
!                 Options for equivalent extinction
      , .FALSE., dummy_ke                                                      &
!                 Spectral region
      , isolir                                                                 &
!                 Infra-red properties
      , planck                                                                 &
!                 Conditions at TOA
      , zen_0, flux_inc_direct, flux_inc_down                                  &
      , i_direct_part                                                          &
!                 Surface properties
      , d_planck_flux_surface                                                  &
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                               &
      , f_brdf, brdf_sol, brdf_hemi                                            &
!                 Spherical geometry
      , sph                                                                    &
!                 Optical properties
      , ss_prop                                                                &
!                 Cloudy properties
      , l_cloud, i_cloud                                                       &
!                 Cloud geometry
      , n_cloud_top, iex                                                       &
      , n_region, k_clr, i_region_cloud, frac_region                           &
      , w_free, cloud_overlap                                                  &
      , n_column_slv, list_column_slv                                          &
      , i_clm_lyr_chn, i_clm_cld_typ, area_column                              &
!                   Levels for calculating radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                           &
!                   Viewing Geometry
      , n_direction, direction                                                 &
!                 Calculated fluxes
      , flux_direct_part, flux_total_part, l_actinic, actinic_flux_part        &
!                   Calculated radiances
      , radiance_part                                                          &
!                   Calculated rate of photolysis
      , photolysis_part                                                        &
!                 Flags for clear-sky calculations
      , l_clear, i_solver_clear                                                &
!                 Clear-sky fluxes calculated
      , flux_direct_clear_part, flux_total_clear_part                          &
      , actinic_flux_clear_part                                                &
!                 Contribution function
      , contrib_funci_part, contrib_funcf_part                                 &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column                   &
      , nd_flux_profile, nd_radiance_profile, nd_j_profile                     &
      , nd_cloud_type, nd_region, nd_overlap_coeff                             &
      , nd_max_order, nd_sph_coeff                                             &
      , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                     &
      , nd_direction, nd_source_coeff, nd_k_term_inner_dummy                   &
      )

!   Increment the fluxes within the band.
    weight_incr=control%weight_band(i_band)*w_esft_target(iex)
    i_abs_major = index_abs(1)
    weight_sub_band_incr = weight_incr / w_abs_esft(iex, i_abs_major)
    iex_minor = 0
    IF (control%l_blue_flux_surf)                                              &
      weight_blue_incr=spectrum%solar%weight_blue(i_band)*w_esft_target(iex)
    
    CALL augment_radiance(control, spectrum, atm, bound, radout                &
      , i_band, iex, iex_minor                                                 &
      , n_profile, n_layer, n_viewing_level, n_direction                       &
      , l_clear, l_initial, l_initial_band, l_initial_channel                  &
      , weight_incr, weight_blue_incr, weight_sub_band_incr                    &
!                   Actual radiances
      , i_direct                                                               &
!                   Increments to radiances
      , flux_direct_part, flux_total_part, actinic_flux_part                   &
      , i_direct_part, radiance_part, photolysis_part                          &
      , flux_direct_clear_part, flux_total_clear_part                          &
      , actinic_flux_clear_part, k_abs_layer                                   &
      , photol, sph, contrib_funci_part, contrib_funcf_part                    &
!                   Dimensions
      , nd_profile, nd_flux_profile, nd_radiance_profile, nd_j_profile         &
      , nd_layer, nd_viewing_level, nd_direction, dimen%nd_channel             &
      , nd_abs, nd_esft_term                                                   &
      )

!   Add in the increments from surface tiles
    IF (l_tile) THEN
      IF ( (i_angular_integration == ip_two_stream).OR.                        &
           (i_angular_integration == ip_ir_gauss) ) THEN
        IF (control%l_spherical_solar) THEN
          DO l=1, n_profile
            flux_direct_ground_part(l)                                         &
              = sph%allsky%flux_direct(l, n_layer+1)
          END DO
        ELSE
          DO l=1, n_profile
            flux_direct_ground_part(l) = flux_direct_part(l, n_layer)
          END DO          
        END IF
      END IF
      CALL augment_tiled_radiance(control, spectrum, radout                    &
        , i_band, iex, iex_minor                                               &
        , n_point_tile, n_tile, list_tile                                      &
        , l_initial_channel_tile                                               &
        , weight_incr, weight_blue_incr, weight_sub_band_incr                  &
!                   Surface characteristics
        , rho_alb_tile                                                         &
!                   Increments to radiances
        , flux_direct_ground_part                                              &
        , flux_total_part(1, 2*n_layer+2)                                      &
        , planck%flux_tile, planck%flux(:, n_layer)                            &
!                   Dimensions
        , nd_flux_profile, nd_point_tile, nd_tile                              &
        , nd_brdf_basis_fnc, dimen%nd_channel, nd_abs                          &
        )
    END IF
    
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solve_band_random_overlap_resort_rebin
END MODULE solve_band_random_overlap_resort_rebin_mod

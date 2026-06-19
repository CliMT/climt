! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the fluxes within the band with no gases.
!
! Method:
!   Gaseous extinction is set to 0 and a monochromatic
!   calculation is performed.
!
!- ---------------------------------------------------------------------
MODULE solve_band_without_gas_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLVE_BAND_WITHOUT_GAS_MOD'
CONTAINS
SUBROUTINE solve_band_without_gas(ierr                                  &
    , control, dimen, spectrum, atm, cld, bound, radout, i_band         &
!                 Atmospheric column
    , n_profile, n_layer, d_mass                                        &
!                 Angular integration
    , i_angular_integration, i_2stream                                  &
    , n_order_phase, l_rescale, n_order_gauss                           &
    , ms_min, ms_max, i_truncation, ls_local_trunc                      &
    , accuracy_adaptive, euler_factor, i_sph_algorithm                  &
    , i_sph_mode                                                        &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Treatment of scattering
    , i_scatter_method_band                                             &
!                 Options for solver
    , i_solver                                                          &
!                 Spectral region
    , isolir                                                            &
!                 Solar properties
    , zen_0, solar_irrad, sph                                           &
!                 Infra-red properties
    , planck                                                            &
!                 Surface properties
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
!                 Tiling of the surface
    , l_tile, n_point_tile, n_tile, list_tile, rho_alb_tile             &
!                 Optical properties
    , ss_prop, photol                                                   &
!                 Cloudy properties
    , l_cloud, i_cloud                                                  &
!                 Cloud geometry
    , n_cloud_top                                                       &
    , n_region, k_clr, i_region_cloud, frac_region                      &
    , w_free, cloud_overlap                                             &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                 Viewing Geometry
    , n_direction, direction                                            &
!                 Calculated radiances
    , i_direct                                                          &
!                 Flags for initialising diagnostics
    , l_initial, l_initial_band                                         &
    , l_initial_channel, l_initial_channel_tile                         &
!                 Flags for fluxes
    , l_actinic, l_clear, i_solver_clear                                &
!                 Dimensions of arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_cloud_type, nd_region, nd_overlap_coeff                        &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                &
    , nd_direction, nd_source_coeff                                     &
    , nd_point_tile, nd_tile                                            &
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
  USE rad_pcf, ONLY: ip_infra_red, ip_solar, ip_two_stream, ip_ir_gauss, &
                     ip_spherical_harmonic, ip_scatter_hybrid, &
                     ip_scatter_full
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE augment_radiance_mod, ONLY: augment_radiance
  USE augment_tiled_radiance_mod, ONLY: augment_tiled_radiance
  USE monochromatic_radiance_mod, ONLY: monochromatic_radiance

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
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_layer_clr                                                      &
!       Size allocated for totally clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_flux_profile                                                   &
!       Size allocated for profiles in arrays of fluxes
    , nd_radiance_profile                                               &
!       Size allocated for profiles in arrays of radiances
    , nd_j_profile                                                      &
!       Size allocated for profiles in arrays of mean radiances
    , nd_column                                                         &
!       Size allocated for columns per point
    , nd_cloud_type                                                     &
!       Size allocated for types of clouds
    , nd_region                                                         &
!       Size allocated for regions of clouds
    , nd_overlap_coeff                                                  &
!       Size allocated for cloud overlap coefficients
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_sph_coeff                                                      &
!       Size allocated for coefficients of spherical harmonics
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for orders of BRDFs
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_point_tile                                                     &
!       Size allocated for points where the surface is tiled
    , nd_tile
!       Size allocated for surface tiles


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      i_band
!       Band being considered

!                 Atmospheric column
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
    , n_order_phase                                                     &
!       Maximum order of terms in the phase function used in
!       the direct calculation of spherical harmonics
    , n_order_gauss                                                     &
!       Order of gaussian integration
    , ms_min                                                            &
!       Lowest azimuthal order used
    , ms_max                                                            &
!       Highest azimuthal order used
    , i_truncation                                                      &
!       Type of spherical truncation used
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient of (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                   &
!       Orders of truncation at each azimuthal order
    , i_sph_mode                                                        &
!       Mode in which the spherical solver is being used
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
    , uplm_sol(nd_radiance_profile, nd_sph_coeff)                       &
!       Values of spherical harmonics in the solar direction
    , accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series

!                 Treatment of scattering
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method_band
!       Method of treating scattering

!                 Options for solver
  INTEGER, INTENT(IN) ::                                                &
      i_solver
!       Two-stream solver used

!                 Spectral region
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Visible or IR

!                 Solar properties
  REAL (RealK), INTENT(IN) ::                                           &
      zen_0(nd_profile)                                                 &
!       Secants (two-stream) or cosines (spherical harmonics)
!       of the solar zenith angle
    , solar_irrad(nd_profile)
!       Incident solar irradiance in the band

  TYPE(StrPlanck), INTENT(INOUT) :: planck
!   Planckian emission fields

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!   Spherical geometry fields
  
!                 Surface properties
  INTEGER, INTENT(IN) ::                                                &
      ls_brdf_trunc                                                     &
!       Order of truncation of BRDFs
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

! Variables related to tiling of the surface
  LOGICAL, INTENT(IN) ::                                                &
      l_tile
!       Logical to allow invoke options
  INTEGER, INTENT(IN) ::                                                &
      n_point_tile                                                      &
!       Number of points to tile
    , n_tile                                                            &
!       Number of tiles used
    , list_tile(nd_point_tile)
!       List of points with surface tiling
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc, nd_tile)
!       Weights for the basis functions of the BRDFs
!       at the tiled points

!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

  TYPE(StrQy), INTENT(IN) :: photol(spectrum%photol%n_pathway)
!   Photolysis quantum yields interpolated to model grid temperatures

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
    , cloud_overlap(nd_profile, id_ct-1: nd_layer, nd_overlap_coeff)    &
!       Coefficients for transfer for energy at interfaces
    , area_column(nd_profile, nd_column)                                &
!       Areas of columns
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region


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

!                   Calculated radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_radiance_profile, 0: nd_layer)
!       Direct solar irradiance on levels

!                   Flags for initialising diagnostics
  LOGICAL, INTENT(INOUT) ::                                             &
      l_initial                                                         &
!       Initialise rather than increment broadband diagnostics
    , l_initial_band(spectrum%dim%nd_band)                              &
!       Initialise rather than increment band-by-band diagnostics
    , l_initial_channel(dimen%nd_channel)                               &
!       Initialise rather than increment channel diagnostics
    , l_initial_channel_tile(dimen%nd_channel)
!       Initialise rather than increment channel diagnostics on tiles
    
!                   Flags for fluxes
  LOGICAL, INTENT(IN) ::                                                &
      l_actinic
!       Flag for calculation of actinic flux
  LOGICAL, INTENT(IN) ::                                                &
      l_clear
!       Calculate net clear-sky properties
  INTEGER, INTENT(IN) ::                                                &
      i_solver_clear
!       Clear solver used



! Local variables.
  INTEGER :: i, l
!       Loop variables
  INTEGER :: iex, iex_minor(1)
!       Dummy integers for k-term
  INTEGER :: i_scatter_method
!       Method of treating scattering
  INTEGER, PARAMETER ::                                                 &
      n_k_term_inner_dummy = 1                                          &
!       Number of monochromatic calculations in inner loop (dummy here)
    , nd_k_term_inner_dummy = 1
!       Maximum number of k-terms in inner loops (dummy here)
  REAL (RealK) ::                                                       &
      flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)                                         &
!       Incident downward flux
    , d_planck_flux_surface(nd_profile)                                 &
!       Ground source function
    , k_null(nd_profile, nd_layer)                                      &
!       Null vector for call to subroutine
    , dummy_ke(nd_profile, nd_layer)
!       Dummy array (not used)

! Monochromatic incrementing radiances:
  REAL (RealK) ::                                                       &
      flux_direct_band(nd_flux_profile, 0: nd_layer)                    &
!       Increment to direct flux
    , flux_direct_ground_band(nd_flux_profile)                          &
!       Increment to direct flux at the surface
    , flux_total_band(nd_flux_profile, 2*nd_layer+2)                    &
!       Increment to total flux
    , actinic_flux_band(nd_flux_profile, nd_layer)                      &
!       Increment to actinic flux
    , flux_direct_clear_band(nd_flux_profile, 0: nd_layer)              &
!       Increment to clear direct flux
    , flux_total_clear_band(nd_flux_profile, 2*nd_layer+2)              &
!       Increment to clear total flux
    , actinic_flux_clear_band(nd_flux_profile, nd_layer)
!       Increment to clear actinic flux
!                   Increments to Radiances
  REAL (RealK) ::                                                       &
      i_direct_band(nd_radiance_profile, 0: nd_layer)                   &
!       Increments to the solar irradiance
    , radiance_band(nd_radiance_profile, nd_viewing_level               &
        , nd_direction)
!       Increments to the radiance
  REAL (RealK) ::                                                       &
      photolysis_band(nd_j_profile, nd_viewing_level)
!       Increments to the rate of photolysis
  REAL (RealK) ::                                                       &
      weight_incr                                                       &
!       Weight applied to increments
    , weight_blue_incr                                                  &
!       Weight applied to blue increments
    , weight_sub_band_incr
!       Weight applied to sub-band increments

  REAL (RealK) ::                                                       &
      contrib_funci_part(nd_flux_profile, nd_layer)
!       Contribution (or weighting) function
  REAL (RealK) ::                                                       &
      contrib_funcf_part(nd_flux_profile, nd_layer)
!       Contribution (or weighting) function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLVE_BAND_WITHOUT_GAS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the appropriate total upward and downward fluxes
! at the boundaries.

  IF ( (i_angular_integration == ip_two_stream).OR.                     &
       (i_angular_integration == ip_ir_gauss) ) THEN
    IF (isolir == ip_solar) THEN
!     Visible region.
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
!     Infra-red region.
      DO l=1, n_profile
        flux_inc_direct(l)=0.0e+00_RealK
        flux_direct_band(l, n_layer)=0.0e+00_RealK
        flux_inc_down(l)=-planck%flux(l, 0)
        d_planck_flux_surface(l)                                        &
          =planck%flux_ground(l)-planck%flux(l, n_layer)
      END DO
      IF (l_clear) THEN
        DO l=1, n_profile
          flux_direct_clear_band(l, n_layer)=0.0e+00_RealK
        END DO
      END IF
    END IF

  ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN

    IF (isolir == ip_solar) THEN
      DO l=1, n_profile
        i_direct_band(l, 0)=solar_irrad(l)
        flux_inc_down(l)=0.0e+00_RealK
      END DO
    ELSE
      DO l=1, n_profile
        flux_inc_down(l)=-planck%flux(l, 0)
        d_planck_flux_surface(l)                                        &
          =planck%flux_ground(l)-planck%flux(l, n_layer)
      END DO
    END IF

  END IF

  DO i=1, n_layer
    DO l=1, n_profile
      k_null(l, i)=0.0e+00_RealK
    END DO
  END DO

  IF (i_scatter_method_band == ip_scatter_hybrid) THEN
    i_scatter_method = ip_scatter_full
  ELSE
    i_scatter_method = i_scatter_method_band
  END IF


  CALL monochromatic_radiance(ierr                                      &
    , control, atm, cld, bound                                          &
!                 Atmospheric properties
    , n_profile, n_layer, d_mass                                        &
!                 Angular integration
    , i_angular_integration, i_2stream                                  &
    , l_rescale, n_order_gauss                                          &
    , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc       &
    , accuracy_adaptive, euler_factor, i_sph_algorithm                  &
    , i_sph_mode                                                        &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Treatment of scattering
    , i_scatter_method                                                  &
!                 Options for solver
    , i_solver                                                          &
!                 Gaseous propreties
    , n_k_term_inner_dummy, k_null                                      &
!                 Options for equivalent extinction
    , .FALSE., dummy_ke                                                 &
!                 Spectral region
    , isolir                                                            &
!                 Infra-red properties
    , planck                                                            &
!                 Conditions at TOA
    , zen_0, flux_inc_direct, flux_inc_down                             &
    , i_direct_band                                                     &
!                 Surface properties
    , d_planck_flux_surface                                             &
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
!                 Spherical geometry
    , sph                                                               &
!                 Optical properties
    , ss_prop                                                           &
!                 Cloudy properties
    , l_cloud, i_cloud                                                  &
!                 Cloud geometry
    , n_cloud_top, 1                                                    &
    , n_region, k_clr, i_region_cloud, frac_region                      &
    , w_free, cloud_overlap                                             &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                   Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                   Viewing Geometry
    , n_direction, direction                                            &
!                 Calculated Flxues
    , flux_direct_band, flux_total_band, l_actinic, actinic_flux_band   &
!                   Calculated radiances
    , radiance_band                                                     &
!                   Calculated rate of photolysis
    , photolysis_band                                                   &
!                 Flags for clear-sky calculations
    , l_clear, i_solver_clear                                           &
!                 Clear-sky fluxes calculated
    , flux_direct_clear_band, flux_total_clear_band                     &
    , actinic_flux_clear_band                                           &
!                 Contribution function
    , contrib_funci_part, contrib_funcf_part                            &
!                 Dimensions of arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_cloud_type, nd_region, nd_overlap_coeff                        &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                &
    , nd_direction, nd_source_coeff, nd_k_term_inner_dummy              &
    )

! Add the increments to the cumulative fluxes.
  weight_incr = control%weight_band(i_band)
  weight_sub_band_incr = control%weight_band(i_band)
  IF (control%l_blue_flux_surf) &
    weight_blue_incr = spectrum%solar%weight_blue(i_band)
  ! Set k-term to 0 to indicate a single sub-band for the band
  iex = 0
  iex_minor = 0

  CALL augment_radiance(control, spectrum, atm, bound, radout           &
    , i_band, iex, iex_minor                                            &
    , n_profile, n_layer, n_viewing_level, n_direction                  &
    , l_clear, l_initial, l_initial_band, l_initial_channel             &
    , weight_incr, weight_blue_incr, weight_sub_band_incr               &
!                   Actual radiances
    , i_direct                                                          &
!                   Increments to radiances
    , flux_direct_band, flux_total_band, actinic_flux_band              &
    , i_direct_band, radiance_band, photolysis_band                     &
    , flux_direct_clear_band, flux_total_clear_band                     &
    , actinic_flux_clear_band, k_null                                   &
    , photol, sph, contrib_funci_part, contrib_funcf_part               &
!                   Dimensions
    , nd_profile, nd_flux_profile, nd_radiance_profile, nd_j_profile    &
    , nd_layer, nd_viewing_level, nd_direction, dimen%nd_channel        &
    , 1, 1)

! Add in the increments from surface tiles
  IF (l_tile) THEN
    IF ( (i_angular_integration == ip_two_stream).OR.                   &
         (i_angular_integration == ip_ir_gauss) ) THEN
      IF (control%l_spherical_solar) THEN
        DO l=1, n_profile
          flux_direct_ground_band(l)                                    &
            = sph%allsky%flux_direct(l, n_layer+1)
        END DO
      ELSE
        DO l=1, n_profile
          flux_direct_ground_band(l) = flux_direct_band(l, n_layer)
        END DO          
      END IF
    END IF
    CALL augment_tiled_radiance(control, spectrum, radout               &
      , i_band, iex, iex_minor                                          &
      , n_point_tile, n_tile, list_tile                                 &
      , l_initial_channel_tile                                          &
      , weight_incr, weight_blue_incr, weight_sub_band_incr             &
!                   Surface characteristics
      , rho_alb_tile                                                    &
!                   Increments to radiances
      , flux_direct_ground_band                                         &
      , flux_total_band(1, 2*n_layer+2)                                 &
      , planck%flux_tile, planck%flux(:, n_layer)                       &
!                   Dimensions
      , nd_flux_profile, nd_point_tile, nd_tile                         &
      , nd_brdf_basis_fnc, dimen%nd_channel, 1                          &
      )
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solve_band_without_gas
END MODULE solve_band_without_gas_mod

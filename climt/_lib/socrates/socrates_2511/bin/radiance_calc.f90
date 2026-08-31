! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the radiance field.
!
! Method:
!   Properties independent of the spectral bands are set.
!   A loop over bands is then entered, grey optical properties
!   are set and an appropriate subroutine is called to treat
!   the gaseous overlaps. The final radiances are assigned.
!
!------------------------------------------------------------------------------
SUBROUTINE radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)

  USE realtype_rd,  ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, &
    ip_solar, ip_two_stream, ip_sph_reduced_iter, &
    ip_spherical_harmonic, ip_infra_red, i_normal, &
    ip_aerosol_param_moist, ip_aerosol_param_phf_moist, &
    ip_overlap_single, ip_overlap_random, &
    ip_overlap_exact_major, ip_overlap_random_resort_rebin, &
    ip_overlap_k_eqv_scl, ip_overlap_k_eqv, &
    ip_overlap_mix_ses2, ip_cloud_column_max, ip_cloud_mcica, &
    ip_cloud_mix_max, ip_cloud_mix_random, ip_cloud_part_corr, &
    ip_cloud_part_corr_cnv, ip_cloud_triple, &
    ip_overlap_hybrid, ip_overlap_k_eqv_mod, &
    ip_region_clear, ip_region_strat, ip_region_conv, &
    ip_scale_ses2, ip_scale_band, ip_scale_term, &
    ip_scale_lookup, ip_scale_t_lookup, ip_scale_null
  USE def_spectrum, ONLY: StrSpecData
  USE def_dimen,    ONLY: StrDim
  USE def_control,  ONLY: StrCtrl
  USE def_atm,      ONLY: StrAtm
  USE def_cld,      ONLY: StrCld
  USE def_aer,      ONLY: StrAer
  USE def_bound,    ONLY: StrBound
  USE def_out,      ONLY: StrOut, allocate_out
  USE def_ss_prop,  ONLY: str_ss_prop, allocate_ss_prop, deallocate_ss_prop
  USE def_qy,       ONLY: StrQy, allocate_qy
  USE gas_list_pcf, ONLY: ip_h2o, ip_air, molar_weight
  USE rad_ccf,      ONLY: mol_weight_air
  USE errormessagelength_mod, ONLY: errormessagelength
  USE yomhook,      ONLY: lhook, dr_hook
  USE parkind1,     ONLY: jprb, jpim
  USE ereport_mod,  ONLY: ereport

  USE def_planck, ONLY: StrPlanck, allocate_planck, deallocate_planck
  USE def_spherical_geometry, ONLY: StrSphGeo, allocate_sph, deallocate_sph
  USE spherical_path_mod,     ONLY: spherical_path
  USE diff_planck_source_mod, ONLY: diff_planck_source
  USE adjust_ir_radiance_mod, ONLY: adjust_ir_radiance
  USE aggregate_cloud_mod,    ONLY: aggregate_cloud
  USE calc_brdf_mod, ONLY: calc_brdf
  USE calc_cg_coeff_mod, ONLY: calc_cg_coeff
  USE calc_uplm_sol_mod, ONLY: calc_uplm_sol
  USE calc_uplm_zero_mod, ONLY: calc_uplm_zero
  USE check_phf_term_mod, ONLY: check_phf_term
  USE cloud_maxcs_split_mod, ONLY: cloud_maxcs_split
  USE diff_albedo_basis_mod, ONLY: diff_albedo_basis
  USE grey_opt_prop_mod, ONLY: grey_opt_prop
  USE inter_k_mod, ONLY: inter_k
  USE inter_pt_mod, ONLY: inter_pt
  USE inter_pt_lookup_mod, ONLY: inter_pt_lookup
  USE inter_t_lookup_mod, ONLY: inter_t_lookup
  USE overlap_coupled_mod, ONLY: overlap_coupled
  USE rescale_continuum_mod, ONLY: rescale_continuum
  USE rescale_phase_fnc_mod, ONLY: rescale_phase_fnc
  USE scale_absorb_mod, ONLY: scale_absorb
  USE ses_rescale_contm_mod, ONLY: ses_rescale_contm
  USE set_cloud_geometry_mod, ONLY: set_cloud_geometry
  USE set_cloud_pointer_mod, ONLY: set_cloud_pointer
  USE set_rad_layer_mod, ONLY: set_rad_layer
  USE set_truncation_mod, ONLY: set_truncation
  USE sol_scat_cos_mod, ONLY: sol_scat_cos
  USE solve_band_k_eqv_scl_mod, ONLY: solve_band_k_eqv_scl
  USE solve_band_random_overlap_mod, ONLY: solve_band_random_overlap
  USE solve_band_random_overlap_resort_rebin_mod, ONLY: &
      solve_band_random_overlap_resort_rebin
  USE solve_band_ses_mod, ONLY: solve_band_ses
  USE solve_band_without_gas_mod, ONLY: solve_band_without_gas
  USE sum_k_mod, ONLY: sum_k

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),     INTENT(IN)  :: control

! Dimensions:
  TYPE(StrDim),      INTENT(IN)  :: dimen

! Spectral data:
  TYPE(StrSpecData), INTENT(IN)  :: spectrum

! Atmospheric properties:
  TYPE(StrAtm),      INTENT(IN)  :: atm

! Cloud properties:
  TYPE(StrCld),      INTENT(IN)  :: cld

! Aerosol properties:
  TYPE(StrAer),      INTENT(IN)  :: aer

! Boundary conditions:
  TYPE(StrBound),    INTENT(IN)  :: bound

! Output fields:
  TYPE(StrOut),      INTENT(OUT) :: radout


! Local arguments.
! General pointers:
  INTEGER                                                                      &
      i_band                                                                   &
!       Spectral band
    , i_gas                                                                    &
!       Local index of gas
    , n_gas                                                                    &
!       Number of active gases
    , i_gas_band, i_gas_band_1, i_gas_band_2                                   &
!       Single variable for gas in band
    , i_cont_band                                                              &
!       Single variable for generalised continuum in band
    , i_gas_band_sb                                                            &
!       Single variable for index of gas in arrays with self-broadening
    , n_continuum                                                              &
!       Number of continua in band
    , i_continuum                                                              &
!       Continuum number
    , i_pointer_water                                                          &
!       Pointer to water vapour
    , n_cont                                                                   &
!       Number of generalised continua
    , i_abs                                                                    &
!       Local index of absorber, including continuum
    , n_abs                                                                    &
!       Total number of absorbers, including continua
    , i_gas_overlap                                                            &
!       Gas overlap method for band
    , i_path                                                                   &
!       Photolysis pathway
    , i_wl
!       Wavelength subscript

! Dimensions:
  INTEGER                                                                      &
      nd_abs                                                                   &
!       Maximum number of absorbers including continua
    , nd_k_term
!       Maximum number of ESFT terms

! Additional variables for angular integration:
  LOGICAL                                                                      &
      l_solar_phf                                                              &
!       Logical to specify a separate treatment of the singly
!       scattered solar beam
    , l_rescale_solar_phf
!       Logical to apply rescaling to the singly scattered
!       solar phase function
  INTEGER                                                                      &
      n_order_phase
!       Order of Legendre polynomials of the phase function

! Pointers to the contents of layers:
  INTEGER                                                                      &
      n_cloud_top                                                              &
!       Topmost cloudy layer
    , n_region                                                                 &
!       Number of cloudy regions
    , n_cloud_profile(dimen%id_cloud_top: dimen%nd_layer)                      &
!       Number of cloudy profiles
    , i_cloud_profile(dimen%nd_profile, dimen%id_cloud_top: dimen%nd_layer)
!       Profiles containing clouds

! Pointers to types of clouds:
  LOGICAL                                                                      &
      l_cloud_cmp(dimen%nd_cloud_component)
!       Logical switches to `include' components
  INTEGER                                                                      &
      i_phase_cmp(dimen%nd_cloud_component)                                    &
!       Phases of components
    , i_cloud_type(dimen%nd_cloud_component)                                   &
!       Pypes of cloud to which each component contributes
    , type_region(dimen%nd_region)                                             &
!       The types of the regions
    , k_clr                                                                    &
!       Index of clear-sky region
    , i_region_cloud(dimen%nd_cloud_type)
!       Regions in which particular type of cloud fall

! Fractional coverage of different regions:
  REAL (RealK) ::                                                              &
      frac_region(dimen%nd_profile,                                            &
                  dimen%id_cloud_top: dimen%nd_layer,                          &
                  dimen%nd_region)
!       Fraction of total cloud occupied by specific regions

! Pointer to table of humidities:
  INTEGER                                                                      &
      i_humidity_pointer(dimen%nd_profile, dimen%nd_layer)
!       Pointer to look-up table for aerosols

! Controlling variables:
  INTEGER                                                                      &
      i, ii                                                                    &
!       Loop variable
    , j, j_cont                                                                &
!       Loop variable
    , k                                                                        &
!       Loop variable
    , l, ll
!       Loop variable

! Logical switches:
  LOGICAL                                                                      &
      l_gas_band                                                               &
!       Flag to `include' gaseous absorption in a particular band
    , l_cont_band                                                              &
!       Flag to include generalised continuum absorption in a particular band
    , l_abs_band                                                               &
!       Flag to include generalised continuum absorption in a particular band
    , l_moist_aerosol                                                          &
!       Flag for moist aerosol
    , l_clear_band                                                             &
!       Flag to calculate clear-sky fluxes for this band
    , l_water                                                                  &
!       Flag for water vapour present in spectral file
    , l_actinic
!       Flag to calculate actinic flux

  REAL (RealK) ::                                                              &
      solar_irrad_band(dimen%nd_profile)                                       &
!       Solar irradiance in the band
    , solar_irrad_band_ses(dimen%nd_profile, spectrum%dim%nd_k_term)
!       Incident solar flux for each k-term
  REAL (RealK) ::                                                              &
      gas_frac_rescaled(dimen%nd_profile,                                      &
                        dimen%nd_layer,                                        &
                        spectrum%dim%nd_species)                               &
!       Rescaled gas mixing ratios
    , gas_mix_amt(dimen%nd_profile, dimen%nd_layer)                            &
!       Mixed gas mixing ratio
    , amount_continuum(dimen%nd_profile,                                       &
                       dimen%nd_layer,                                         &
                       spectrum%dim%nd_continuum)                              &
!       Amounts of continua
    , k_continuum_mono(spectrum%dim%nd_continuum)
!       Monochromatic continuum components

! Surface BRDF terms
  LOGICAL                                                                      &
      l_diff_alb
!       Flag to calculate diffuse albedos
  REAL (RealK) ::                                                              &
      brdf_sol(dimen%nd_profile, dimen%nd_brdf_basis_fnc, dimen%nd_direction)  &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(dimen%nd_profile, dimen%nd_brdf_basis_fnc, dimen%nd_direction) &
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction
    , diffuse_alb_basis(dimen%nd_brdf_basis_fnc)
!       The diffuse albedo of isotropic radiation for each
!       basis function

! Fields for moist aerosols:
  REAL (RealK) ::                                                              &
      delta_humidity
!       Increment in look-up table for hum.
  INTEGER ::                                                                   &
      nhumidity_common
!       Common number of humidities for moist aerosols

! Fundamental optical properties of layers:
  TYPE(str_ss_prop) :: ss_prop
!   Single scattering properties of the atmosphere

  TYPE(StrPlanck) :: planck
!   Planckian emission fields

  TYPE(StrSphGeo) :: sph
!   Spherical geometry fields

  TYPE(StrQy), ALLOCATABLE :: photol(:)
!   Photolysis quantum yields interpolated to model grid temperatures

  REAL (RealK) ::                                                              &
      k_esft_layer(dimen%nd_profile, dimen%nd_layer, spectrum%dim%nd_k_term,   &
                   spectrum%dim%nd_species)                                    &
!       Exponential ESFT terms at actual pressure layer
    , k_mix_gas_layer(dimen%nd_profile, spectrum%dim%nd_k_term, dimen%nd_layer)&
!       Exponential ESFT terms at actual pressure layer
    , k_contm_layer(dimen%nd_profile, spectrum%dim%nd_k_term, dimen%nd_layer,  &
                    spectrum%dim%nd_continuum)                                 &
!       Continuum absorption coefficients at layer pressure
    , k_esft_cont_layer(dimen%nd_profile, dimen%nd_layer,                      &
                        spectrum%dim%nd_k_term_cont, spectrum%dim%nd_cont)
!       Continuum exponential ESFT terms at actual pressure layer

  REAL (RealK), ALLOCATABLE ::                                                 &
      k_abs_layer(:, :, :, :)                                                  &
!       Exponential ESFT terms for gases and generalised continua at actual
!       pressure layer multiplied by gas mixing ratios
    , k_cont_layer(:, :, :)                                                    &
!       Exponential ESFT terms for generalised continua at actual
!       pressure layer multiplied by gas mixing ratios
    , w_abs_esft(:, :)
!       ESFT weights for gases in band
  INTEGER, ALLOCATABLE ::                                                      &
      index_abs(:)                                                             &
!       Local indexing numbers of gases and continua
    , n_abs_esft(:)                                                            &
!       Number of ESFT terms for gases and continua in a band
    , i_scatter_method_term(:, :)
!       Method of treating scattering for each k-term
  LOGICAL, ALLOCATABLE ::                                                      &
      l_cont_added(:)
!       Flag for having added continuum k-distribution to that of a gas

! Local variables for spherical harmonic integration
  INTEGER                                                                      &
      ls_max_order                                                             &
!       Maximum order of terms required
    , ls_local_trunc(0: dimen%nd_max_order)                                    &
!       Actual truncation for each particular value of m
    , ms_trunc(0: dimen%nd_max_order)                                          &
!       Maximum azimuthal quantum number for each order
    , ia_sph_mm(0: dimen%nd_max_order)
!       Address of spherical coefficient of (m, m) for each m
  REAL (RealK) ::                                                              &
      cg_coeff(dimen%nd_sph_coeff)                                             &
!       Clebsch-gordon coefficients
    , uplm_zero(dimen%nd_sph_coeff)                                            &
!       Upsilon terms
    , uplm_sol(dimen%nd_radiance_profile, dimen%nd_sph_coeff)                  &
!       Upsilon terms for solar radiation
    , cos_sol_view(dimen%nd_radiance_profile, dimen%nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing direction

! Specification of the grid for radiances:
  INTEGER                                                                      &
      i_rad_layer(dimen%nd_viewing_level)
!       Layers in which to intercept radiances
  REAL (RealK) ::                                                              &
      frac_rad_layer(dimen%nd_viewing_level)
!       Fractions below the tops of the layers

  REAL (RealK) ::                                                              &
      i_direct(dimen%nd_radiance_profile, 0: dimen%nd_layer)
!       Direct solar irradiance on levels (not split by
!       diagnostic bands or returned, but retained for
!       future use)

  LOGICAL ::                                                                   &
      l_initial                                                                &
!       Initialise rather than increment broadband diagnostics
    , l_initial_band(spectrum%dim%nd_band)                                     &
!       Initialise rather than increment band-by-band diagnostics
    , l_initial_channel(dimen%nd_channel)                                      &
!       Initialise rather than increment channel diagnostics
    , l_initial_channel_tile(dimen%nd_channel)
!       Initialise rather than increment channel diagnostics on tiles


! Coefficients for the transfer of energy between
! Partially cloudy layers:
  REAL (RealK) ::                                                              &
      cloud_overlap(dimen%nd_profile, dimen%id_cloud_top-1: dimen%nd_layer     &
        , dimen%nd_overlap_coeff)                                              &
!       Coefficients defining overlapping options for clouds:
!       these also depend on the solver selected.
    , w_free(dimen%nd_profile, dimen%id_cloud_top: dimen%nd_layer)
!       Clear-sky fraction

! Cloud geometry
  INTEGER                                                                      &
      n_column_cld(dimen%nd_profile)                                           &
!       Number of columns in each profile (including those of
!       zero width)
    , n_column_slv(dimen%nd_profile)                                           &
!       Number of columns to be solved in each profile
    , list_column_slv(dimen%nd_profile, dimen%nd_column)                       &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(dimen%nd_profile, dimen%nd_column)                         &
!       Layer in the current column to change
    , i_clm_cld_typ(dimen%nd_profile, dimen%nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK) ::                                                              &
      area_column(dimen%nd_profile, dimen%nd_column)
!       Areas of columns

! Secondary arrays for diagnostics
  REAL (RealK) ::                                                              &
      cloud_absorptivity_band(dimen%nd_profile, dimen%nd_layer)                &
!       Absorptivity of cloud in a particular band
    , cloud_extinction_band(dimen%nd_profile, dimen%nd_layer)                  &
!       Extinction of cloud in a particular band
    , ls_cloud_absorptivity_band(dimen%nd_profile, dimen%nd_layer)             &
!       Absorptivity of large-scale cloud in a particular band
    , ls_cloud_extinction_band(dimen%nd_profile, dimen%nd_layer)               &
!       Extinction of large-scale cloud in a particular band
    , cnv_cloud_absorptivity_band(dimen%nd_profile, dimen%nd_layer)            &
!       Absorptivity of convective cloud in a particular band
    , cnv_cloud_extinction_band(dimen%nd_profile, dimen%nd_layer)
!       Extinction of convective cloud in a particular band

  INTEGER ::                                                                   &
      jp(dimen%nd_profile, dimen%nd_layer), jp1                                &
!       Index for pressure interpolation of absorption coefficient
    , jph2oc(dimen%nd_profile, dimen%nd_layer)                                 &
!       Same as JP but for water vapour pressure
    , jt(dimen%nd_profile, dimen%nd_layer), jt1                                &
!       Index for temperature interpolation of absorption coeff
    , jtt(dimen%nd_profile, dimen%nd_layer), jtt1                              &
!       Index of reference temperature at level i+1
!       such that the actual temperature is between JTT and JTT+1
    , jto2c(dimen%nd_profile, dimen%nd_layer)                                  &
!       Index of o2 continuum  reference temperature at level I
!       such that the actual temperature is between JTO2C and JTO2C+1
    , jtswo3(dimen%nd_profile, dimen%nd_layer)                                 &
!       Index of sw o3 reference temp
    , jt_ct(dimen%nd_profile, dimen%nd_layer)                                  &
!       Index for temperature interpolation of generalised continuum
    , jgf(dimen%nd_profile, dimen%nd_layer, spectrum%dim%nd_species_sb)        &
!       Index of reference gas fraction such that the actual gas
!       fraction is between jgf and jgf+1
    , jgfp1(dimen%nd_profile, dimen%nd_layer, spectrum%dim%nd_species_sb)      &
!       Equal to jgf + 1, except if n_gas_frac == 1, in which case it
!       is equal to jgf to avoid out of bounds errors in interpolation.
    , jgf0, jgf1
!       Temporary storage of jgf and jgfp1

  REAL (RealK) ::                                                              &
      fac00(dimen%nd_profile, dimen%nd_layer),                                 &
      fac01(dimen%nd_profile, dimen%nd_layer),                                 &
      fac10(dimen%nd_profile, dimen%nd_layer),                                 &
      fac11(dimen%nd_profile, dimen%nd_layer),                                 &
!       Multiplication factors for P & T interpolation
      fac00c(dimen%nd_profile, dimen%nd_layer),                                &
      fac01c(dimen%nd_profile, dimen%nd_layer),                                &
      fac10c(dimen%nd_profile, dimen%nd_layer),                                &
      fac11c(dimen%nd_profile, dimen%nd_layer),                                &
!       Multiplication factors for H2O cont P & T interpolation
      facc00(dimen%nd_profile, dimen%nd_layer),                                &
      facc01(dimen%nd_profile, dimen%nd_layer),                                &
!       Multiplication factors for O2 continuum T interpolation
      wt_ct(dimen%nd_profile, dimen%nd_layer),                                 &
!       Weight of jt_ct-term in generalised continuum T interpolation
      fgf(dimen%nd_profile, dimen%nd_layer, spectrum%dim%nd_species_sb)
!       Multiplication factors for gas fraction interpolation

! Temperature dependent gas absorption interpolation
  INTEGER :: n_gas_lk, i_gas_lk(spectrum%dim%nd_species)
!   Indexing of gas species for lookup tables
  INTEGER, ALLOCATABLE :: jt_gas(:, :, :)
!   Index for temperature interpolation of gas absorption
  REAL (RealK), ALLOCATABLE :: wt_gas(:, :, :)
!   Weight of jt_gas term in gas absorption temperature interpolation

! Temperature dependent quantum yield interpolation
  INTEGER :: n_t_lookup
!   Number of temperatures in look-up table
  INTEGER, ALLOCATABLE :: jt_qy(:, :)
!   Index for temperature interpolation of quantum yields
  REAL (RealK), ALLOCATABLE :: wt_qy(:, :)
!   Weight of jt_qy term in quantum yield temperature interpolation

  LOGICAL :: l_grey_cont
!   Flag to add continuum in grey_opt_prop

  LOGICAL, ALLOCATABLE :: l_photol_only(:)
!   Only use gas for photolysis, ignoring affect on flux

  INTEGER :: nd_esft_max
!   Maximum number of ESFT terms needed in each band (for
!   arrays when using random overlap with resorting and rebinning)

  INTEGER :: major_abs
!   Variable to temporarily hold index of major absorber

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  INTEGER                       :: ierr = i_normal
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RADIANCE_CALC'


  IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


  CALL allocate_out(radout, control, dimen, spectrum)


! Initial determination of flags and switches:
  IF (control%i_angular_integration == ip_two_stream) THEN

!   Only one term in the phase function is required.
    n_order_phase=1

    l_solar_phf=.FALSE.
    l_rescale_solar_phf=.FALSE.

  ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

!   Set limits on ranges of harmonics and set pointers to arrays.
    CALL set_truncation(ierr                                                   &
      , control%i_truncation, control%ls_global_trunc                          &
      , ls_max_order, ls_local_trunc                                           &
      , control%ms_min, control%ms_max, ms_trunc                               &
      , ia_sph_mm, n_order_phase                                               &
      , dimen%nd_max_order                                                     &
      )

!   Determine whether special treatment of the solar
!   beam is required.
    l_solar_phf=(control%isolir == ip_solar).AND.                              &
                (control%i_sph_algorithm == ip_sph_reduced_iter)
    l_rescale_solar_phf=control%l_rescale.AND.l_solar_phf
!   Calculate the solar scattering angles if treating the
!   solar beam separately.
    IF (l_solar_phf) THEN
      CALL sol_scat_cos(atm%n_profile, atm%n_direction                         &
        , bound%zen_0, atm%direction, cos_sol_view                             &
        , dimen%nd_profile, dimen%nd_direction)
    END IF

!   Calculate Clebsch-Gordan coefficients once and for all.
    CALL calc_cg_coeff(ls_max_order                                            &
      , ia_sph_mm, control%ms_min, ms_trunc                                    &
      , cg_coeff                                                               &
      , dimen%nd_max_order, dimen%nd_sph_coeff)

!   Calculate spherical harmonics at polar angles of pi/2 for
!   use in Marshak's boundary conditions.
    CALL calc_uplm_zero(control%ms_min, control%ms_max, ia_sph_mm              &
      , ls_local_trunc, uplm_zero                                              &
      , dimen%nd_max_order, dimen%nd_sph_coeff)

    IF (control%isolir == ip_solar) THEN
!     Calculate the spherical harmonics of the solar direction.
      CALL calc_uplm_sol(atm%n_profile, control%ms_min, control%ms_max         &
        , ia_sph_mm                                                            &
        , ls_local_trunc, bound%zen_0, uplm_sol                                &
        , dimen%nd_profile, dimen%nd_max_order, dimen%nd_sph_coeff)
    END IF

    IF (control%i_sph_algorithm == ip_sph_reduced_iter) THEN
!     Calcuate some arrays of terms for the BRDF.
      CALL calc_brdf(control%isolir, control%ms_min, control%ms_max            &
        , ia_sph_mm, uplm_sol, uplm_zero                                       &
        , bound%n_brdf_basis_fnc, control%ls_brdf_trunc, bound%f_brdf          &
        , atm%n_profile, atm%n_direction, atm%direction                        &
        , brdf_sol, brdf_hemi                                                  &
        , dimen%nd_profile, dimen%nd_radiance_profile, dimen%nd_direction      &
        , dimen%nd_max_order, dimen%nd_sph_coeff                               &
        , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc)
    END IF

!   For the calculation of equivalent extinction in the IR
!   we need the diffuse albedo for each basis function.
    l_diff_alb=.FALSE.
    DO i_band=1, spectrum%basic%n_band
      l_diff_alb=l_diff_alb.OR.                                                &
        (control%i_gas_overlap_band(i_band) == ip_overlap_k_eqv).OR.           &
        (control%i_gas_overlap_band(i_band) == ip_overlap_k_eqv_scl).OR.       &
        (control%i_gas_overlap_band(i_band) == ip_overlap_k_eqv_mod).OR.       &
        (control%i_gas_overlap_band(i_band) == ip_overlap_hybrid .AND.         &
         spectrum%gas%i_overlap(i_band) == ip_overlap_k_eqv_scl)
    END DO
    IF ( (control%isolir == ip_infra_red).AND.l_diff_alb ) THEN
      CALL diff_albedo_basis(bound%n_brdf_basis_fnc                            &
        , control%ls_brdf_trunc, bound%f_brdf                                  &
        , uplm_zero(ia_sph_mm(0))                                              &
        , diffuse_alb_basis                                                    &
        , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc, dimen%nd_sph_coeff     &
        )
    END IF

!   Determine which layers will be required to give radiances.
    CALL set_rad_layer(ierr                                                    &
      , atm%n_layer, atm%n_viewing_level, atm%viewing_level                    &
      , i_rad_layer, frac_rad_layer                                            &
      , dimen%nd_viewing_level                                                 &
      )


  END IF

! Set the pointer for water vapour to a legal value: this must be done
! for cases where water vapour is not included in the spectral file.
  i_pointer_water=MAX(spectrum%cont%index_water, 1)

! Set flag for the presence of water vapor in the spectral file.
  l_water=ANY(spectrum%gas%type_absorb(1:spectrum%gas%n_absorb) == ip_h2o)


! Initial calculations for aerosols:
  IF (control%l_aerosol) THEN

!   Set the spectrally independent properties of moist aerosols.
    l_moist_aerosol  = .FALSE.
    DO j = 1, spectrum%aerosol%n_aerosol
      SELECT CASE ( spectrum%aerosol%i_aerosol_parm(j) )
      CASE (ip_aerosol_param_moist, ip_aerosol_param_phf_moist)
        l_moist_aerosol  = .TRUE.
        nhumidity_common = spectrum%aerosol%nhumidity(j)
      END SELECT
    END DO

    IF (l_moist_aerosol) THEN
!     Currently all aerosols use the same lookup table of humidities
      delta_humidity = 1.0e+00_RealK                                           &
        / ( REAL(nhumidity_common, RealK) - 1.0e+00_RealK )
      DO i = 1, atm%n_layer
        DO l = 1, atm%n_profile
          i_humidity_pointer(l, i) = 1 +                                       &
            INT( aer%mean_rel_humidity(l, i)*(nhumidity_common-1) )
        END DO
      END DO
    END IF

  END IF


! Initial calculations for clouds:
  IF (control%l_cloud) THEN

!   Set pointers to the types of cloud.
    CALL set_cloud_pointer(ierr                                                &
      , cld%n_condensed, cld%type_condensed, control%i_cloud_representation    &
      , control%l_drop, control%l_ice                                          &
      , i_phase_cmp, i_cloud_type, l_cloud_cmp                                 &
      , dimen%nd_cloud_component                                               &
      )


!   Set the geometry of the clouds.
    CALL set_cloud_geometry(atm%n_profile, atm%n_layer                         &
      , control%l_global_cloud_top, cld%w_cloud                                &
      , n_cloud_top, n_cloud_profile, i_cloud_profile                          &
      , dimen%nd_profile, dimen%nd_layer, dimen%id_cloud_top                   &
      )

    k_clr=1
    IF ( (control%i_cloud == ip_cloud_triple).OR.                              &
         (control%i_cloud == ip_cloud_part_corr_cnv) ) THEN
!     Aggregate clouds into regions for solving.
!     Three regions are used with this option. Additionally,
!     flag the clear-sky region.
      n_region=3
      type_region(1)=ip_region_clear
      type_region(2)=ip_region_strat
      type_region(3)=ip_region_conv
      CALL aggregate_cloud(ierr                                                &
        , atm%n_profile, atm%n_layer, n_cloud_top                              &
        , control%i_cloud, control%i_cloud_representation                      &
        , cld%n_cloud_type, cld%frac_cloud                                     &
        , i_region_cloud, frac_region                                          &
        , dimen%nd_profile, dimen%nd_layer, dimen%nd_cloud_type                &
        , dimen%nd_region, dimen%id_cloud_top                                  &
        )
    ELSE IF ( (control%i_cloud == ip_cloud_mix_max).OR.                        &
              (control%i_cloud == ip_cloud_mix_random).OR.                     &
              (control%i_cloud == ip_cloud_part_corr) ) THEN
!     There will be only one cloudy region.
      n_region=2
      type_region(1)=ip_region_clear
      type_region(2)=ip_region_strat
      DO i=n_cloud_top, atm%n_layer
        DO l=1, atm%n_profile
          frac_region(l, i, 2)=1.0e+00_RealK
        END DO
      END DO
    END IF

!   Calculate energy transfer coefficients in a mixed column,
!   or split the atmosphere into columns with a column model:
    IF ( (control%i_cloud == ip_cloud_mix_max).OR.                             &
         (control%i_cloud == ip_cloud_mix_random).OR.                          &
         (control%i_cloud == ip_cloud_triple).OR.                              &
         (control%i_cloud == ip_cloud_part_corr).OR.                           &
         (control%i_cloud == ip_cloud_part_corr_cnv) ) THEN

      CALL overlap_coupled(atm%n_profile, atm%n_layer, n_cloud_top             &
        , cld%w_cloud, w_free, n_region, type_region, frac_region, atm%p       &
        , control%i_cloud                                                      &
        , cloud_overlap                                                        &
        , dimen%nd_profile, dimen%nd_layer, dimen%nd_overlap_coeff             &
        , dimen%nd_region, dimen%id_cloud_top                                  &
        , cld%dp_corr_strat, cld%dp_corr_conv, radout%tot_cloud_cover          &
        )

    ELSE IF (control%i_cloud == ip_cloud_column_max) THEN

        CALL cloud_maxcs_split(ierr, atm%n_profile, atm%n_layer                &
          , n_cloud_top, cld%w_cloud, cld%frac_cloud                           &
          , cld%n_cloud_type                                                   &
          , n_column_cld, n_column_slv, list_column_slv                        &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                          &
          , dimen%nd_profile, dimen%nd_layer, dimen%id_cloud_top               &
          , dimen%nd_column, dimen%nd_cloud_type                               &
          )

    END IF

  ELSE

    n_cloud_top=atm%n_layer+1

  END IF


! Calculate temperature and pressure interpolation factor for gas k-terms
  IF ( ANY(spectrum%gas%i_scale_fnc(                                           &
    control%first_band : control%last_band, 1) == ip_scale_ses2) ) THEN
    CALL inter_pt(dimen%nd_profile, dimen%nd_layer                             &
      , atm%n_profile, atm%n_layer, atm%gas_mix_ratio(1,1,i_pointer_water)     &
      , atm%p, atm%t, fac00, fac01, fac10, fac11                               &
      , fac00c, fac01c, fac10c, fac11c                                         &
      , facc00, facc01, jp, jph2oc, jt, jtt, jto2c, jtswo3)
  ELSE IF ( ANY(spectrum%gas%i_scale_fnc(                                      &
    control%first_band : control%last_band, :) == ip_scale_lookup) ) THEN
    CALL inter_pt_lookup(dimen%nd_profile, dimen%nd_layer, spectrum%dim%nd_pre &
      , spectrum%dim%nd_tmp, spectrum%dim%nd_gas_frac, spectrum%dim%nd_species &
      , spectrum%dim%nd_species_sb, atm%n_profile, atm%n_layer                 &
      , spectrum%gas%n_gas_frac, spectrum%gas%n_absorb                         &
      , spectrum%gas%type_absorb, spectrum%gas%l_self_broadening               &
      , spectrum%gas%index_sb, l_water, control%l_mixing_ratio                 &
      , i_pointer_water, atm%p, atm%t, atm%gas_mix_ratio                       &
      , spectrum%gas%p_lookup, spectrum%gas%t_lookup, spectrum%gas%gf_lookup   &
      , fac00, fac01, fac10, fac11, jp, jt, jtt, fgf, jgf, jgfp1)
  END IF

! Calculate temperature interpolation factors for gas k-terms
  n_gas_lk = 0
  i_gas_lk(:) = 0
  DO i_gas=1, Spectrum%Gas%n_absorb
    DO i_band=control%first_band, control%last_band
      IF (spectrum%gas%i_scale_fnc(i_band, i_gas) == ip_scale_t_lookup) THEN
        n_gas_lk = n_gas_lk + 1
        i_gas_lk(i_gas) = n_gas_lk
        EXIT
      END IF
    END DO
  END DO
  ALLOCATE(jt_gas(atm%n_profile, atm%n_layer, n_gas_lk))
  ALLOCATE(wt_gas(atm%n_profile, atm%n_layer, n_gas_lk))
  DO i_gas=1, Spectrum%Gas%n_absorb
    ii = i_gas_lk(i_gas)
    IF (ii > 0) THEN
      CALL inter_t_lookup( atm%n_profile, atm%n_layer, &
        Spectrum%Gas%n_t_lookup_gas(i_gas), atm%n_profile, atm%n_layer, &
        atm%t, Spectrum%Gas%t_lookup_gas(:, i_gas), &
        wt_gas(:, :, ii), jt_gas(:, :, ii) )
    END IF
  END DO

! Calculate temperature interpolation factor for continuum k-terms
  IF (ANY(spectrum%contgen%i_band_k_cont(                                      &
            control%first_band : control%last_band, :) > 0)) THEN
    CALL inter_t_lookup(dimen%nd_profile, dimen%nd_layer &
        , spectrum%dim%nd_t_lookup_cont, atm%n_profile, atm%n_layer, atm%t     &
        , spectrum%contgen%t_lookup_cont                                       &
        , wt_ct, jt_ct)
  END IF

! Interpolate quantum yields to model temperatures
  ALLOCATE(photol(spectrum%photol%n_pathway))
  DO i_path=1, spectrum%photol%n_pathway
    n_t_lookup = spectrum%photol%n_t_lookup_photol(i_path)
    IF (n_t_lookup > 1) THEN
      CALL allocate_qy(photol(i_path), atm%n_profile, atm%n_layer, &
                       spectrum%photol%n_wl_lookup_photol(i_path))
      ALLOCATE (jt_qy(atm%n_profile, atm%n_layer))
      ALLOCATE (wt_qy(atm%n_profile, atm%n_layer))
      CALL inter_t_lookup(atm%n_profile, atm%n_layer, &
        n_t_lookup, atm%n_profile, atm%n_layer, atm%t, &
        spectrum%photol%t_lookup_photol(1:n_t_lookup, i_path), &
        wt_qy, jt_qy)
      DO i_wl=1, spectrum%photol%n_wl_lookup_photol(i_path)
        DO i=1, atm%n_layer
          DO l=1, atm%n_profile
            photol(i_path)%qy(l, i, i_wl) = wt_qy(l, i) &
              * spectrum%photol%quantum_yield(jt_qy(l, i), i_wl, i_path) &
              + (1.0_RealK - wt_qy(l, i)) &
              * spectrum%photol%quantum_yield(jt_qy(l, i)+1, i_wl, i_path)
          END DO
        END DO
      END DO
      DEALLOCATE (wt_qy)
      DEALLOCATE (jt_qy)
    END IF
  END DO

! Initialise solar_tail_flux if necessary
  IF (control%l_solar_tail_flux) radout%solar_tail_flux=0.0

! Check that there is enough information in the case of spherical
! harmonics. This check is rather late in the logical order of
! things, but we had to wait for certain other calculations to be
! made.
  IF (control%i_angular_integration == ip_spherical_harmonic) THEN
    CALL check_phf_term(ierr                                                   &
      , control%l_aerosol, spectrum%aerosol%n_aerosol                          &
      , spectrum%aerosol%i_aerosol_parm                                        &
      , spectrum%aerosol%n_aerosol_phf_term                                    &
      , aer%n_phase_term_prsc                                                  &
      , control%l_cloud, cld%n_condensed, cld%i_condensed_param, i_phase_cmp   &
      , cld%condensed_n_phf                                                    &
      , cld%n_phase_term_drop_prsc, cld%n_phase_term_ice_prsc                  &
      , n_order_phase, control%l_henyey_greenstein_pf                          &
      , control%l_rescale, control%n_order_forward                             &
      , l_solar_phf, control%n_order_phase_solar                               &
      , spectrum%dim%nd_aerosol_species, dimen%nd_cloud_component)
  END IF


! Initialization of cloud diagnostics
  IF (control%l_cloud_extinction) THEN
     radout%cloud_extinction        = 0.0
     radout%cloud_weight_extinction = 0.0
  END IF

  IF (control%l_ls_cloud_extinction) THEN
     radout%ls_cloud_extinction        = 0.0
     radout%ls_cloud_weight_extinction = 0.0
  END IF

  IF (control%l_cnv_cloud_extinction) THEN
     radout%cnv_cloud_extinction        = 0.0
     radout%cnv_cloud_weight_extinction = 0.0
  END IF

  IF (control%l_cloud_absorptivity) THEN
     radout%cloud_absorptivity        = 0.0
     radout%cloud_weight_absorptivity = 0.0
  END IF

  IF (control%l_ls_cloud_absorptivity) THEN
     radout%ls_cloud_absorptivity        = 0.0
     radout%ls_cloud_weight_absorptivity = 0.0
  END IF

  IF (control%l_cnv_cloud_absorptivity) THEN
     radout%cnv_cloud_absorptivity        = 0.0
     radout%cnv_cloud_weight_absorptivity = 0.0
  END IF


! Allocate fields and calculate the path for spherical geometry
  IF (control%l_spherical_solar) THEN
    CALL allocate_sph(sph, dimen)
    CALL spherical_path(dimen, atm, bound, sph)
    IF (control%l_spherical_path_diag) THEN
      radout%spherical_path = sph%common%path
    END IF
    radout%flux_direct = 0.0_RealK
    radout%flux_direct_clear = 0.0_RealK
  ELSE
    ! path_div is passed as an argument so must always be allocated
    ALLOCATE(sph%common%path_div(dimen%nd_profile,dimen%nd_layer))
  END IF

! Set flag for calculation of actinic flux
  l_actinic = spectrum%photol%n_pathway > 0 .OR. &
    control%l_actinic_flux .OR. &
    control%l_actinic_flux_band .OR. &
    control%l_actinic_flux_clear .OR. &
    control%l_actinic_flux_clear_band

! On first setting, diagnostics should be initialised
  l_initial=.TRUE.
  l_initial_band(:)=.TRUE.
  l_initial_channel(:)=.TRUE.
  l_initial_channel_tile(:)=.TRUE.

! Solve the equation of transfer in each band and increment the fluxes.
  DO i_band=control%first_band, control%last_band

!   Determine whether clear-sky fluxes are required for this band
    l_clear_band = control%l_clear .OR. control%l_clear_band(i_band)

!   Determine the gas overlap method for the band
    IF (control%i_gas_overlap_band(i_band) == ip_overlap_hybrid) THEN
      i_gas_overlap = spectrum%gas%i_overlap(i_band)
    ELSE
      i_gas_overlap = control%i_gas_overlap_band(i_band)
    END IF

!   Rescale amounts of continua.
    l_grey_cont = .FALSE.
    IF (control%l_continuum) THEN
      n_continuum=spectrum%cont%n_band_continuum(i_band)
      DO i=1, n_continuum
        i_continuum=spectrum%cont%index_continuum(i_band, i)
        IF (spectrum%cont%i_scale_fnc_cont(i_band, i) == ip_scale_ses2) THEN
          k_continuum_mono(i)=0.0
          CALL ses_rescale_contm(dimen%nd_profile, dimen%nd_layer              &
            , i_continuum, atm%n_profile, atm%n_layer                          &
            , atm%p, atm%t, atm%gas_mix_ratio(1,1,i_pointer_water)             &
            , amount_continuum(1, 1, i)                                        &
            )
        ELSE
          l_grey_cont = .TRUE.
          k_continuum_mono(i) = spectrum%cont%k_cont(i_band, i)
          CALL rescale_continuum(control, atm%n_profile, atm%n_layer           &
            , i_continuum, atm%p, atm%t                                        &
            , atm%density, atm%gas_mix_ratio(1, 1, i_pointer_water)            &
            , amount_continuum(1, 1, i)                                        &
            , spectrum%cont%i_scale_fnc_cont(i_band, i)                        &
            , spectrum%cont%p_ref_cont(i, i_band)                              &
            , spectrum%cont%t_ref_cont(i, i_band)                              &
            , spectrum%cont%scale_cont(1, i_band, i)                           &
            , dimen%nd_profile, dimen%nd_layer                                 &
            , spectrum%dim%nd_scale_variable)
        END IF
      END DO
    END IF

!   Allocate the single scattering propeties.
    CALL allocate_ss_prop(ss_prop, dimen)

!   Calculate the grey extinction within the band.

    CALL grey_opt_prop(ierr, control, radout, i_band                           &
      , atm%n_profile, atm%n_layer, atm%p, atm%t, atm%density                  &
      , n_order_phase, l_solar_phf, atm%n_direction, cos_sol_view              &
      , spectrum%rayleigh%i_rayleigh_scheme, spectrum%rayleigh%n_gas_rayleigh  &
      , spectrum%rayleigh%index_rayleigh                                       &
      , spectrum%rayleigh%rayleigh_coeff(i_band)                               &
      , spectrum%rayleigh%rayleigh_coeff_gas(1, i_band)                        &
      , atm%gas_mix_ratio                                                      &
      , l_grey_cont, n_continuum                                               &
      , k_continuum_mono, amount_continuum                                     &
      , spectrum%aerosol%n_aerosol_mr, aer%mix_ratio                           &
      , aer%mr_source, aer%mr_type_index                                       &
      , spectrum%aerosol%i_aerosol_parm                                        &
      , i_humidity_pointer, spectrum%aerosol%humidities, delta_humidity        &
      , aer%mean_rel_humidity                                                  &
      , spectrum%aerosol%abs(1, 1, i_band)                                     &
      , spectrum%aerosol%scat(1, 1, i_band)                                    &
      , spectrum%aerosol%phf_fnc(1, 1, 1, i_band)                              &
      , aer%n_mode, aer%mode_mix_ratio                                         &
      , aer%mode_absorption(1, 1, 1, i_band)                                   &
      , aer%mode_scattering(1, 1, 1, i_band)                                   &
      , aer%mode_asymmetry(1, 1, 1, i_band)                                    &
      , aer%n_opt_level_prsc, aer%pressure_prsc                                &
      , aer%absorption_prsc(1, 1, 1, i_band)                                   &
      , aer%scattering_prsc(1, 1, 1, i_band)                                   &
      , aer%phase_fnc_prsc(1, 1, 1, 1, i_band)                                 &
      , n_cloud_profile, i_cloud_profile                                       &
      , n_cloud_top, cld%n_condensed, l_cloud_cmp, i_phase_cmp                 &
      , cld%i_condensed_param                                                  &
      , cld%condensed_param_list(1, 1, i_band)                                 &
      , cld%condensed_mix_ratio, cld%condensed_dim_char                        &
      , cld%condensed_rel_var_dens                                             &
      , cld%n_cloud_type, i_cloud_type                                         &
      , cld%n_opt_level_drop_prsc, cld%drop_pressure_prsc                      &
      , cld%drop_absorption_prsc(1, 1, i_band)                                 &
      , cld%drop_scattering_prsc(1, 1, i_band)                                 &
      , cld%drop_phase_fnc_prsc(1, 1, 1, i_band)                               &
      , cld%n_opt_level_ice_prsc, cld%ice_pressure_prsc                        &
      , cld%ice_absorption_prsc(1, 1, i_band)                                  &
      , cld%ice_scattering_prsc(1, 1, i_band)                                  &
      , cld%ice_phase_fnc_prsc(1, 1, 1, i_band)                                &
      , ss_prop                                                                &
      , cld%frac_cloud                                                         &
      , cloud_extinction_band, cloud_absorptivity_band                         &
      , ls_cloud_extinction_band, ls_cloud_absorptivity_band                   &
      , cnv_cloud_extinction_band, cnv_cloud_absorptivity_band                 &
      , dimen%nd_profile, dimen%nd_radiance_profile, dimen%nd_layer            &
      , dimen%nd_layer_clr, dimen%id_cloud_top                                 &
      , spectrum%dim%nd_continuum, spectrum%dim%nd_species                     &
      , spectrum%dim%nd_aerosol_species                                        &
      , spectrum%dim%nd_aerosol_mr, spectrum%dim%nd_humidity                   &
      , spectrum%dim%nd_cloud_parameter, dimen%nd_cloud_component              &
      , dimen%nd_cloud_type, spectrum%dim%nd_phase_term                        &
      , dimen%nd_phf_term_aerosol_prsc, dimen%nd_phf_term_cloud_prsc           &
      , dimen%nd_max_order, dimen%nd_direction, dimen%nd_aerosol_mode          &
      , dimen%nd_profile_aerosol_prsc, dimen%nd_profile_cloud_prsc             &
      , dimen%nd_opt_level_aerosol_prsc, dimen%nd_opt_level_cloud_prsc         &
      )


    IF ( (control%i_angular_integration == ip_two_stream).OR.                  &
         (control%i_angular_integration == ip_spherical_harmonic) ) THEN

!     Rescale the phase function and calculate the scattering
!     fractions. (These are grey and may be calculated outside
!     a loop over gases).

      IF (control%l_rescale) THEN

!       Rescale clear-sky phase function:

!       The section above clouds.
        CALL rescale_phase_fnc(atm%n_profile, 1, n_cloud_top-1                 &
          , atm%n_direction, cos_sol_view                                      &
          , n_order_phase                                                      &
          , ss_prop%phase_fnc_clr, ss_prop%forward_scatter_clr                 &
          , ss_prop%forward_solar_clr                                          &
          , l_rescale_solar_phf, control%n_order_phase_solar                   &
          , ss_prop%phase_fnc_solar_clr                                        &
          , dimen%nd_profile, dimen%nd_radiance_profile, dimen%nd_layer_clr, 1 &
          , dimen%nd_direction, dimen%nd_max_order                             &
          )
!       The section including clouds.
        CALL rescale_phase_fnc(atm%n_profile, n_cloud_top                      &
          , atm%n_layer, atm%n_direction, cos_sol_view                         &
          , n_order_phase                                                      &
          , ss_prop%phase_fnc(:, :, :, 0)                                      &
          , ss_prop%forward_scatter(:, :, 0)                                   &
          , ss_prop%forward_solar(:, :, 0)                                     &
          , l_rescale_solar_phf, control%n_order_phase_solar                   &
          , ss_prop%phase_fnc_solar(:, :, :, 0)                                &
          , dimen%nd_profile, dimen%nd_radiance_profile, dimen%nd_layer        &
          , dimen%id_cloud_top, dimen%nd_direction, dimen%nd_max_order         &
          )


        IF (control%l_cloud .AND. control%i_cloud /= ip_cloud_mcica) THEN

!         Rescale cloudy phase functions:
          DO k=1, cld%n_cloud_type
            CALL rescale_phase_fnc(atm%n_profile, n_cloud_top                  &
              , atm%n_layer, atm%n_direction, cos_sol_view                     &
              , n_order_phase                                                  &
              , ss_prop%phase_fnc(:, :, :, k)                                  &
              , ss_prop%forward_scatter(:, :, k)                               &
              , ss_prop%forward_solar(:, :, k)                                 &
              , l_rescale_solar_phf, control%n_order_phase_solar               &
              , ss_prop%phase_fnc_solar(:, :, :, k)                            &
              , dimen%nd_profile, dimen%nd_radiance_profile, dimen%nd_layer    &
              , dimen%id_cloud_top, dimen%nd_direction, dimen%nd_max_order     &
              )
          END DO

        END IF

      END IF

    END IF


!   Preliminary calculations for source terms:
    IF (i_gas_overlap == ip_overlap_mix_ses2) THEN

!     Interpolate absorption coefficients onto model grid
      CALL inter_k(atm%n_profile, atm%n_layer, spectrum%gas%n_band_absorb      &
        , spectrum%gas%mix_gas_band, spectrum%gas%n_mix_gas(i_band)            &
        , spectrum%gas%index_mix_gas                                           &
        , spectrum%gas%i_band_k_ses(i_band), i_band                            &
        , n_continuum, spectrum%cont%k_cont_ses, control%l_continuum           &
        , spectrum%cont%index_continuum                                        &
        , spectrum%cont%k_h2oc(1,1,1,i_band)                                   &
        , fac00, fac01, fac10, fac11                                           &
        , fac00c, fac01c, fac10c, fac11c, facc00, facc01                       &
        , jp, jph2oc, jt, jtt, jto2c, jtswo3                                   &
        , spectrum%gas%k_lookup(1,1,1,1,i_band), spectrum%gas%k_mix_gas        &
        , spectrum%gas%f_mix(i_band), atm%gas_mix_ratio, gas_mix_amt           &
        , k_esft_layer, k_mix_gas_layer, k_contm_layer                         &
!       Dimensions
        , dimen%nd_profile, dimen%nd_layer                                     &
        , spectrum%dim%nd_band, spectrum%dim%nd_species                        &
        , spectrum%dim%nd_continuum                                            &
        , spectrum%dim%nd_k_term, spectrum%dim%nd_mix                          &
        , spectrum%dim%nd_tmp, spectrum%dim%nd_pre                             &
        , spectrum%dim%nd_band_mix_gas)

      IF ((control%isolir == ip_solar) .OR. control%l_solar_tail_flux) THEN
!       Convert normalized band fluxes to actual energy fluxes.
        DO k=1,spectrum%gas%i_band_k_ses(i_band)
          DO l=1, atm%n_profile
            solar_irrad_band_ses(l,k)=bound%solar_irrad(l)                     &
              *spectrum%solar%solar_flux_band_ses(k, i_band)
          END DO
        END DO
      END IF

      IF (control%l_solar_tail_flux) THEN
!       Calculate solar tail flux and add it to the solar region
!       for diagnostic output
        DO k=1,spectrum%gas%i_band_k_ses(i_band)
          DO l=1, atm%n_profile
            radout%solar_tail_flux(l)=radout%solar_tail_flux(l)                &
              +spectrum%gas%w_ses(k,i_band)*solar_irrad_band_ses(l,k)
          END DO
        END DO
      END IF

    ELSE

      DO j=1, spectrum%gas%n_band_absorb(i_band)
        i_gas_band=spectrum%gas%index_absorb(j, i_band)
        IF (spectrum%gas%i_scale_fnc(i_band, i_gas_band)                       &
          == ip_scale_lookup) THEN
          IF (spectrum%gas%l_self_broadening(i_gas_band)) THEN
            i_gas_band_sb=spectrum%gas%index_sb(i_gas_band)
            DO k=1, spectrum%gas%i_band_k(i_band, i_gas_band)
              DO i=1, atm%n_layer
                DO l=1, atm%n_profile
                  jp1=jp(l,i)+1
                  jt1=jt(l,i)+1
                  jtt1=jtt(l,i)+1
                  jgf0=jgf(l,i,i_gas_band_sb)
                  jgf1=jgfp1(l,i,i_gas_band_sb)
                  k_esft_layer(l,i,k,i_gas_band) = MAX(0.0_RealK,              &
                      fgf(l,i,i_gas_band_sb)                                   &
                    *(fac00(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jt(l,i),  jp(l,i), jgf0, k, i_gas_band_sb, i_band )      &
                    + fac10(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jtt(l,i), jp1,     jgf0, k, i_gas_band_sb, i_band )      &
                    + fac01(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jt1,      jp(l,i), jgf0, k, i_gas_band_sb, i_band )      &
                    + fac11(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jtt1,     jp1,     jgf0, k, i_gas_band_sb, i_band ) )    &
                    +(1.0_RealK - fgf(l,i,i_gas_band_sb))                      &
                    *(fac00(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jt(l,i),  jp(l,i), jgf1, k, i_gas_band_sb, i_band )      &
                    + fac10(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jtt(l,i), jp1,     jgf1, k, i_gas_band_sb, i_band )      &
                    + fac01(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jt1,      jp(l,i), jgf1, k, i_gas_band_sb, i_band )      &
                    + fac11(l,i)*spectrum%gas%k_lookup_sb(                     &
                      jtt1,     jp1,     jgf1, k, i_gas_band_sb, i_band ) ) )
                END DO
              END DO
            END DO
          ELSE
            DO k=1, spectrum%gas%i_band_k(i_band, i_gas_band)
              DO i=1, atm%n_layer
                DO l=1, atm%n_profile
                  jp1=jp(l,i)+1
                  jt1=jt(l,i)+1
                  jtt1=jtt(l,i)+1
                  k_esft_layer(l,i,k,i_gas_band) = MAX(0.0_RealK,              &
                      fac00(l,i)*spectrum%gas%k_lookup(                        &
                      jt(l,i),  jp(l,i), k, i_gas_band, i_band )               &
                    + fac10(l,i)*spectrum%gas%k_lookup(                        &
                      jtt(l,i), jp1,     k, i_gas_band, i_band )               &
                    + fac01(l,i)*spectrum%gas%k_lookup(                        &
                      jt1,      jp(l,i), k, i_gas_band, i_band )               &
                    + fac11(l,i)*spectrum%gas%k_lookup(                        &
                      jtt1,     jp1,     k, i_gas_band, i_band ) )
                END DO
              END DO
            END DO
          END IF
        ELSE IF (spectrum%gas%i_scale_fnc(i_band, i_gas_band) &
          == ip_scale_t_lookup) THEN
          ii = i_gas_lk(i_gas_band)
          DO k=1, spectrum%gas%i_band_k(i_band, i_gas_band)
            DO i=1, atm%n_layer
              DO l=1, atm%n_profile
                k_esft_layer(l,i,k,i_gas_band) = MAX(0.0_RealK, &
                  wt_gas(l,i,ii)*Spectrum%Gas%k_t_lookup_gas( &
                    jt_gas(l,i,ii), k, i_gas_band, i_band ) &
                  + (1.0_RealK - wt_gas(l,i,ii))*Spectrum%Gas%k_t_lookup_gas( &
                    jt_gas(l,i,ii)+1, k, i_gas_band, i_band) )
              END DO
            END DO
          END DO
        END IF
      END DO

      IF (control%isolir == ip_solar) THEN
!       Convert normalized band fluxes to actual energy fluxes.
        IF (control%l_spherical_solar) THEN
          DO l=1, atm%n_profile
            ! For the surface and TOA it's OK to use cos_zen, which may equal
            ! zero, as in that case zero flux would be appropriate.
            sph%common%flux_inc_direct(l,0) = bound%solar_irrad(l)             &
              * bound%lit(l,0) * spectrum%solar%solar_flux_band(i_band)        &
              * bound%cos_zen(l,0)
          END DO
          DO i=1, atm%n_layer
            DO l=1, atm%n_profile
              ! For the flux arriving at the layers we solve directly along
              ! the beam direction but scale the area normal to the beam to
              ! maintain a constant volume for the column element.
              sph%common%flux_inc_direct(l,i) = bound%solar_irrad(l)           &
                * bound%lit(l,i) * spectrum%solar%solar_flux_band(i_band)      &
                / sph%common%path_div(l,i)
            END DO
          END DO
          i=atm%n_layer+1
          DO l=1, atm%n_profile
            sph%common%flux_inc_direct(l,i) = bound%solar_irrad(l)             &
              * bound%lit(l,i) * spectrum%solar%solar_flux_band(i_band)        &
              * bound%cos_zen(l,i)
            solar_irrad_band(l)=0.0_RealK
          END DO
        ELSE
          DO l=1, atm%n_profile
            solar_irrad_band(l)=bound%solar_irrad(l)                           &
              *spectrum%solar%solar_flux_band(i_band)
          END DO
        END IF
      END IF

    END IF


!   Perform interpolation of generalised continuum ESFT/k-terms look-up table
!   to model grid.
    DO j=1, spectrum%contgen%n_band_cont(i_band)
      i_cont_band=spectrum%contgen%index_cont(j, i_band)
      DO k=1, spectrum%contgen%i_band_k_cont(i_band, i_cont_band)
        DO i=1, atm%n_layer
          DO l=1, atm%n_profile
            k_esft_cont_layer(l,i,k,i_cont_band) = MAX(0.0_RealK,              &
                wt_ct(l,i)*spectrum%contgen%k_lookup_cont(                     &
                jt_ct(l,i), k, i_cont_band, i_band )                           &
              + (1.0_RealK - wt_ct(l,i))*spectrum%contgen%k_lookup_cont(       &
                jt_ct(l,i)+1, k, i_cont_band, i_band) )
          END DO
        END DO
      END DO
    END DO

!   The foreign broadened water vapour continuum is treated as a special case:
!   The mean molar weight of dry air was used to derive the continuum
!   ESFT/k-terms, however, this mean molecular weight may have been changed
!   here. The terms are therefore scaled to the correct mean molar weight of
!   dry air here. This ensures a broadening of water vapour by all other gases
!   even if the atmosphere is not composed of pure air.
    DO j=1, spectrum%contgen%n_band_cont(i_band)
      i_cont_band=spectrum%contgen%index_cont(j, i_band)
      i_gas_band=spectrum%contgen%index_cont_gas_1(i_cont_band)
      IF (spectrum%gas%type_absorb(i_gas_band) == ip_air) THEN
        DO k=1, spectrum%contgen%i_band_k_cont(i_band, i_cont_band)
          DO i=1, atm%n_layer
            DO l=1, atm%n_profile
              k_esft_cont_layer(l,i,k,i_cont_band)                             &
                =k_esft_cont_layer(l,i,k,i_cont_band)                          &
                *molar_weight(ip_air)*1.0E-03_RealK/mol_weight_air
            END DO
          END DO
        END DO
      END IF
      i_gas_band=spectrum%contgen%index_cont_gas_2(i_cont_band)
      IF (spectrum%gas%type_absorb(i_gas_band) == ip_air) THEN
        DO k=1, spectrum%contgen%i_band_k_cont(i_band, i_cont_band)
          DO i=1, atm%n_layer
            DO l=1, atm%n_profile
              k_esft_cont_layer(l,i,k,i_cont_band)                             &
                =k_esft_cont_layer(l,i,k,i_cont_band)                          &
                *molar_weight(ip_air)*1.0E-03_RealK/mol_weight_air
            END DO
          END DO
        END DO
      END IF
    END DO

!   Determine whether gaseous absorption is included in this band.
    IF ((control%l_gas).AND.(spectrum%gas%n_band_absorb(i_band) > 0)) THEN
      n_gas=spectrum%gas%n_band_absorb(i_band)
      l_gas_band=.TRUE.
    ELSE
      n_gas=0
      l_gas_band=.FALSE.
    END IF

!   Determine whether generalised continuum absorption is included in
!   this band.
    IF (control%l_cont_gen .AND. spectrum%contgen%n_band_cont(i_band) > 0) THEN
      n_cont=spectrum%contgen%n_band_cont(i_band)
      IF (n_cont > 0) THEN
!       Set the flag for generalised continuum absorption in the band.
        l_cont_band=.TRUE.

!       Generalised continuum absorption requires the complete scaling of
!       ESFT terms to be performed here.
        IF (i_gas_overlap == ip_overlap_mix_ses2) THEN
          cmessage =                                                           &
            '*** Error: The selected gaseous overlap method is invalid.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

      ELSE
        l_cont_band=.FALSE.
      END IF
    ELSE
      n_cont=0
      l_cont_band=.FALSE.
    END IF

!   Absorbers includes both gases and generalised continua.
    l_abs_band=l_gas_band.OR.l_cont_band

!   Determine the number of absorbers and the maximum number of ESFT
!   terms needed.
    n_abs=n_gas+n_cont
    IF (l_cont_band) THEN
      DO j=1, n_gas
        i_gas_band=spectrum%gas%index_absorb(j, i_band)
        DO j_cont=1, n_cont
          i_cont_band=spectrum%contgen%index_cont(j_cont, i_band)
          IF (spectrum%contgen%i_cont_overlap_band(i_band, i_cont_band) ==     &
              i_gas_band) THEN
            n_abs=n_abs-1
          END IF
        END DO
      END DO
    END IF
    nd_abs=n_abs
    nd_k_term=MAX(spectrum%dim%nd_k_term, spectrum%dim%nd_k_term_cont)

!   Allocate absorber arrays
    ALLOCATE(k_abs_layer(dimen%nd_profile, dimen%nd_layer, nd_k_term, nd_abs))
    ALLOCATE(k_cont_layer(dimen%nd_profile, dimen%nd_layer,                    &
                          spectrum%dim%nd_k_term_cont))
    ALLOCATE(w_abs_esft(nd_k_term, nd_abs))
    ALLOCATE(index_abs(nd_abs))
    ALLOCATE(n_abs_esft(nd_abs))
    ALLOCATE(i_scatter_method_term(nd_k_term, nd_abs))
    ALLOCATE(l_cont_added(spectrum%dim%nd_cont))
    ALLOCATE(l_photol_only(nd_abs))

!   Get gaseous absorption data for this band.
    l_cont_added=.FALSE.
    l_photol_only=.FALSE.
    IF (l_gas_band) THEN
      DO j=1, n_gas

        i_gas_band=spectrum%gas%index_absorb(j, i_band)

        index_abs(j)=j
        n_abs_esft(j)=spectrum%gas%i_band_k(i_band, i_gas_band)
        DO k=1, n_abs_esft(j)
          i_scatter_method_term(k, j)                                          &
            = spectrum%gas%i_scat(k, i_band, i_gas_band)
          w_abs_esft(k, j) = spectrum%gas%w(k, i_band, i_gas_band)
        END DO

        l_photol_only(j) &
          = control%l_photol_only(spectrum%gas%type_absorb(i_gas_band))

!       Perform the rescaling of the amount of gas for this band.
        IF (.NOT.l_photol_only(j) .AND. &
            i_gas_overlap == ip_overlap_single .AND. &
            spectrum%gas%type_absorb(i_gas_band) /= control%i_gas) THEN
          ! If only considering one gas zero the absorption for other species
          DO k=1, n_abs_esft(j)
            DO i=1, atm%n_layer
              DO l=1, atm%n_profile
                k_abs_layer(l, i, k, j) = 0.0_RealK
              END DO
            END DO
          END DO

        ELSE IF (spectrum%gas%i_scale_k(i_band, i_gas_band)                    &
            == ip_scale_band) THEN
          CALL scale_absorb(ierr, atm%n_profile, atm%n_layer                   &
            , atm%p, atm%t                                                     &
            , gas_frac_rescaled(1, 1, i_gas_band)                              &
            , spectrum%gas%i_scale_fnc(i_band, i_gas_band)                     &
            , spectrum%gas%p_ref(i_gas_band, i_band)                           &
            , spectrum%gas%t_ref(i_gas_band, i_band)                           &
            , spectrum%gas%scale(1, 1, i_band, i_gas_band)                     &
            , 1, i_band, spectrum%gas%l_doppler(i_gas_band)                    &
            , spectrum%gas%doppler_cor(i_gas_band)                             &
            , dimen%nd_profile, dimen%nd_layer                                 &
            , spectrum%dim%nd_scale_variable)
          IF (l_photol_only(j)) THEN
            DO i=1, atm%n_layer
              DO l=1, atm%n_profile
                gas_frac_rescaled(l, i, i_gas_band) = MAX(0.0_RealK, &
                  gas_frac_rescaled(l, i, i_gas_band))
              END DO
            END DO
          ELSE
            DO i=1, atm%n_layer
              DO l=1, atm%n_profile
                gas_frac_rescaled(l, i, i_gas_band) = MAX(0.0_RealK, &
                  gas_frac_rescaled(l, i, i_gas_band) &
                  * atm%gas_mix_ratio(l, i, i_gas_band))
              END DO
            END DO
          END IF
          DO k=1, n_abs_esft(j)
            DO i=1, atm%n_layer
              DO l=1, atm%n_profile
                k_abs_layer(l, i, k, j)                                        &
                  = spectrum%gas%k(k, i_band, i_gas_band)                      &
                  * gas_frac_rescaled(l, i, i_gas_band)
              END DO
            END DO
          END DO

        ELSE IF (spectrum%gas%i_scale_k(i_band, i_gas_band)                    &
            == ip_scale_null) THEN
!         Copy across the unscaled array.
          IF (l_photol_only(j)) THEN
            DO k=1, n_abs_esft(j)
              DO i=1, atm%n_layer
                DO l=1, atm%n_profile
                  k_abs_layer(l, i, k, j)                                      &
                    = spectrum%gas%k(k, i_band, i_gas_band)
                END DO
              END DO
            END DO
          ELSE
            DO k=1, n_abs_esft(j)
              DO i=1, atm%n_layer
                DO l=1, atm%n_profile
                  k_abs_layer(l, i, k, j)                                      &
                    = spectrum%gas%k(k, i_band, i_gas_band)                    &
                    * atm%gas_mix_ratio(l, i, i_gas_band)
                END DO
              END DO
            END DO
          END IF

        ELSE IF (spectrum%gas%i_scale_k(i_band, i_gas_band)                    &
            == ip_scale_term) THEN
          IF (spectrum%gas%i_scale_fnc(i_band, i_gas_band) &
                == ip_scale_lookup .OR. &
              spectrum%gas%i_scale_fnc(i_band, i_gas_band) &
                == ip_scale_t_lookup) THEN
            IF (l_photol_only(j)) THEN
              DO k=1, n_abs_esft(j)
                DO i=1, atm%n_layer
                  DO l=1, atm%n_profile
                    k_abs_layer(l, i, k, j) = k_esft_layer(l, i, k, i_gas_band)
                  END DO
                END DO
              END DO
            ELSE
              DO k=1, n_abs_esft(j)
                DO i=1, atm%n_layer
                  DO l=1, atm%n_profile
                    k_abs_layer(l, i, k, j) = MAX(0.0_RealK, &
                      k_esft_layer(l, i, k, i_gas_band) &
                      * atm%gas_mix_ratio(l, i, i_gas_band))
                  END DO
                END DO
              END DO
            END IF
          ELSE
            DO k=1, n_abs_esft(j)
              CALL scale_absorb(ierr, atm%n_profile, atm%n_layer              &
                , atm%p, atm%t                                                &
                , gas_frac_rescaled(1, 1, i_gas_band)                         &
                , spectrum%gas%i_scale_fnc(i_band, i_gas_band)                &
                , spectrum%gas%p_ref(i_gas_band, i_band)                      &
                , spectrum%gas%t_ref(i_gas_band, i_band)                      &
                , spectrum%gas%scale(1, k, i_band, i_gas_band)                &
                , k, i_band, spectrum%gas%l_doppler(i_gas_band)               &
                , spectrum%gas%doppler_cor(i_gas_band)                        &
                , dimen%nd_profile, dimen%nd_layer                            &
                , spectrum%dim%nd_scale_variable)
              IF (l_photol_only(j)) THEN
                DO i=1, atm%n_layer
                  DO l=1, atm%n_profile
                    k_abs_layer(l, i, k, j) &
                      = spectrum%gas%k(k, i_band, i_gas_band) &
                      * MAX(0.0_RealK, gas_frac_rescaled(l, i, i_gas_band))
                  END DO
                END DO
              ELSE
                DO i=1, atm%n_layer
                  DO l=1, atm%n_profile
                    k_abs_layer(l, i, k, j) &
                      = spectrum%gas%k(k, i_band, i_gas_band) &
                      * MAX(0.0_RealK, gas_frac_rescaled(l, i, i_gas_band) &
                      * atm%gas_mix_ratio(l, i, i_gas_band))
                  END DO
                END DO
              END IF
            END DO
          END IF
        END IF

!       Add k-distribution of continua that are perfectly correlated with the
!       k-distribution of this gas.
        IF (l_cont_band) THEN
          DO j_cont=1, n_cont
            i_cont_band=spectrum%contgen%index_cont(j_cont, i_band)
            IF (spectrum%contgen%i_cont_overlap_band(i_band, i_cont_band) ==   &
                i_gas_band) THEN
              l_cont_added(i_cont_band)=.TRUE.
              i_gas_band_1=spectrum%contgen%index_cont_gas_1(i_cont_band)
              i_gas_band_2=spectrum%contgen%index_cont_gas_2(i_cont_band)
              DO k=1, spectrum%contgen%i_band_k_cont(i_band, i_cont_band)
                DO i=1, atm%n_layer
                  DO l=1, atm%n_profile
                    k_cont_layer(l,i,k)                                        &
                      =k_esft_cont_layer(l,i,k,i_cont_band)*atm%density(l,i)   &
                      *atm%gas_mix_ratio(l,i,i_gas_band_1)                     &
                      *atm%gas_mix_ratio(l,i,i_gas_band_2)
                  END DO
                END DO
              END DO
              CALL sum_k(ierr, atm%n_profile, atm%n_layer, n_abs_esft(j)       &
                , k_abs_layer(1, 1, 1, j), w_abs_esft(1, j)                    &
                , spectrum%contgen%i_band_k_cont(i_band, i_cont_band)          &
                , k_cont_layer                                                 &
                , spectrum%contgen%w_cont(1, i_band, i_cont_band)              &
                , dimen%nd_profile, dimen%nd_layer, nd_k_term                  &
                , spectrum%dim%nd_k_term_cont)
            END IF
          END DO
        END IF
      END DO
    END IF

!   Get generalised continuum absorption data for this band.
    IF (l_cont_band) THEN
!     Scale the continuum ESFT/k-terms
      i_abs=n_gas
      DO j=1, n_cont
        i_cont_band=spectrum%contgen%index_cont(j, i_band)
        IF (.NOT. l_cont_added(i_cont_band)) THEN
!         The overlap treatment for the continuum is the same as that
!         for gases.
          i_abs = i_abs + 1
          index_abs(i_abs)=i_abs
          n_abs_esft(i_abs)=spectrum%contgen%i_band_k_cont(i_band, i_cont_band)
          i_scatter_method_term(1:n_abs_esft(i_abs), i_abs)                    &
            =spectrum%contgen%i_scat_cont(1:n_abs_esft(i_abs), i_band,         &
                                          i_cont_band)
          w_abs_esft(1:n_abs_esft(i_abs), i_abs)                               &
            =spectrum%contgen%w_cont(1:n_abs_esft(i_abs), i_band, i_cont_band)
          i_gas_band_1=spectrum%contgen%index_cont_gas_1(i_cont_band)
          i_gas_band_2=spectrum%contgen%index_cont_gas_2(i_cont_band)
          DO k=1, n_abs_esft(i_abs)
            DO i=1, atm%n_layer
              DO l=1, atm%n_profile
                k_abs_layer(l,i,k,i_abs)                                       &
                  =k_esft_cont_layer(l,i,k,i_cont_band)*atm%density(l,i)       &
                  *atm%gas_mix_ratio(l,i,i_gas_band_1)                         &
                  *atm%gas_mix_ratio(l,i,i_gas_band_2)
              END DO
            END DO
          END DO
          IF (j ==1 .AND. spectrum%contgen%l_cont_major(i_band)) THEN
!           The continuum is the major absorber in this band, swap the index of
!           continuum 1 with the index of gas 1.
            major_abs=index_abs(i_abs)
            index_abs(i_abs)=index_abs(1)
            index_abs(1)=major_abs
          END IF
        END IF
      END DO
    END IF

    IF (i_gas_overlap == ip_overlap_single) THEN
      ! Only the selected gas is active in the band. If different, the major
      ! gas is also included but with zero absorption (set above).
      n_gas=0
      DO i_abs=1, n_abs
        i = index_abs(i_abs)
        i_gas_band = spectrum%gas%index_absorb(i, i_band)
        IF (i <= spectrum%gas%n_band_absorb(i_band)) THEN
          IF (i == 1 .OR. &
              spectrum%gas%type_absorb(i_gas_band) == control%i_gas) THEN
            n_gas=n_gas+1
            index_abs(n_gas)=i
          END IF
        END IF
      END DO
      n_abs=n_gas
      l_abs_band=(n_abs > 0)
    END IF

!   Allocate space for the Planckian emission fields
    CALL allocate_planck(planck, dimen)

    IF (control%isolir == ip_infra_red) THEN
!     Calculate the change in the thermal source function
!     across each layer for the infra-red part of the spectrum.
      CALL diff_planck_source(control, dimen, spectrum, atm, bound,          &
        i_band, i_rad_layer, frac_rad_layer, planck)
    END IF


!   Call a solver appropriate to the presence of gases and
!   the overlap assumed:

    IF (.NOT.l_abs_band) THEN

!     There is no gaseous absorption. Solve for the radiances directly.

      CALL solve_band_without_gas(ierr                                         &
        , control, dimen, spectrum, atm, cld, bound, radout, i_band            &
!                 Atmospheric properties
        , atm%n_profile, atm%n_layer, atm%mass                                 &
!                 Angular integration
        , control%i_angular_integration, control%i_2stream                     &
        , n_order_phase, control%l_rescale, control%n_order_gauss              &
        , control%ms_min, control%ms_max, control%i_truncation                 &
        , ls_local_trunc                                                       &
        , control%accuracy_adaptive, control%euler_factor                      &
        , control%i_sph_algorithm, control%i_sph_mode                          &
!                 Precalculated angular arrays
        , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                             &
!                 Treatment of scattering
        , control%i_scatter_method_band(i_band)                                &
!                 Options for solver
        , control%i_solver                                                     &
!                 Spectral region
        , control%isolir                                                       &
!                 Solar properties
        , bound%zen_0, solar_irrad_band, sph                                   &
!                 Infra-red properties
        , planck                                                               &
!                 Surface properties
        , control%ls_brdf_trunc, bound%n_brdf_basis_fnc                        &
        , bound%rho_alb(1, 1, i_band)                                          &
        , bound%f_brdf, brdf_sol, brdf_hemi                                    &
!                 Tiling of the surface
        , control%l_tile, bound%n_point_tile, bound%n_tile, bound%list_tile    &
        , bound%rho_alb_tile(1, 1, 1, i_band)                                  &
!                 Optical Properties
        , ss_prop, photol                                                      &
!                 Cloudy properties
        , control%l_cloud, control%i_cloud                                     &
!                 Cloudy geometry
        , n_cloud_top                                                          &
        , n_region, k_clr, i_region_cloud, frac_region                         &
        , w_free, cloud_overlap                                                &
        , n_column_slv, list_column_slv                                        &
        , i_clm_lyr_chn, i_clm_cld_typ, area_column                            &
!                 Levels for calculating radiances
        , atm%n_viewing_level, i_rad_layer, frac_rad_layer                     &
!                 Viewing Geometry
        , atm%n_direction, atm%direction                                       &
!                 Radiances
        , i_direct                                                             &
!                 Flags for initialising diagnostics
        , l_initial, l_initial_band                                            &
        , l_initial_channel, l_initial_channel_tile                            &
!                 Flags for fluxes
        , l_actinic, l_clear_band, control%i_solver_clear                      &
!                 Dimensions of arrays
        , dimen%nd_profile, dimen%nd_layer, dimen%nd_layer_clr                 &
        , dimen%id_cloud_top, dimen%nd_column                                  &
        , dimen%nd_flux_profile, dimen%nd_radiance_profile, dimen%nd_j_profile &
        , dimen%nd_cloud_type, dimen%nd_region, dimen%nd_overlap_coeff         &
        , dimen%nd_max_order, dimen%nd_sph_coeff                               &
        , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc                         &
        , dimen%nd_viewing_level, dimen%nd_direction                           &
        , dimen%nd_source_coeff, dimen%nd_point_tile, dimen%nd_tile            &
        )


    ELSE IF (l_abs_band) THEN

!     Gases are included.

!     Treat the gaseous overlaps as directed by the overlap switch.
      SELECT CASE (i_gas_overlap)

      CASE (ip_overlap_random, ip_overlap_exact_major, ip_overlap_single)
        CALL solve_band_random_overlap(ierr                                    &
          , control, dimen, spectrum, atm, cld, bound, radout, i_band          &
!                 Atmospheric properties
          , atm%n_profile, atm%n_layer, atm%mass                               &
!                 Angular integration
          , control%i_angular_integration, control%i_2stream                   &
          , n_order_phase, control%l_rescale, control%n_order_gauss            &
          , control%ms_min, control%ms_max, control%i_truncation               &
          , ls_local_trunc                                                     &
          , control%accuracy_adaptive, control%euler_factor                    &
          , control%i_sph_algorithm, control%i_sph_mode                        &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                           &
!                 Treatment of scattering
          , control%i_scatter_method_band(i_band)                              &
!                 Options for solver
          , control%i_solver, i_gas_overlap                                    &
!                 Gaseous properties
          , n_abs, index_abs, n_abs_esft                                       &
          , k_abs_layer, w_abs_esft                                            &
!                 Spectral region
          , control%isolir                                                     &
!                 Solar properties
          , bound%zen_0, solar_irrad_band, sph                                 &
!                 Infra-red properties
          , planck                                                             &
!                 Surface properties
          , control%ls_brdf_trunc, bound%n_brdf_basis_fnc                      &
          , bound%rho_alb(1, 1, i_band)                                        &
          , bound%f_brdf, brdf_sol, brdf_hemi                                  &
!                 Tiling of the surface
          , control%l_tile, bound%n_point_tile, bound%n_tile, bound%list_tile  &
          , bound%rho_alb_tile(1, 1, 1, i_band)                                &
!                 Optical Properties
          , ss_prop, photol, l_photol_only                                     &
!                 Cloudy properties
          , control%l_cloud, control%i_cloud                                   &
!                 Cloud geometry
          , n_cloud_top                                                        &
          , n_region, k_clr, i_region_cloud, frac_region                       &
          , w_free, cloud_overlap                                              &
          , n_column_slv, list_column_slv                                      &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                          &
!                 Additional variables required for mcica
          , l_cloud_cmp, n_cloud_profile, i_cloud_profile                      &
          , i_cloud_type, dimen%nd_cloud_component                             &
          , control%i_cloud_representation                                     &
!                 Levels for calculating radiances
          , atm%n_viewing_level, i_rad_layer, frac_rad_layer                   &
!                 Viewing Geometry
          , atm%n_direction, atm%direction                                     &
!                 Radiances
          , i_direct                                                           &
!                 Flags for initialising diagnostics
          , l_initial, l_initial_band                                          &
          , l_initial_channel, l_initial_channel_tile                          &
!                 Flags for flux calculations
          , l_actinic, l_clear_band, control%i_solver_clear                    &
!                 Dimensions of arrays
          , dimen%nd_profile, dimen%nd_layer, dimen%nd_layer_clr               &
          , dimen%id_cloud_top, dimen%nd_column                                &
          , dimen%nd_flux_profile, dimen%nd_radiance_profile                   &
          , dimen%nd_j_profile                                                 &
          , nd_abs, nd_k_term, dimen%nd_k_term_inner                           &
          , dimen%nd_cloud_type, dimen%nd_region, dimen%nd_overlap_coeff       &
          , dimen%nd_max_order, dimen%nd_sph_coeff                             &
          , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc                       &
          , dimen%nd_viewing_level, dimen%nd_direction                         &
          , dimen%nd_source_coeff, dimen%nd_point_tile, dimen%nd_tile          &
          )

      CASE (ip_overlap_random_resort_rebin)
!       Set maximum number of ESFT terms needed
        nd_esft_max = MAX(control%n_esft_red, nd_k_term)
        CALL solve_band_random_overlap_resort_rebin(ierr                       &
          , control, dimen, spectrum, atm, cld, bound, radout, i_band          &
!                 Atmospheric properties
          , atm%n_profile, atm%n_layer, atm%mass                               &
!                 Angular integration
          , control%i_angular_integration, control%i_2stream                   &
          , n_order_phase, control%l_rescale, control%n_order_gauss            &
          , control%ms_min, control%ms_max, control%i_truncation               &
          , ls_local_trunc                                                     &
          , control%accuracy_adaptive, control%euler_factor                    &
          , control%i_sph_algorithm, control%i_sph_mode                        &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                           &
!                 Treatment of scattering
          , control%i_scatter_method_band(i_band)                              &
!                 Options for solver
          , control%i_solver                                                   &
          , control%n_esft_red, control%gpnt_split                             &
!                 Gaseous properties
          , n_abs, index_abs, n_abs_esft                                       &
          , k_abs_layer, w_abs_esft                                            &
!                 Spectral region
          , control%isolir                                                     &
!                 Solar properties
          , bound%zen_0, solar_irrad_band, sph                                 &
!                 Infra-red properties
          , planck                                                             &
!                 Surface properties
          , control%ls_brdf_trunc, bound%n_brdf_basis_fnc                      &
          , bound%rho_alb(1, 1, i_band)                                        &
          , bound%f_brdf, brdf_sol, brdf_hemi                                  &
!                 Tiling of the surface
          , control%l_tile, bound%n_point_tile, bound%n_tile, bound%list_tile  &
          , bound%rho_alb_tile(1, 1, 1, i_band)                                &
!                 Optical Properties
          , ss_prop, photol                                                    &
!                 Cloudy properties
          , control%l_cloud, control%i_cloud                                   &
!                 Cloud geometry
          , n_cloud_top                                                        &
          , n_region, k_clr, i_region_cloud, frac_region                       &
          , w_free, cloud_overlap                                              &
          , n_column_slv, list_column_slv                                      &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                          &
!                 Levels for calculating radiances
          , atm%n_viewing_level, i_rad_layer, frac_rad_layer                   &
!                 Viewing Geometry
          , atm%n_direction, atm%direction                                     &
!                 Radiances
          , i_direct                                                           &
!                 Flags for initialising diagnostics
          , l_initial, l_initial_band                                          &
          , l_initial_channel, l_initial_channel_tile                          &
!                 Flags for flux calculations
          , l_actinic, l_clear_band, control%i_solver_clear                    &
!                 Dimensions of arrays
          , dimen%nd_profile, dimen%nd_layer, dimen%nd_layer_clr               &
          , dimen%id_cloud_top, dimen%nd_column                                &
          , dimen%nd_flux_profile, dimen%nd_radiance_profile                   &
          , dimen%nd_j_profile                                                 &
          , nd_abs, nd_k_term, nd_esft_max                                     &
          , dimen%nd_cloud_type, dimen%nd_region, dimen%nd_overlap_coeff       &
          , dimen%nd_max_order, dimen%nd_sph_coeff                             &
          , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc                       &
          , dimen%nd_viewing_level, dimen%nd_direction                         &
          , dimen%nd_source_coeff, dimen%nd_point_tile, dimen%nd_tile          &
          )

      CASE (ip_overlap_k_eqv_scl, ip_overlap_k_eqv, ip_overlap_k_eqv_mod)
        CALL solve_band_k_eqv_scl(ierr                                         &
          , control, dimen, spectrum, atm, cld, bound, radout                  &
!                 Atmospheric properties
          , atm%n_profile, atm%n_layer, atm%mass                               &
!                 Angular integration
          , control%i_angular_integration, control%i_2stream                   &
          , n_order_phase, control%l_rescale, control%n_order_gauss            &
          , control%ms_min, control%ms_max, control%i_truncation               &
          , ls_local_trunc                                                     &
          , control%accuracy_adaptive, control%euler_factor                    &
          , control%i_sph_algorithm, control%i_sph_mode                        &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                           &
!                 Treatment of scattering
          , control%i_scatter_method_band(i_band), i_scatter_method_term       &
!                 Options for solver
          , control%i_solver                                                   &
!                 Gaseous properties
          , i_band, n_abs, index_abs, n_abs_esft, k_abs_layer, w_abs_esft      &
!                 Spectral region
          , control%isolir                                                     &
!                 Solar properties
          , bound%zen_0, solar_irrad_band, sph                                 &
!                 Infra-red properties
          , planck                                                             &
!                 Surface properties
          , control%ls_brdf_trunc, bound%n_brdf_basis_fnc                      &
          , bound%rho_alb(1, 1, i_band)                                        &
          , bound%f_brdf, brdf_sol, brdf_hemi                                  &
          , diffuse_alb_basis                                                  &
!                 Tiling of the surface
          , control%l_tile, bound%n_point_tile, bound%n_tile, bound%list_tile  &
          , bound%rho_alb_tile(1, 1, 1, i_band)                                &
!                 Optical Properties
          , ss_prop, photol, l_photol_only                                     &
!                 Cloudy properties
          , control%l_cloud, control%i_cloud                                   &
!                 Cloud geometry
          , n_cloud_top                                                        &
          , n_region, k_clr, i_region_cloud, frac_region                       &
          , w_free, cloud_overlap                                              &
          , n_column_slv, list_column_slv                                      &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                          &
!                 Additional variables required for mcica
          , l_cloud_cmp, n_cloud_profile, i_cloud_profile                      &
          , i_cloud_type, dimen%nd_cloud_component                             &
          , control%i_cloud_representation                                     &
!                 Levels for calculating radiances
          , atm%n_viewing_level, i_rad_layer, frac_rad_layer                   &
!                 Viewing Geometry
          , atm%n_direction, atm%direction                                     &
!                 Radiances
          , i_direct                                                           &
!                 Flags for initialising diagnostics
          , l_initial, l_initial_band                                          &
          , l_initial_channel, l_initial_channel_tile                          &
!                 Flags for flux calculations
          , l_actinic, l_clear_band, control%i_solver_clear                    &
!                 Dimensions of arrays
          , dimen%nd_profile, dimen%nd_layer, dimen%nd_layer_clr               &
          , dimen%id_cloud_top, dimen%nd_column                                &
          , dimen%nd_flux_profile, dimen%nd_radiance_profile                   &
          , dimen%nd_j_profile                                                 &
          , nd_abs, nd_k_term                                                  &
          , dimen%nd_cloud_type, dimen%nd_region, dimen%nd_overlap_coeff       &
          , dimen%nd_max_order, dimen%nd_sph_coeff                             &
          , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc                       &
          , dimen%nd_viewing_level, dimen%nd_direction                         &
          , dimen%nd_source_coeff, dimen%nd_point_tile, dimen%nd_tile          &
          )

      CASE (ip_overlap_mix_ses2)
        CALL solve_band_ses(ierr                                               &
          , control, dimen, spectrum, atm, cld, bound, radout                  &
!                 Atmospheric properties
          , atm%n_profile, atm%n_layer, atm%mass                               &
!                 Angular integration
          , control%i_angular_integration, control%i_2stream                   &
          , n_order_phase, control%l_rescale, control%n_order_gauss            &
          , control%ms_min, control%ms_max, control%i_truncation               &
          , ls_local_trunc                                                     &
          , control%accuracy_adaptive, control%euler_factor                    &
          , control%i_sph_algorithm, control%i_sph_mode                        &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                           &
!                 Treatment of scattering
          , control%i_scatter_method_band(i_band)                              &
!                 Options for solver
          , control%i_solver                                                   &
!                 Gaseous properties
          , i_band, n_gas, spectrum%gas%index_absorb                           &
          , spectrum%gas%i_band_k_ses, k_esft_layer, spectrum%gas%w_ses        &
          , k_mix_gas_layer, spectrum%gas%n_mix_gas(i_band)                    &
          , atm%gas_mix_ratio, gas_mix_amt                                     &
!                 Continuum absorption
          , k_contm_layer, control%l_continuum                                 &
          , n_continuum, amount_continuum                                      &
!                 Spectral region
          , control%isolir                                                     &
!                 Solar properties
          , bound%zen_0, solar_irrad_band_ses, control%l_solar_tail_flux, sph  &
!                 Infra-red properties
          , planck                                                             &
!                 Surface properties
          , control%ls_brdf_trunc, bound%n_brdf_basis_fnc                      &
          , bound%rho_alb(1, 1, i_band)                                        &
          , bound%f_brdf, brdf_sol, brdf_hemi                                  &
!                 Tiling of the surface
          , control%l_tile, bound%n_point_tile, bound%n_tile, bound%list_tile  &
          , bound%rho_alb_tile(1, 1, 1, i_band)                                &
!                 Optical Properties
          , ss_prop, photol                                                    &
!                 Cloudy properties
          , control%l_cloud, control%i_cloud                                   &
!                 Cloud geometry
          , n_cloud_top                                                        &
          , n_region, k_clr, i_region_cloud, frac_region                       &
          , w_free, cloud_overlap                                              &
          , n_column_slv, list_column_slv                                      &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                          &
!                 Additional variables required for mcica
          , l_cloud_cmp, n_cloud_profile, i_cloud_profile                      &
          , i_cloud_type, dimen%nd_cloud_component                             &
          , control%i_cloud_representation                                     &
!                 Levels for calculating radiances
          , atm%n_viewing_level, i_rad_layer, frac_rad_layer                   &
!                 Viewing Geometry
          , atm%n_direction, atm%direction                                     &
!                 Radiances
          , i_direct                                                           &
!                 Flags for initialising diagnostics
          , l_initial, l_initial_band                                          &
          , l_initial_channel, l_initial_channel_tile                          &
!                 Flags for flux calculations
          , l_actinic, l_clear_band, control%i_solver_clear                    &
!                 Dimensions of arrays
          , dimen%nd_profile, dimen%nd_layer, dimen%nd_layer_clr               &
          , dimen%id_cloud_top, dimen%nd_column                                &
          , dimen%nd_flux_profile, dimen%nd_radiance_profile                   &
          , dimen%nd_j_profile                                                 &
          , spectrum%dim%nd_band, spectrum%dim%nd_species                      &
          , spectrum%dim%nd_k_term, spectrum%dim%nd_continuum                  &
          , dimen%nd_cloud_type, dimen%nd_region, dimen%nd_overlap_coeff       &
          , dimen%nd_max_order, dimen%nd_sph_coeff                             &
          , dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc                       &
          , dimen%nd_viewing_level, dimen%nd_direction                         &
          , dimen%nd_source_coeff, dimen%nd_point_tile, dimen%nd_tile          &
          )

      CASE DEFAULT
        cmessage = '*** Error: An appropriate gaseous overlap ' //             &
          'has not been specified, even though gaseous ' //                    &
          'absorption is to be included.'
        ierr=i_err_fatal
        GO TO 9999
      END SELECT

    END IF

!   Spectral diagnostics for clouds:
!     The extinction and absorptivity calculated in the band are mean values
!     weighted with the fractions of individual types of cloud (which sum to 1).
!     The extinction is then weighted with the clear-sky direct solar flux in
!     the band at the top of the current layer and the total amount of cloud
!     in the grid-box. This definition has the advantage of convenience, but
!     there appears to be no optimal definition of an average extinction.
!     The absorptivity is weighted with the modulus of the clear_sky
!     differential flux in the band at the top of the current layer and the
!     total amount of cloud in the grid-box, as the diagnostic is a measure
!     of the effect of introducing an infinitesimal layer at the top of the
!     current layer into a clear atmosphere on the upward flux at the top of
!     the cloud.
    IF (control%l_cloud_extinction) THEN
      IF (control%l_spherical_solar) THEN
        DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
             l=i_cloud_profile(ll, i)
             radout%cloud_weight_extinction(l, i)                              &
                =radout%cloud_weight_extinction(l, i) + cld%w_cloud(l, i)      &
                *radout%flux_direct_clear_sph_band(l, i, i_band)
             radout%cloud_extinction(l, i)                                     &
                =radout%cloud_extinction(l, i) + cld%w_cloud(l, i)             &
                *radout%flux_direct_clear_sph_band(l, i, i_band)               &
                *cloud_extinction_band(l, i)
          END DO
        END DO
      ELSE
        DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
             l=i_cloud_profile(ll, i)
             radout%cloud_weight_extinction(l, i)                              &
                =radout%cloud_weight_extinction(l, i) + cld%w_cloud(l, i)      &
                *radout%flux_direct_clear_band(l, i-1, i_band)
             radout%cloud_extinction(l, i)                                     &
                =radout%cloud_extinction(l, i) + cld%w_cloud(l, i)             &
                *radout%flux_direct_clear_band(l, i-1, i_band)                 &
                *cloud_extinction_band(l, i)
          END DO
        END DO
      END IF
    END IF

    IF (control%l_cloud_absorptivity) THEN
      DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           radout%cloud_weight_absorptivity(l, i)                              &
              =radout%cloud_weight_absorptivity(l, i) + cld%w_cloud(l, i)      &
              *ABS(radout%flux_up_clear_band(l, i-1, i_band))
           radout%cloud_absorptivity(l, i)                                     &
              =radout%cloud_absorptivity(l, i) + cld%w_cloud(l, i)             &
              *ABS(radout%flux_up_clear_band(l, i-1, i_band))                  &
              *cloud_absorptivity_band(l, i)
        END DO
      END DO
    END IF

    IF (control%l_ls_cloud_extinction) THEN
      IF (control%l_spherical_solar) THEN
        DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
             l=i_cloud_profile(ll, i)
             radout%ls_cloud_weight_extinction(l, i)                           &
                =radout%ls_cloud_weight_extinction(l, i)+cld%w_cloud(l, i)     &
                *(cld%frac_cloud(l,i,1)+cld%frac_cloud(l,i,2))                 &
                *radout%flux_direct_clear_sph_band(l, i, i_band)
             radout%ls_cloud_extinction(l, i)                                  &
                =radout%ls_cloud_extinction(l, i)+cld%w_cloud(l, i)            &
                *(cld%frac_cloud(l,i,1)+cld%frac_cloud(l,i,2))                 &
                *radout%flux_direct_clear_sph_band(l, i, i_band)               &
                *ls_cloud_extinction_band(l, i)
          END DO
        END DO
      ELSE
        DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
             l=i_cloud_profile(ll, i)
             radout%ls_cloud_weight_extinction(l, i)                           &
                =radout%ls_cloud_weight_extinction(l, i)+cld%w_cloud(l, i)     &
                *(cld%frac_cloud(l,i,1)+cld%frac_cloud(l,i,2))                 &
                *radout%flux_direct_clear_band(l, i-1, i_band)
             radout%ls_cloud_extinction(l, i)                                  &
                =radout%ls_cloud_extinction(l, i)+cld%w_cloud(l, i)            &
                *(cld%frac_cloud(l,i,1)+cld%frac_cloud(l,i,2))                 &
                *radout%flux_direct_clear_band(l, i-1, i_band)                 &
                *ls_cloud_extinction_band(l, i)
          END DO
        END DO
      END IF
    END IF

    IF (control%l_ls_cloud_absorptivity) THEN
      DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           radout%ls_cloud_weight_absorptivity(l, i)                           &
              =radout%ls_cloud_weight_absorptivity(l, i)                       &
              +cld%w_cloud(l, i)*(cld%frac_cloud(l,i,1)+cld%frac_cloud(l,i,2)) &
              *ABS(radout%flux_up_clear_band(l, i-1, i_band))
           radout%ls_cloud_absorptivity(l, i)                                  &
              =radout%ls_cloud_absorptivity(l, i)                              &
              +cld%w_cloud(l, i)*(cld%frac_cloud(l,i,1)+cld%frac_cloud(l,i,2)) &
              *ABS(radout%flux_up_clear_band(l, i-1, i_band))                  &
              *ls_cloud_absorptivity_band(l, i)
        END DO
      END DO
    END IF

    IF (control%l_cnv_cloud_extinction) THEN
      IF (control%l_spherical_solar) THEN
        DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
             l=i_cloud_profile(ll, i)
             radout%cnv_cloud_weight_extinction(l, i)                          &
                =radout%cnv_cloud_weight_extinction(l, i)+cld%w_cloud(l, i)    &
                *(cld%frac_cloud(l,i,3)+cld%frac_cloud(l,i,4))                 &
                *radout%flux_direct_clear_sph_band(l, i, i_band)
             radout%cnv_cloud_extinction(l, i)                                 &
                =radout%cnv_cloud_extinction(l, i)+cld%w_cloud(l, i)           &
                *(cld%frac_cloud(l,i,3)+cld%frac_cloud(l,i,4))                 &
                *radout%flux_direct_clear_sph_band(l, i, i_band)               &
                *cnv_cloud_extinction_band(l, i)
          END DO
        END DO
      ELSE
        DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
             l=i_cloud_profile(ll, i)
             radout%cnv_cloud_weight_extinction(l, i)                          &
                =radout%cnv_cloud_weight_extinction(l, i)+cld%w_cloud(l, i)    &
                *(cld%frac_cloud(l,i,3)+cld%frac_cloud(l,i,4))                 &
                *radout%flux_direct_clear_band(l, i-1, i_band)
             radout%cnv_cloud_extinction(l, i)                                 &
                =radout%cnv_cloud_extinction(l, i)+cld%w_cloud(l, i)           &
                *(cld%frac_cloud(l,i,3)+cld%frac_cloud(l,i,4))                 &
                *radout%flux_direct_clear_band(l, i-1, i_band)                 &
                *cnv_cloud_extinction_band(l, i)
          END DO
        END DO
      END IF
    END IF

    IF (control%l_cnv_cloud_absorptivity) THEN
      DO i=n_cloud_top, atm%n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           radout%cnv_cloud_weight_absorptivity(l, i)                          &
              =radout%cnv_cloud_weight_absorptivity(l, i)                      &
              +cld%w_cloud(l, i)*(cld%frac_cloud(l,i,3)+cld%frac_cloud(l,i,4)) &
              *ABS(radout%flux_up_clear_band(l, i-1, i_band))
           radout%cnv_cloud_absorptivity(l, i)                                 &
              =radout%cnv_cloud_absorptivity(l, i)                             &
              +cld%w_cloud(l, i)*(cld%frac_cloud(l,i,3)+cld%frac_cloud(l,i,4)) &
              *ABS(radout%flux_up_clear_band(l, i-1, i_band))                  &
              *cnv_cloud_absorptivity_band(l, i)
        END DO
      END DO
    END IF

!   Deallocate the single scattering propeties.
    CALL deallocate_ss_prop(ss_prop)

!   Deallocate absorber arrays
    DEALLOCATE(l_photol_only)
    DEALLOCATE(l_cont_added)
    DEALLOCATE(i_scatter_method_term)
    DEALLOCATE(n_abs_esft)
    DEALLOCATE(index_abs)
    DEALLOCATE(w_abs_esft)
    DEALLOCATE(k_cont_layer)
    DEALLOCATE(k_abs_layer)


!   Make any adjustments to fluxes and radiances to convert to actual values.
!   This is done inside the loop over bands to allow for division of the
!   output fluxes between separate diagnostic bands.
    IF (control%isolir == ip_infra_red) THEN
      CALL adjust_ir_radiance(control, spectrum, atm, radout, &
        planck, i_band, l_clear_band)
    END IF

!   Deallocate Planckian emission fields
    CALL deallocate_planck(planck)

  END DO ! i_band

  CALL deallocate_sph(sph)
  DEALLOCATE(photol)
  DEALLOCATE(wt_gas)
  DEALLOCATE(jt_gas)

  9999 CONTINUE
  IF (ierr /= i_normal) THEN
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE radiance_calc

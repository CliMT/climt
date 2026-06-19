! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Radiation Program Configuration File

MODULE rad_pcf

! Description:
!   Module containing settings of parameters concerned with the operation
!   of the radiation scheme. This replaces a number of separate program
!   configuration files. The original file names are indicated at the
!   head of each section.

IMPLICIT NONE

! aerosol_parametrization_pcf, aerprm3a
! ------------------------------------------------------------------
! Module to set the parametrizations available for aerosols.
INTEGER, PARAMETER :: ip_aerosol_param_dry       = 1
!   Parametrization for dry aerosols
INTEGER, PARAMETER :: ip_aerosol_param_moist     = 2
!   Parametrization for moist aerosols
INTEGER, PARAMETER :: ip_aerosol_unparametrized  = 3
!   Observational aerosol data
INTEGER, PARAMETER :: ip_aerosol_param_phf_dry   = 4
!   Parametrization of the phase function for dry aerosols
INTEGER, PARAMETER :: ip_aerosol_param_phf_moist = 5
!   Parametrization of the phase function for moist aerosols

! ------------------------------------------------------------------
! aerosol_component_pcf, aercmp3a
! ------------------------------------------------------------------
! Module to set indices of aerosol components.
INTEGER, PARAMETER :: npd_aerosol_component = 32
!   Size allocated for identifiers for aerosols. N.B: this must
!   be at least as large as the largest value in the list below.

! SRA Climatological Aerosols
INTEGER, PARAMETER :: ip_water_soluble      = 1
!   Water soluble aerosol
INTEGER, PARAMETER :: ip_dust_like          = 2
!   Dust-like aerosol
INTEGER, PARAMETER :: ip_oceanic            = 3
!   Oceanic aerosol
INTEGER, PARAMETER :: ip_soot               = 4
!   Soot aerosol
INTEGER, PARAMETER :: ip_ash                = 5
!   Volcanic ash
INTEGER, PARAMETER :: ip_sulphuric          = 6
!   Sulphuric acid

INTEGER, PARAMETER :: ip_ammonium_sulphate  = 7
!   Generic ammonium sulphate aerosol
INTEGER, PARAMETER :: ip_aerosol_uncharacterized = 8
!   Uncharacterized aerosol (for observations)
INTEGER, PARAMETER :: ip_saharan_dust       = 9
!   Saharan dust

! Aerosols for the sulphur cycle
INTEGER, PARAMETER :: ip_accum_sulphate     = 10
!   Accumulation mode sulphate
INTEGER, PARAMETER :: ip_aitken_sulphate    = 11
!   Aitken mode sulphate

! Aerosols for the standard soot model
INTEGER, PARAMETER :: ip_fresh_soot         = 12
!   Fresh soot
INTEGER, PARAMETER :: ip_aged_soot          = 13
!   Aged soot

! Aerosols for sea-salt modelling
INTEGER, PARAMETER :: ip_sodium_chloride    = 14
!   Sodium chloride (generic aerosol)
INTEGER, PARAMETER :: ip_seasalt_film       = 15
!   Sodium chloride (film mode)
INTEGER, PARAMETER :: ip_seasalt_jet        = 16
!   Sodium chloride (jet mode)

! Aerosols for the dust scheme
INTEGER, PARAMETER :: ip_dust_div1          = 17
!   Dust, division 1
INTEGER, PARAMETER :: ip_dust_div2          = 18
!   Dust, division 2
INTEGER, PARAMETER :: ip_dust_div3          = 19
!   Dust, division 3
INTEGER, PARAMETER :: ip_dust_div4          = 20
!   Dust, division 4
INTEGER, PARAMETER :: ip_dust_div5          = 21
!   Dust, division 5
INTEGER, PARAMETER :: ip_dust_div6          = 22
!   Dust, division 6

! Biomass aerosols
INTEGER, PARAMETER :: ip_biomass_1          = 23
!   Biomass (division 1)
INTEGER, PARAMETER :: ip_biomass_2          = 24
!   Biomass (division 2)

! Biogenic aerosol
INTEGER, PARAMETER :: ip_biogenic           = 25

! Fossil-fuel organic carbon aerosols
INTEGER, PARAMETER :: ip_ocff_fresh         = 26
INTEGER, PARAMETER :: ip_ocff_aged          = 27

! "Delta" aerosol
! (unspecified aerosol used to complement NWP climatologies)
INTEGER, PARAMETER :: ip_delta              = 28

! Murk aerosol
INTEGER, PARAMETER :: ip_murk               = 29

! Ammonium nitrate aerosol
INTEGER, PARAMETER :: ip_nitrate            = 30

! Aerosols for the two-bin dust scheme
INTEGER, PARAMETER :: ip_twobindust_1       = 31
INTEGER, PARAMETER :: ip_twobindust_2       = 32

CHARACTER (LEN=20), PARAMETER :: &
  name_aerosol_component(npd_aerosol_component) = (/ &
                                   "Water soluble       ", &
                                   "Dust-like           ", &
                                   "Oceanic             ", &
                                   "Soot                ", &
                                   "Volcanic Ash        ", &
                                   "Sulphuric Acid      ", &
                                   "Ammonium Sulphate   ", &
                                   "Uncharacterized     ", &
                                   "Saharan Dust        ", &
                                   "Accum. Sulphate     ", &
                                   "Aitken Sulphate     ", &
                                   "Fresh Soot          ", &
                                   "Aged Soot           ", &
                                   "Generic NaCl        ", &
                                   "NaCl film mode      ", &
                                   "NaCl jet mode       ", &
                                   "Dust Division 1     ", &
                                   "Dust Division 2     ", &
                                   "Dust Division 3     ", &
                                   "Dust Division 4     ", &
                                   "Dust Division 5     ", &
                                   "Dust Division 6     ", &
                                   "Biomass Division 1  ", &
                                   "Biomass Division 2  ", &
                                   "Biogenic            ", &
                                   "Fresh fossil-fuel OC", &
                                   "Aged fossil-fuel OC ", &
                                   "Delta aerosol       ", &
                                   "Murk                ", &
                                   "Ammonium nitrate    ", &
                                   "Two-bin Dust Div 1  ", &
                                   "Two-bin Dust Div 2  " /)

! Aerosol data can be provided to radiaiton etc from different
! sources, either prognostic, or different climatologies.
! These can either effect radiative fluxes (RON), or are there
! just for diagnostics (ROFF)
! Simple land/sea climatology from Cusack(98):
INTEGER, PARAMETER :: ip_aersrc_cusack_ron   = 0
INTEGER, PARAMETER :: ip_aersrc_cusack_roff  = 10
! CLASSIC prognostic aerosols:
INTEGER, PARAMETER :: ip_aersrc_classic_ron  = 1
INTEGER, PARAMETER :: ip_aersrc_classic_roff = 11
! Aerosol climatolgies, derived from CLASSIC:
INTEGER, PARAMETER :: ip_aersrc_arcl_ron     = 2
INTEGER, PARAMETER :: ip_aersrc_arcl_roff    = 12

! ------------------------------------------------------------------
! angular_integration_pcf, angint3a
! ------------------------------------------------------------------
! Module to set the types of angular integration.
INTEGER, PARAMETER :: ip_two_stream          = 1
!   Two-stream scheme
INTEGER, PARAMETER :: ip_ir_gauss            = 2
!   Gaussian integration in the IR
INTEGER, PARAMETER :: ip_spherical_harmonic  = 3
!   Integration by spherical harmonics

! ------------------------------------------------------------------
! cloud_component_pcf, cldcmp3a
! ------------------------------------------------------------------
! Module to set components of clouds.
INTEGER, PARAMETER :: ip_clcmp_st_water    = 1
!   Stratiform water droplets
INTEGER, PARAMETER :: ip_clcmp_st_ice      = 2
!   Stratiform ice crystals
INTEGER, PARAMETER :: ip_clcmp_cnv_water   = 3
!   Convective water droplets
INTEGER, PARAMETER :: ip_clcmp_cnv_ice     = 4
!   Convective ice crystals

! ------------------------------------------------------------------
! cloud_parametrization_pcf, wclprm3a
! ------------------------------------------------------------------
! Module to set numbers for water cloud schemes.
INTEGER, PARAMETER :: ip_slingo_schrecker    = 1
!   Parametrization of Slingo & Schrecker
INTEGER, PARAMETER :: ip_ackerman_stephens   = 2
!   Parametrization of Ackerman & Stephens
INTEGER, PARAMETER :: ip_drop_unparametrized = 3
!   Unparametrized droplet data
INTEGER, PARAMETER :: ip_drop_pade_2         = 5
!   Pade approximation of the second order
!   (third order for the extinction)
INTEGER, PARAMETER :: ip_slingo_schr_phf     = 6
!   Parameterization of Slingo & Schrecker +
!   Moments of Phase Functions
INTEGER, PARAMETER :: IP_drop_Pade_2_PHF     = 7
!   Pade approximation of the second order (third order for the
!    extinction) extended to higher moments of the phase function
INTEGER, PARAMETER :: IP_ps_size_PHF         = 8
!   Pade approximation of the second order (third order for the
!    extinction) extended to higher moments of the phase function
INTEGER, PARAMETER :: ip_drop_sun_shine_vis  = 10
!   Sun and Shine's parametrization in the visible
INTEGER, PARAMETER :: ip_drop_stamnes        = 11
!   Parametrization of Stamnes

! ------------------------------------------------------------------
! cloud_region_pcf, cldreg3a
! ------------------------------------------------------------------
! Module to define reference numbers for regions of clouds.
INTEGER, PARAMETER :: ip_region_clear = 1
!   Reference number for clear-sky region
INTEGER, PARAMETER :: ip_region_strat = 2
!   Reference number for stratiform cloudy region
INTEGER, PARAMETER :: ip_region_conv  = 3
!   Reference number for convective cloudy region

! ------------------------------------------------------------------
! cloud_representation_pcf, clrepp3a (part)
! ------------------------------------------------------------------
! Module to set representations of clouds.
INTEGER, PARAMETER :: ip_cloud_homogen           = 1
!   Stratiform components are mixed homogeneously
INTEGER, PARAMETER :: ip_cloud_ice_water         = 2
!   Ice and water clouds are treated separately
!   for stratiform cloud
INTEGER, PARAMETER :: ip_cloud_conv_strat        = 3
!   Clouds are divided into homogeneously mixed
!   stratiform and convective parts
INTEGER, PARAMETER :: ip_cloud_csiw              = 4
!   Clouds divided into ice and water phases and
!   into stratiform and convective components.
INTEGER, PARAMETER :: ip_cloud_off               = 5
!   Clear column
INTEGER, PARAMETER :: ip_cloud_combine_homogen   = 6
!   All phases of stratiform and convective
!   components are mixed homogeneously
INTEGER, PARAMETER :: ip_cloud_combine_ice_water = 7
!   Stratiform and convective clouds are combined
!   but divided into ice and water components
INTEGER, PARAMETER :: ip_cloud_split_homogen     = 8
!   Cloud is split into optically thick and thin
!   regions of mixed phase cloud.
INTEGER, PARAMETER :: ip_cloud_split_ice_water   = 9
!   Cloud is split into optically thick and thin
!   regions divided into ice and water components

! ------------------------------------------------------------------
! treatment of cloud entrapment of flux between layers
! ------------------------------------------------------------------
INTEGER, PARAMETER :: ip_max_entrapment  = 1
!   Maximum entrapment
INTEGER, PARAMETER :: ip_zero_entrapment = 2
!   Zero entrapment using 'Hogan' solvers

! ------------------------------------------------------------------
! treatment of in-cloud horizontal water content inhomogeneity
! ------------------------------------------------------------------
INTEGER, PARAMETER :: ip_homogeneous       = 0
!   Treat clouds as horizontally homogeneous
INTEGER, PARAMETER :: ip_scaling           = 1
!   Represent inhomogeneity with a scaling factor
INTEGER, PARAMETER :: ip_mcica             = 2
!   Represent inhomogeneity with the McICA scheme.
INTEGER, PARAMETER :: ip_cairns            = 3
!   Represent inhomogeneity by calculating corrections to
!   single scattering parameters as described in Cairns et al.,
!   Journal of Atmospheric Sciences, Volume 57, 2000.
!   Only supported with two-stream approximation.
INTEGER, PARAMETER :: ip_tripleclouds_2019 = 4
!   Split cloud into optically thick and thin regions using
!   the parametrisation from Hogan et al., 2019
!   DOI:10.1175/JAS-D-18-0366.1

! ------------------------------------------------------------------
! treatment of cloud vertical overlap
! ------------------------------------------------------------------
INTEGER, PARAMETER :: ip_max_rand         = 0
!   Maximum/random overlap
INTEGER, PARAMETER :: ip_rand             = 1
!   Random overlap
INTEGER, PARAMETER :: ip_exponential      = 2
!   Exponential overlap
INTEGER, PARAMETER :: ip_exponential_rand = 3
!   Exponential-random overlap

! ------------------------------------------------------------------
! cloud_scheme_pcf, clschm3a
! ------------------------------------------------------------------
! Module to define reference numbers for cloud schemes.
INTEGER, PARAMETER :: ip_cloud_mix_max         = 2
!   Maximum/random overlap in a mixed column
INTEGER, PARAMETER :: ip_cloud_mix_random      = 4
!   Random overlap in a mixed column
INTEGER, PARAMETER :: ip_cloud_column_max      = 3
!   Maximum overlap in a column model
INTEGER, PARAMETER :: ip_cloud_clear           = 5
!   Clear column
INTEGER, PARAMETER :: ip_cloud_triple          = 6
!   Mixed column with split between convective and layer cloud
INTEGER, PARAMETER :: ip_cloud_part_corr       = 7
!   Coupled overlap with partial correlation of cloud
INTEGER, PARAMETER :: ip_cloud_part_corr_cnv   = 8
!   Coupled overlap with partial correlation of cloud
!   with a separate treatment of convective cloud
INTEGER, PARAMETER :: ip_cloud_mcica           = 10
!   Performs McICA on generated cloud

! ------------------------------------------------------------------
! cloud_type_pcf, cldtyp3a
! ------------------------------------------------------------------
! Module to set types of clouds.
INTEGER, PARAMETER :: ip_cloud_type_homogen = 1
!   Cloud composed of mixed water and ice

INTEGER, PARAMETER :: ip_cloud_type_water   = 1
!   Cloud composed only of water
INTEGER, PARAMETER :: ip_cloud_type_ice     = 2
!   Cloud composed only of ice

INTEGER, PARAMETER :: ip_cloud_type_strat   = 1
!   Mixed-phase stratiform cloud
INTEGER, PARAMETER :: ip_cloud_type_conv    = 2
!   Mixed-phase convective cloud

INTEGER, PARAMETER :: ip_cloud_type_sw      = 1
!   Stratiform water cloud
INTEGER, PARAMETER :: ip_cloud_type_si      = 2
!   Stratiform ice cloud
INTEGER, PARAMETER :: ip_cloud_type_cw      = 3
!   Convective water cloud
INTEGER, PARAMETER :: ip_cloud_type_ci      = 4
!   Convective ice cloud

! ------------------------------------------------------------------
! setting of droplet effective radii
! ------------------------------------------------------------------
INTEGER, PARAMETER :: ip_re_default  = 1
!   Use default effective radius of 7um
INTEGER, PARAMETER :: ip_re_external = 2
!   Use input values of effective radius
INTEGER, PARAMETER :: ip_re_liu      = 3
!   Liu spectral dispersion
INTEGER, PARAMETER :: ip_re_constant = 4
!   Use input scalar value of effective radius

! ------------------------------------------------------------------
! continuum_pcf, cntuum3a
! ------------------------------------------------------------------
! Module setting parameters for continuum data.
INTEGER, PARAMETER :: ip_self_continuum = 1
!   Self-broadened continuum
INTEGER, PARAMETER :: ip_frn_continuum  = 2
!   Foreign-broadened continuum
INTEGER, PARAMETER :: ip_n2_continuum   = 3
!   Nitrogen continuum

! ------------------------------------------------------------------
! continuum_index (SES2)
! ------------------------------------------------------------------
! Module to set the defining numbers for the cotinuum types
INTEGER, PARAMETER :: ip_cont_h2o    = 1
!   h2o continuum index number
INTEGER, PARAMETER :: ip_cont_o4     = 2
!   o2-o2 continuum
INTEGER, PARAMETER :: ip_cont_o2n2_1 = 3
!   o2n2 continuum at 1.06 mu
INTEGER, PARAMETER :: ip_cont_o2n2_2 = 4
!   o2n2 continuum at 1.27 mu

! ------------------------------------------------------------------
! error_pcf, error3a
! ------------------------------------------------------------------
! Module to set error flags in the radiation code.
INTEGER, PARAMETER :: i_normal            = 0
!   Error free condition
INTEGER, PARAMETER :: i_err_fatal         = 1
!   Fatal error: immediate return
INTEGER, PARAMETER :: i_abort_calculation = 2
!   Calculation aborted
INTEGER, PARAMETER :: i_missing_data      = 3
!   Missing data error: conditional
INTEGER, PARAMETER :: i_err_io            = 4
!   I/O error
INTEGER, PARAMETER :: i_err_range         = 5
!   Interpolation range error
INTEGER, PARAMETER :: i_err_exist         = 6
!   Existence error
INTEGER, PARAMETER :: i_warning           = -1
!   Non-fatal warning

! ------------------------------------------------------------------
! k_scale_pcf, esft_scale_pcf, esftsc3a
! ------------------------------------------------------------------
! Module setting types of ESFT scaling
INTEGER, PARAMETER :: ip_scale_null = 0
!   No scaling at all
INTEGER, PARAMETER :: ip_scale_band = 1
!   Same scaling throughout band
INTEGER, PARAMETER :: ip_scale_term = 2
!   Different scaling for each ESFT

! ------------------------------------------------------------------
! gas_overlap_pcf, gasovl3a
! ------------------------------------------------------------------
! Module to set treatments of overlapping gaseous absorption.
INTEGER, PARAMETER :: ip_overlap_hybrid              = 0
!   Overlap treatment depends on band (specified in spectral file)
INTEGER, PARAMETER :: ip_overlap_single              = 1
!   One species only
INTEGER, PARAMETER :: ip_overlap_random              = 2
!   Random overlap
INTEGER, PARAMETER :: ip_overlap_exact_major         = 3
!   Exact overlap between minor and major species,
!   random overlap between minor species
INTEGER, PARAMETER :: ip_overlap_k_eqv_scl           = 4
!   Equivalent extinction using correct scaling
INTEGER, PARAMETER :: ip_overlap_k_eqv               = 5
!   Equivalent extinction
INTEGER, PARAMETER :: ip_overlap_k_eqv_mod           = 6
!   Equivalent extinction using modulus of flux
INTEGER, PARAMETER :: ip_overlap_mix_ses2            = 7
!   Mixed gases used in the SES2 scheme
INTEGER, PARAMETER :: ip_overlap_random_resort_rebin = 8
!   Random overlap with resorting and rebinning

! ------------------------------------------------------------------
! ice_cloud_param_pcf, iclprm3a
! ------------------------------------------------------------------
! Module to set numbers for ice cloud schemes.
INTEGER, PARAMETER :: ip_slingo_schrecker_ice     = 1
!   Parametrization of Slingo and Schrecker.
INTEGER, PARAMETER :: ip_ice_unparametrized       = 3
!   Unparametrized ice crystal data
INTEGER, PARAMETER :: ip_sun_shine_vn2_vis        = 4
!   Sun and Shine's parametrization in the visible (version 2)
INTEGER, PARAMETER :: ip_sun_shine_vn2_ir         = 5
!   Sun and Shine's parametrization in the IR (version 2)
INTEGER, PARAMETER :: ip_ice_adt                  = 6
!   ADT-based scheme for ice crystals
INTEGER, PARAMETER :: ip_ice_adt_10               = 7
!   ADT-based scheme for ice crystals using 10th order polynomials
INTEGER, PARAMETER :: ip_ice_fu_solar             = 9
!   Fu's parametrization in the solar region of the spectrum
INTEGER, PARAMETER :: ip_ice_fu_ir                = 10
!   Fu's parametrization in the infra-red region of the spectrum
INTEGER, PARAMETER :: ip_slingo_schr_ice_phf      = 11
!   Parametrization of Slingo and Schrecker
!   (Moments of phase function).
INTEGER, PARAMETER :: ip_ice_agg_de               = 12
!   Provisional agregate parametrization.
INTEGER, PARAMETER :: ip_ice_fu_phf               = 12
!   Parametrization like Fu (Moments of phase function).
!   ip_ice_agg_de is equivalent to ip_ice_fu_phf with 1 moment.
INTEGER, PARAMETER :: ip_ice_sun_fu               = 14
!   Parametrization of Sun and Fu
INTEGER, PARAMETER :: ip_ice_chou_vis             = 15
!   Parametrization of M-D Chou, et al.
!   (2002, JGR, 107, D21, 10.1029/2002JD002061)
INTEGER, PARAMETER :: ip_ice_agg_de_sun           = 16
!   This scheme is the same as IP_ICE_AGG_DE except that the
!   number of parameters is 11 and is used for SES2
INTEGER, PARAMETER :: ip_ice_t_iwc                = 17
!   Fit single scattering properties directly from T and IWC
INTEGER, PARAMETER :: ip_ice_iwc_only             = 18
!   Fit single scattering properties directly from IWC
INTEGER, PARAMETER :: ip_ice_baran                = 19
!   Optical properties derived from ensemble of ice crystals
!   similar to Baran 2014, with additional temperature dependence. 
INTEGER, PARAMETER :: ip_ice_pade_2_phf           = 20
!   Pade approximation of the second order (third order for the
!   extinction) including higher moments of the phase function

! ------------------------------------------------------------------
! phase_pcf, phase3a
! ------------------------------------------------------------------
! Module to set indices for phases.
INTEGER, PARAMETER :: ip_phase_water = 1
!   Liquid phase
INTEGER, PARAMETER :: ip_phase_ice   = 2
!   Ice phase

! ------------------------------------------------------------------
! scale_fnc_pcf, sclfnc3a
! ------------------------------------------------------------------
! Module to set types of scaling for absorber amounts
INTEGER, PARAMETER :: npd_scale_fnc = 10
!   Number of scaling functions
INTEGER, PARAMETER :: n_scale_variable(0: npd_scale_fnc) = &
  (/ 0, 2, 3, 4, 0, 0, 6, 8, 10, 0, 0 /)
!   Number of parameters in scaling functions

INTEGER, PARAMETER :: ip_scale_fnc_null     = 0
!   Null scaling function
INTEGER, PARAMETER :: ip_scale_power_law    = 1
!   Power law scaling function
INTEGER, PARAMETER :: ip_scale_power_quad   = 2
!   Power law for p; quadratic for T
INTEGER, PARAMETER :: ip_scale_doppler_quad = 3
!   Power law for p; quadratic for T with implicit
!   Doppler correction
INTEGER, PARAMETER :: ip_scale_wenyi        = 4
!   Wenyi scaling for pressure and temperature
INTEGER, PARAMETER :: ip_scale_ses2         = 5
!   SES2 scaling for pressure and temperature
INTEGER, PARAMETER :: ip_scale_dbl_pow_law  = 6
!   Two power law scaling functions either side
!   of the reference pressure/temperature.
INTEGER, PARAMETER :: ip_scale_dbl_pow_quad = 7
!   Power law for p; quadratic for T (two)
INTEGER, PARAMETER :: ip_scale_dbl_dop_quad = 8
!   Power law for p; quadratic for T with implicit
!   Doppler correction (two)
INTEGER, PARAMETER :: ip_scale_lookup       = 9
!   Scaling uses a look-up table in P, T
INTEGER, PARAMETER :: ip_scale_t_lookup     = 10
!   Scaling uses a look-up table in T only

! -----------------------------------------------------------------
! scatter_method_pcf, sctmth3a
! ------------------------------------------------------------------
! Module to set the methods of treating scattering.
INTEGER, PARAMETER :: ip_scatter_full   = 1
!   Full treatment of scattering
INTEGER, PARAMETER :: ip_no_scatter_abs = 2
!   Scattering ignored completely.
INTEGER, PARAMETER :: ip_no_scatter_ext = 3
!   Scattering treated as absorption
INTEGER, PARAMETER :: ip_scatter_approx = 4
!   Approximate treatment of scattering
INTEGER, PARAMETER :: ip_scatter_hybrid = 5
!   Scattering treatment depends on k-term (specified in spectral file)

! -----------------------------------------------------------------
! Scaling optical depth for direct solar flux calculation
! ------------------------------------------------------------------
INTEGER, PARAMETER :: ip_direct_noscaling = 0
!   No scaling optical depth for direct flux calculation
INTEGER, PARAMETER :: ip_direct_delta_scaling = 1
!   Delta-Eddington scaling for direct flux calculation
INTEGER, PARAMETER :: ip_direct_csr_scaling = 2
!   Circumsolar fraction scaling for direct flux calculation

! ------------------------------------------------------------------
! solver_pcf, solver3a
! ------------------------------------------------------------------
! Module to define reference numbers for solvers.
INTEGER, PARAMETER :: ip_solver_pentadiagonal    = 1
!   Pentadiagonal scheme
INTEGER, PARAMETER :: ip_solver_no_scat          = 2
!   Direct solver ignoring scattering
INTEGER, PARAMETER :: ip_solver_mix_app_scat     = 9
!   Mixed column scheme with approximate scattering
INTEGER, PARAMETER :: ip_solver_mix_direct       = 11
!   Direct mixed column scheme for full fluxes
INTEGER, PARAMETER :: ip_solver_homogen_direct   = 13
!   Direct solver for a homogeneous column
INTEGER, PARAMETER :: ip_solver_triple           = 14
!   Direct solver for triple column
INTEGER, PARAMETER :: ip_solver_triple_app_scat  = 15
!   Direct solver for triple column approximating scattering
INTEGER, PARAMETER :: ip_solver_mix_direct_hogan = 16
!   Direct mixed column scheme for full fluxes (modified
!   for correct treatment of shadowing by Robin Hogan)
INTEGER, PARAMETER :: ip_solver_triple_hogan     = 17
!   Direct solver for triple column (modified for correct
!   treatment of shadowing by Robin Hogan)

! ------------------------------------------------------------------
! source_coeff_pointer_pcf, scfpt3a
! ------------------------------------------------------------------
! Module to set pointers to source coefficients.
INTEGER, PARAMETER :: ip_scf_solar_up   = 1
!   Pointer to source coeficient for upward solar beam
INTEGER, PARAMETER :: ip_scf_solar_down = 2
!   Pointer to source coeficient for downward solar beam
INTEGER, PARAMETER :: ip_scf_ir_1d      = 1
!   Pointer to source coeficient for 1st difference of planckian
INTEGER, PARAMETER :: ip_scf_ir_2d      = 2
!   Pointer to source coeficient for 2nd difference of planckian

! -----------------------------------------------------------------
! spectral_region_pcf, spcrg3a
! ------------------------------------------------------------------
! Module to set flags for different portions of the spectrum.
INTEGER, PARAMETER :: ip_solar     = 1
!   Solar region
INTEGER, PARAMETER :: ip_infra_red = 2
!   Infra-red region

! ------------------------------------------------------------------
! sph_algorithm_pcf
! ------------------------------------------------------------------
! Module to set the allowed spherical harmonic algorithms.
INTEGER, PARAMETER :: ip_sph_direct       = 1
!   Direct solution using spherical harmonics
INTEGER, PARAMETER :: ip_sph_reduced_iter = 2
!   The spherical harmonic solution a reduced order of
!   truncation is used to define a source
!   term for integration along a line: this can be combined
!   with a higher order of truncation for the solar beam
!   to yield a solution almost identical to the TMS method
!   of Nakajima and Tanaka.

! ------------------------------------------------------------------
! sph_mode_pcf
! ------------------------------------------------------------------
! Module to set modes in which the spherical harmonic algorithm
! can be used.
INTEGER, PARAMETER :: ip_sph_mode_rad  = 1
!   Spherical harmonics are used to calculate radiances
INTEGER, PARAMETER :: ip_sph_mode_flux = 2
!   Spherical harmonics are used to calculate fluxes
INTEGER, PARAMETER :: ip_sph_mode_j    = 3
!   Spherical harmonics are used to calculate mean
!   radiances (actinic flux/4 pi)

! ------------------------------------------------------------------
! sph_truncation_pcf
! ------------------------------------------------------------------
! Module to set the types of spherical truncation.
INTEGER, PARAMETER :: ip_trunc_triangular   = 1
!   Trapezoidal truncation
INTEGER, PARAMETER :: ip_trunc_rhombohedral = 2
!   Rhombohedral truncation
INTEGER, PARAMETER :: ip_trunc_azim_sym     = 3
!   Truncation with azimuthal symmetry
INTEGER, PARAMETER :: ip_trunc_adaptive     = 4
!   Truncation set adaptively

! ------------------------------------------------------------------
! surface_spec_pcf, srfsp3a
! ------------------------------------------------------------------
! Module to set permitted methods of specifying the characteristics
! of the surface.

! Ways of specifiying the surface
INTEGER, PARAMETER :: ip_surface_specified           = 1
!   Properties specified by surface type
INTEGER, PARAMETER :: ip_surface_internal            = 2
!   Properties passed into code
INTEGER, PARAMETER :: ip_surface_polynomial          = 3
!   Direct albedo fitted as polynomial
INTEGER, PARAMETER :: ip_surface_payne               = 4
!   Fit in the functional form used by Payne
INTEGER, PARAMETER :: ip_surface_lambertian          = 5
!   BRDF represented by a Lambertian
INTEGER, PARAMETER :: ip_surface_roujean             = 6
!   BRDF represented by a Roujean's basis
INTEGER, PARAMETER :: ip_surface_lommel_seeliger_axi = 7
!   BRDF represented by an axisymmetric Lommel-Seeliger function

! Pointers to specific components of arrays
INTEGER, PARAMETER :: ip_surf_alb_diff = 1
!   Pointer to diffuse surface albedo
INTEGER, PARAMETER :: ip_surf_alb_dir  = 2
!   Pointer to direct surface albedo

! ------------------------------------------------------------------
! tileid3z
! ------------------------------------------------------------------
! This module defines identifiers for different surface types
! as used in the radiation scheme.
INTEGER, PARAMETER :: npd_tile_type = 4
!   Number identifiers for surface types (highest value in the list below).

INTEGER, PARAMETER :: ip_ocean_tile   = 1
!   Identifier for open sea
INTEGER, PARAMETER :: ip_seaice_tile  = 2
!   Idenitifer for ice
INTEGER, PARAMETER :: ip_land_tile    = 3
!   Identifer for land
INTEGER, PARAMETER :: ip_landice_tile = 4

! ------------------------------------------------------------------
! two_stream_scheme_pcf, twostr3a
! ------------------------------------------------------------------
! Module to set the defining numbers for the two-stream schemes.
INTEGER, PARAMETER :: ip_eddington        = 2
!   Eddington approximation
INTEGER, PARAMETER :: ip_discrete_ord     = 4
!   Discrete ordinate method
INTEGER, PARAMETER :: ip_ifm              = 5
!   Improved flux method
INTEGER, PARAMETER :: ip_pifm85           = 6
!   Practical improved flux method
!   (version of Zdunkowski et al. 1985)
INTEGER, PARAMETER :: ip_zdk_flux         = 7
!   Zdunkowski's flux method
INTEGER, PARAMETER :: ip_krschg_flux      = 8
!   Kerschgen's flux method
INTEGER, PARAMETER :: ip_coakley_chylek_1 = 9
!   Coakley & chylek's 1st method
INTEGER, PARAMETER :: ip_coakley_chylek_2 = 10
!   Coakley & chylek's 2nd method
INTEGER, PARAMETER :: ip_meador_weaver    = 11
!   Meador & weaver's method
INTEGER, PARAMETER :: ip_elsasser         = 12
!   Elsasser's diffusivity scheme
INTEGER, PARAMETER :: ip_2s_test          = 14
!   User's defined test approximation.
INTEGER, PARAMETER :: ip_hemi_mean        = 15
!   Hemispheric mean approximation.
INTEGER, PARAMETER :: ip_pifm80           = 16
!   Practical improved flux method
!   (version of Zdunkowski et al. 1980)

! ------------------------------------------------------------------
! Rayleigh scattering
! ------------------------------------------------------------------
! Module with pointers to Rayleigh scattering coefficient
! calculation scheme.
INTEGER, PARAMETER :: ip_rayleigh_total  = 1
!   Rayleigh scattering coefficients tabulated for total gas
INTEGER, PARAMETER :: ip_rayleigh_custom = 2
!   Rayleigh scattering coefficients tabulated for gases separately

! ------------------------------------------------------------------

END MODULE rad_pcf

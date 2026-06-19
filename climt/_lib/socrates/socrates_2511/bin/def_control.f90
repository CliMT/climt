! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure of control options.
!
! Description:
!   This module defines the elements of the structure defining
!   algorithmic control of the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_control

USE filenamelength_mod, ONLY: filenamelength
USE gas_list_pcf, ONLY: npd_gases  
USE missing_data_mod, ONLY: rmdi, imdi
USE realtype_rd, ONLY: RealK

IMPLICIT NONE

TYPE StrCtrl

! Name of spectral file
  CHARACTER (LEN=filenamelength) :: spectral_file                 = ''


! Spectral region and bands
  INTEGER :: isolir                                               = imdi
!   Spectral region
  INTEGER :: first_band                                           = imdi
!   First band to use in the calculation
  INTEGER :: last_band                                            = imdi
!   Last band to use in the calculation


! Physical processes
  LOGICAL :: l_microphysics                                       = .FALSE.
!   Flag for microphysics
  LOGICAL :: l_gas                                                = .FALSE.
!   Flag for gaseous absorption
  LOGICAL :: l_rayleigh                                           = .FALSE.
!   Flag for Rayleigh scattering
  LOGICAL :: l_continuum                                          = .FALSE.
!   Flag for the continuum
  LOGICAL :: l_cont_gen                                           = .FALSE.
!   Flag for the generalised continua
  LOGICAL :: l_cloud                                              = .FALSE.
!   Flag for clouds
  LOGICAL :: l_drop                                               = .FALSE.
!   Flag for droplets
  LOGICAL :: l_ice                                                = .FALSE.
!   Flag for ice crystals
  LOGICAL :: l_aerosol                                            = .FALSE.
!   Flag for aerosols
  LOGICAL :: l_aerosol_mode                                       = .FALSE.
!   Flag for modal aerosols
  LOGICAL :: l_aerosol_ccn                                        = .FALSE.
!   Flag for aerosols as CCN
  LOGICAL :: l_solar_tail_flux                                    = .FALSE.
!   Flag for adding solar tail flux to LW ragion
  LOGICAL :: l_orog                                               = .FALSE.
!   Correct the direct solar flux at the surface for sloping terrain
  LOGICAL :: l_orog_fix                                           = .FALSE.
!   Fix a bug in the orographic correction for the column solver
  LOGICAL :: l_solvar                                             = .FALSE.
!   Time variation of the solar spectrum

! Gaseous absorption:
  INTEGER :: i_gas_overlap                                        = imdi
!   Treatment of gaseous overlaps
  INTEGER :: n_esft_red                                           = imdi
!   Number of reduced (resorted and rebinned) ESFT terms in
!   each band using random overlap with resorting and rebinning
  REAL (RealK) :: gpnt_split                                      = rmdi
!   g-coordinate point for splitting into two subintervals using
!   random overlap with resorting and rebinning
  INTEGER :: i_gas                                                = imdi
!   Gas to be considered (if only one gas)
  LOGICAL :: l_h2o                                                = .FALSE.
!   Flag for absorption by water vapour
  LOGICAL :: l_o2                                                 = .FALSE.
!   Flag for absorption by oxygen
  LOGICAL :: l_n2o                                                = .FALSE.
!   Flag for absorption by nitrous oxide
  LOGICAL :: l_ch4                                                = .FALSE.
!   Flag for absorption by methane
  LOGICAL :: l_cfc11                                              = .FALSE.
!   Flag for absorption by CFC11
  LOGICAL :: l_cfc12                                              = .FALSE.
!   Flag for absorption by CFC12
  LOGICAL :: l_cfc113                                             = .FALSE.
!   Flag for absorption by CFC113
  LOGICAL :: l_cfc114                                             = .FALSE.
!   Flag for absorption by CFC114
  LOGICAL :: l_hcfc22                                             = .FALSE.
!   Flag for absorption by HCFC22
  LOGICAL :: l_hfc125                                             = .FALSE.
!   Flag for absorption by HFC125
  LOGICAL :: l_hfc134a                                            = .FALSE.
!   Flag for absorption by HFC134A
  LOGICAL :: l_co                                                 = .FALSE.
!   Flag for absorption by carbon monoxide
  LOGICAL :: l_nh3                                                = .FALSE.
!   Flag for absorption by ammonia
  LOGICAL :: l_tio                                                = .FALSE.
!   Flag for absorption by titanium oxide
  LOGICAL :: l_vo                                                 = .FALSE.
!   Flag for absorption by vanadium oxide
  LOGICAL :: l_h2                                                 = .FALSE.
!   Flag for absorption by hydrogen
  LOGICAL :: l_he                                                 = .FALSE.
!   Flag for absorption by helium
  LOGICAL :: l_na                                                 = .FALSE.
!   Flag for absorption by sodium
  LOGICAL :: l_k                                                  = .FALSE.
!   Flag for absorption by potassium
  LOGICAL :: l_co2                                                = .FALSE.
!   Flag for absorption by carbon dioxide
  LOGICAL :: l_feh                                                = .FALSE.
!   Flag for absorption by iron hydride
  LOGICAL :: l_crh                                                = .FALSE.
!   Flag for absorption by chromium hydride
  LOGICAL :: l_li                                                 = .FALSE.
!   Flag for absorption by lithium
  LOGICAL :: l_rb                                                 = .FALSE.
!   Flag for absorption by rubidium
  LOGICAL :: l_cs                                                 = .FALSE.
!   Flag for absorption by cesium
  LOGICAL :: l_ph3                                                = .FALSE.
!   Flag for absorption by phosphine
  LOGICAL :: l_c2h2                                               = .FALSE.
!   Flag for absorption by acetylene
  LOGICAL :: l_hcn                                                = .FALSE.
!   Flag for absorption by hydrogen cyanide
  LOGICAL :: l_o3                                                 = .FALSE.
!   Flag for absorption by ozone
  LOGICAL :: l_so2                                                = .FALSE.
!   Flag for absorption by sulfur dioxide
  LOGICAL :: l_no2                                                = .FALSE.
!   Flag for absorption by nitrogen dioxide
  LOGICAL :: l_n2                                                 = .FALSE.
!   Flag for absorption by nitrogen
  LOGICAL :: l_ar                                                 = .FALSE.
!   Flag for absorption by argon
  LOGICAL :: l_o                                                  = .FALSE.
!   Flag for absorption by atomic oxygen
  LOGICAL :: l_n                                                  = .FALSE.
!   Flag for absorption by atomic nitrogen
  LOGICAL :: l_hno3                                               = .FALSE.
!   Flag for absorption by nitric acid
  LOGICAL :: l_no                                                 = .FALSE.
!   Flag for absorption by nitric oxide
  LOGICAL :: l_no3                                                = .FALSE.
!   Flag for absorption by nitrate radical
  LOGICAL :: l_n2o5                                               = .FALSE.
!   Flag for absorption by nitrogen pentoxide
  LOGICAL :: l_hono                                               = .FALSE.
!   Flag for absorption by nitrous acid
  LOGICAL :: l_ho2no2                                             = .FALSE.
!   Flag for absorption by peroxynitric acid
  LOGICAL :: l_h2o2                                               = .FALSE.
!   Flag for absorption by hydrogen peroxide
  LOGICAL :: l_c2h6                                               = .FALSE.
!   Flag for absorption by ethane
  LOGICAL :: l_ch3                                                = .FALSE.
!   Flag for absorption by methyl radical
  LOGICAL :: l_h2co                                               = .FALSE.
!   Flag for absorption by formaldehyde
  LOGICAL :: l_ho2                                                = .FALSE.
!   Flag for absorption by hydroperoxy radical
  LOGICAL :: l_hdo                                                = .FALSE.
!   Flag for absorption by semiheavy water
  LOGICAL :: l_hcl                                                = .FALSE.
!   Flag for absorption by hydrogen chloride
  LOGICAL :: l_hf                                                 = .FALSE.
!   Flag for absorption by hrdrogen peroxide
  LOGICAL :: l_cosso                                              = .FALSE.
!   Flag for absorption by cis-OSSO
  LOGICAL :: l_tosso                                              = .FALSE.
!   Flag for absorption by trans-OSSO
  LOGICAL :: l_yosos                                              = .FALSE.
!   Flag for absorption by OSO-S
  LOGICAL :: l_h2s                                                = .FALSE.
!   Flag for absorption by Hydrogen sulphide
  LOGICAL :: l_cocs                                               = .FALSE.
!   Flag for absorption by Carbonyl sulphide
  LOGICAL :: l_ch3cho                                             = .FALSE.
!   Flag for absorption by Acetaldehyde
  LOGICAL :: l_ch3ooh                                             = .FALSE.
!   Flag for absorption by Methylhydroperoxide
  LOGICAL :: l_ch3coch3                                           = .FALSE.
!   Flag for absorption by acetone
  LOGICAL :: l_ch3cocho                                           = .FALSE.
!   Flag for absorption by methylglyoxal (MGLY)
  LOGICAL :: l_chocho                                             = .FALSE.
!   Flag for absorption by glyoxal
  LOGICAL :: l_c2h5cho                                            = .FALSE.
!   Flag for absorption by propanal
  LOGICAL :: l_hoch2cho                                           = .FALSE.
!   Flag for absorption by glycolaldehyde
  LOGICAL :: l_c2h5coch3                                          = .FALSE.
!   Flag for absorption by methyl ethyl ketone (MEK)
  LOGICAL :: l_mvk                                                = .FALSE.
!   Flag for absorption by methyl vinyl ketone (MVK)
  LOGICAL :: l_macr                                               = .FALSE.
!   Flag for absorption by methacrolein (MACR)
  LOGICAL :: l_pan                                                = .FALSE.
!   Flag for absorption by peroxyacetyl nitrate (PAN)
  LOGICAL :: l_ch3ono2                                            = .FALSE.
!   Flag for absorption by methylnitrate

  LOGICAL :: l_include_gas(npd_gases)                             = .FALSE.
!   Flags to treat radiative effect of gases
  LOGICAL :: l_photol_only(npd_gases)                             = .FALSE.
!   Flags to treat gases for photolysis only, ignoring affect on flux

! Properties of clouds:
  INTEGER :: i_cloud                                              = imdi
!   Cloud scheme
  INTEGER :: i_cloud_representation                               = imdi
!   Representation of clouds
  INTEGER :: i_st_water                                           = imdi
!   Type of water droplet in stratiform clouds
  INTEGER :: i_cnv_water                                          = imdi
!   Type of water droplet in convective clouds
  INTEGER :: i_st_ice                                             = imdi
!   Type of ice crystal in stratiform clouds
  INTEGER :: i_cnv_ice                                            = imdi
!   Type of ice crystal in convective clouds
  INTEGER :: i_inhom                                              = imdi
!   Method of treating cloud water content variability
  INTEGER :: i_cloud_entrapment                                   = imdi
!   Method of treating cloud entrapment of flux between layers
  INTEGER :: i_mcica_sampling                                     = imdi
!   Method of sampling sub-grid cloud using MCICA
  INTEGER :: i_overlap                                            = imdi
!   Method of treating cloud vertical overlap
  INTEGER :: i_drop_re                                            = imdi
!   Method of setting the effective radii of cloud droplets
  LOGICAL :: l_local_cnv_partition                                = .FALSE.
!   Flag to partition convective clouds between water and
!   ice using the local temperature
  LOGICAL :: l_global_cloud_top                                   = .FALSE.
!   Flag to use a global value for the topmost cloudy layer
!   (This is used to obtained bit-reproducible results
!   across different configurations of PEs on MPP systems)
  LOGICAL :: l_avg_phase_fnc                                      = .FALSE.
!   Use a grid-box average cloud phase function for sub-grid cloud

! Angular integration (including algorithmic options):
  INTEGER :: n_channel                                            = imdi
!   Number of channels in output
  INTEGER :: i_angular_integration                                = imdi
!   Method of angular integration
  INTEGER :: i_2stream                                            = imdi
!   Two-stream scheme
  INTEGER :: i_solver                                             = imdi
!   Two-stream solver
  INTEGER :: i_solver_clear                                       = imdi
!   Clear-sky solver
  INTEGER :: n_order_gauss                                        = imdi
!   Order of Gaussian quadrature
  INTEGER :: i_truncation                                         = imdi
!   Type of truncation for spherical harmonics
  INTEGER :: i_sph_algorithm                                      = imdi
!   Algorithm used for spherical harmonic calculations
  INTEGER :: n_order_phase_solar                                  = imdi
!   Order of truncation of the solar phase function
  INTEGER :: ls_global_trunc                                      = imdi
!   Global order of truncation
  INTEGER :: ms_min                                               = imdi
!   Minimum azimuthal order
  INTEGER :: ms_max                                               = imdi
!   Maximum azimuthal order
  INTEGER :: ls_brdf_trunc                                        = imdi
!   Order of truncation of BRDFs
  INTEGER :: n_order_forward                                      = imdi
!   Order of the term used to `define' the forward scattering fraction.
  INTEGER :: i_sph_mode                                           = imdi
!   Mode of operation of spherical harmonic code
  INTEGER :: i_scatter_method                                     = imdi
!   Method of treating scattering
  INTEGER :: i_solar_src                                          = imdi
!   Index of solar source function
!   i_solar_src = 1  Original Kurucz function
!   i_solar_src = 2  AER version of Kurucz
!   i_solar_src = 3  Kurucz function used in UK model
!   i_solar_src = 4  UK reduced version
!   i_solar_src = 5  Kurucz modtran
!   i_solar_src = 6  Labs Neckel
!                         Z. Sun
  INTEGER :: i_direct_tau                                         = imdi
!   Method of scaling optical depth for direct flux calculation
  LOGICAL :: l_ir_source_quad                                     = .FALSE.
!   Flag to use a quadratic source function in the IR
  LOGICAL :: l_rescale                                            = .FALSE.
!   Flag for rescaling
  LOGICAL :: l_spherical_solar                                    = .FALSE.
!   Flag to use spherical geometry for the direct beam path
  LOGICAL :: l_grey_single                                        = .FALSE.
!   Prevent unnecessary calculation of equivalent extinction with one term
  LOGICAL :: l_henyey_greenstein_pf                               = .FALSE.
!   Flag to use Henyey-Greenstein phase functions
  LOGICAL :: l_lanczos                                            = .FALSE.
!   Flag to use Lanczos smoothing of solar phf
  LOGICAL :: l_euler_trnf                                         = .FALSE.
!   Flag to apply Euler's transformation to alternating series
  REAL (RealK) :: accuracy_adaptive                               = rmdi
!   Accuracy for adaptive truncation
  REAL (RealK) :: euler_factor                                    = rmdi
!   Factor applied to the last term of an alternating series
  REAL (RealK) :: half_angle                                      = rmdi
!   Half angle of FOV used for CSR fraction calculation

! Miscallaneous options
  LOGICAL :: l_tile                                               = .FALSE.
!   Allow tiling of the surface
  LOGICAL :: l_tile_emissivity                                    = .FALSE.
!   Use tile emissivities to calculate the ground source function
  LOGICAL :: l_flux_ground                                        = .FALSE.
!   Use emission from the surface
  LOGICAL, ALLOCATABLE :: l_flux_tile(:)
!   Use emission from selected tiles
  LOGICAL :: l_extra_top                                          = .FALSE.
!   Flag to insert an extra layer into radiation above the
!   top of the model (this is sometimes desirable to ensure
!   proper radiative heating in the top resolved layer).
  LOGICAL :: l_rad_deg                                            = .FALSE.
!   Flag to apply spatial degradation in radiation: in this case
!   radiation quantities are interpolated to all grid-points,
!   whereas subsampling refers to selecting a portion of the area
!   on the PE and returning radiative quentatities only at those
!   points
  LOGICAL :: l_subsample                                          = .FALSE.
!   Flag to apply spatial subsampling (for satellite footprints)
!   in radiation.
  LOGICAL :: l_mixing_ratio                                       = .FALSE.
!   True if mixing ratios are with respect to dry mass
  LOGICAL :: l_map_sub_bands                                      = .FALSE.
!   Map sub-bands to channels rather than bands


! Band-by-band control options
  INTEGER, ALLOCATABLE :: i_scatter_method_band(:)
!   Method of treating scattering in each band
  INTEGER, ALLOCATABLE :: i_gas_overlap_band(:)
!   Gas overlap assumption in each band
  INTEGER, ALLOCATABLE :: map_channel(:)
!   Mapping of actual bands or sub-bands to the output channels
  REAL (RealK), ALLOCATABLE :: weight_band(:)
!   Weighting function for bands
  REAL (RealK), ALLOCATABLE :: weight_diag(:)
!   Weighting function for bands for diagnostic fluxes
  LOGICAL, ALLOCATABLE :: l_clear_band(:)
!   Calculate clear-sky fluxes for band


! Switches for diagnostic output
  LOGICAL :: l_clear                                              = .FALSE.
!   Calculate clear-sky fluxes
  LOGICAL :: l_clean                                              = .FALSE.
!   Calculate clean-air fluxes
  LOGICAL :: l_clear_clean                                        = .FALSE.
!   Calculate clear-clean fluxes
  LOGICAL :: l_flux_div                                           = .FALSE.
!   Calculate flux divergence
  LOGICAL :: l_blue_flux_surf                                     = .FALSE.
!   Calculate blue surface fluxes
  LOGICAL :: l_cloud_absorptivity                                 = .FALSE.
!   Calculate absorptivity of clouds (only infra-red)
  LOGICAL :: l_cloud_extinction                                   = .FALSE.
!   Calculate extinction of clouds (only solar)
  LOGICAL :: l_ls_cloud_absorptivity                              = .FALSE.
!   Calculate absorptivity of layer clouds (only infra-red)
  LOGICAL :: l_ls_cloud_extinction                                = .FALSE.
!   Calculate extinction of layer clouds (only solar)
  LOGICAL :: l_cnv_cloud_absorptivity                             = .FALSE.
!   Calculate absorptivity of conv.clouds (only infra-red)
  LOGICAL :: l_cnv_cloud_extinction                               = .FALSE.
!   Calculate extinction of conv.clouds (only solar)
  LOGICAL :: l_flux_direct_band                                   = .FALSE.
!   Calculate direct flux per band
  LOGICAL :: l_flux_direct_div_band                               = .FALSE.
!   Calculate direct flux divergence per band
  LOGICAL :: l_flux_direct_sph_band                               = .FALSE.
!   Calculate direct flux for spherical path per band
  LOGICAL :: l_flux_down_band                                     = .FALSE.
!   Calculate total downward flux per band
  LOGICAL :: l_flux_up_band                                       = .FALSE.
!   Calculate upward flux per band
  LOGICAL :: l_flux_div_band                                      = .FALSE.
!   Calculate flux divergence per band
  LOGICAL :: l_flux_direct_clear_band                             = .FALSE.
!   Calculate clear-sky direct flux per band
  LOGICAL :: l_flux_direct_clear_div_band                         = .FALSE.
!   Calculate clear-sky direct flux divergence per band
  LOGICAL :: l_flux_direct_clear_sph_band                         = .FALSE.
!   Calculate clear-sky direct flux for spherical path per band
  LOGICAL :: l_flux_down_clear_band                               = .FALSE.
!   Calculate clear-sky downward flux per band
  LOGICAL :: l_flux_up_clear_band                                 = .FALSE.
!   Calculate clear-sky upward flux per band
  LOGICAL :: l_flux_div_clear_band                                = .FALSE.
!   Calculate clear-sky flux divergence per band
  LOGICAL :: l_actinic_flux                                       = .FALSE.
!   Calculate actinic flux per channel
  LOGICAL :: l_actinic_flux_clear                                 = .FALSE.
!   Calculate clear-sky actinic flux per channel
  LOGICAL :: l_actinic_flux_band                                  = .FALSE.
!   Calculate actinic flux per band
  LOGICAL :: l_actinic_flux_clear_band                            = .FALSE.
!   Calculate clear-sky actinic flux per band
  LOGICAL :: l_photolysis_rate                                    = .FALSE.
!   Calculate photolysis rates per channel
  LOGICAL :: l_photolysis_rate_clear                              = .FALSE.
!   Calculate clear-sky photolysis rates per channel
  LOGICAL :: l_photolysis_rate_clean                              = .FALSE.
!   Calculate clean-air photolysis rates per channel
  LOGICAL :: l_photolysis_rate_clear_clean                        = .FALSE.
!   Calculate clear-clean photolysis rates per channel
  LOGICAL :: l_photolysis_div                                     = .FALSE.
!   Calculate flux divergence for photolysis
  LOGICAL :: l_photolysis_div_clear                               = .FALSE.
!   Calculate clear-sky flux divergence for photolysis
  LOGICAL :: l_aerosol_absorption_band                            = .FALSE.
!   Calculate total aerosol absorption per band
  LOGICAL :: l_aerosol_scattering_band                            = .FALSE.
!   Calculate total aerosol scattering per band
  LOGICAL :: l_aerosol_asymmetry_band                             = .FALSE.
!   Calculate total aerosol asymmetry (weighted by scattering) per band
  LOGICAL :: l_spherical_path_diag                                = .FALSE.
!   Output the direct beam path through spherical layers
  LOGICAL :: l_contrib_func                                       = .FALSE.
!   Output the contribution function
  LOGICAL :: l_contrib_func_band                                  = .FALSE.
!   Output the contribution function per band

! Satellite Data:

! Current provisions are largely intended for geostationary satellites.
! Note that in accordance with the UM's conventions these must be
! given in SI units; i.e. angles will be in radians etc.

  LOGICAL :: l_geostationary                                      = .FALSE.
!   Flag to signal that a geostationary satellite is assumed.
  CHARACTER (LEN=80) :: sat_desc                                  = ''
!   String for description of satellite
  REAL (RealK) :: sat_hgt                                         = rmdi
!   Height of the orbit above the Earth's surface
  REAL (RealK) :: sat_lon                                         = rmdi
!   Longitude of the (geostationary) satellite
  REAL (RealK) :: sat_lat                                         = rmdi
!   Latitude of the (geostationary) satellite (in practice, for
!   a geostationary satellite this must be 0.0)

! Viewing domain:
  REAL (RealK) :: max_view_lon                                    = rmdi
!   Maximum longitude of viewing domain
  REAL (RealK) :: min_view_lon                                    = rmdi
!   Minimum longitude of viewing domain
  REAL (RealK) :: max_view_lat                                    = rmdi
!   Maximum latitude of viewing domain
  REAL (RealK) :: min_view_lat                                    = rmdi
!   Minimum latitude of viewing domain

END TYPE StrCtrl


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_control(control, sp, n_tile)

USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrCtrl),     INTENT(INOUT) :: control
TYPE (StrSpecData), INTENT(IN), OPTIONAL :: sp
INTEGER,            INTENT(IN), OPTIONAL :: n_tile

IF (PRESENT(sp)) THEN

  IF (.NOT. ALLOCATED(control%i_scatter_method_band))       &
    ALLOCATE(control%i_scatter_method_band ( sp%dim%nd_band ))

  IF (.NOT. ALLOCATED(control%i_gas_overlap_band))          &
    ALLOCATE(control%i_gas_overlap_band    ( sp%dim%nd_band ))

  IF (.NOT. ALLOCATED(control%map_channel))                 &
    ALLOCATE(control%map_channel ( MAX(sp%dim%nd_sub_band, sp%dim%nd_band) ))

  IF (.NOT. ALLOCATED(control%weight_band))                 &
    ALLOCATE(control%weight_band           ( sp%dim%nd_band ))

  IF (.NOT. ALLOCATED(control%weight_diag))                 &
    ALLOCATE(control%weight_diag           ( sp%dim%nd_band ))

  IF (.NOT. ALLOCATED(control%l_clear_band)) THEN
    ALLOCATE(control%l_clear_band          ( sp%dim%nd_band ))
    control%l_clear_band(1:sp%dim%nd_band) = .FALSE.
  END IF

END IF

IF (PRESENT(n_tile)) THEN

  IF (.NOT. ALLOCATED(control%l_flux_tile)) THEN
    ALLOCATE(control%l_flux_tile ( n_tile ))
    control%l_flux_tile = .FALSE.
  END IF

END IF

END SUBROUTINE allocate_control
!------------------------------------------------------------------------------
SUBROUTINE deallocate_control(control)

IMPLICIT NONE

TYPE (StrCtrl), INTENT(INOUT) :: control

IF (ALLOCATED(control%l_flux_tile))  DEALLOCATE(control%l_flux_tile)
IF (ALLOCATED(control%l_clear_band)) DEALLOCATE(control%l_clear_band)
IF (ALLOCATED(control%weight_diag))  DEALLOCATE(control%weight_diag)
IF (ALLOCATED(control%weight_band))  DEALLOCATE(control%weight_band)
IF (ALLOCATED(control%map_channel))  DEALLOCATE(control%map_channel)
IF (ALLOCATED(control%i_gas_overlap_band))    &
                                     DEALLOCATE(control%i_gas_overlap_band)
IF (ALLOCATED(control%i_scatter_method_band)) &
                                     DEALLOCATE(control%i_scatter_method_band)

END SUBROUTINE deallocate_control
!------------------------------------------------------------------------------

END MODULE def_control

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to define headers and suffixes for field input

MODULE input_head_pcf

! Description:
!
!   This module defines suffixes for files, identifying the contents
!   and headers for columns of input data for generating files.
!   It also defines information on physical units.
!
!   IMPORTANT: Initialization by array constructors is used here.
!   Great care must be taken if physical identifiers (IP_...) are
!   changed to ensure that the changes are mirrored in all arrays.

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: npd_aerosol_component

  IMPLICIT NONE


  INTEGER, Parameter :: NP_max_length_line   = 132
!   Maximum length of input line
  INTEGER, Parameter :: NPD_in_profile       = 9
!   Maximum number of input profiles
  INTEGER, Parameter :: NPD_data_column      = 35
!   Maximum number of columns of data
  INTEGER, Parameter :: NPD_phys_type        = 62
!   Maximum number of types of data
  INTEGER, Parameter :: NPD_unit             = 20
!   Number of physical units
!
!
! Lengths of character strings
  INTEGER, Parameter :: len_col_header       = 10
!   Length of character strings for column headers
  INTEGER, Parameter :: len_long_title       = 30
!   Length of character strings for file titles 
  INTEGER, Parameter :: len_file_suffix      = 12
!   Length of character strings for file suffixes
!
!
!
  INTEGER, Parameter :: IP_physical_data = 1
!   Index for physical data group
  INTEGER, Parameter :: IP_gaseous_data  = 2
!   Index for gaseous data group
  INTEGER, Parameter :: IP_aerosol_data  = 3
!   Index for aerosol data group
!
!
! Physical Data:
!
  INTEGER, Parameter :: IP_pressure                              = 1
!   Pressure
  INTEGER, Parameter :: IP_height                                = 2
!   Height
  INTEGER, Parameter :: IP_temperature                           = 3
!   Temperature
  INTEGER, Parameter :: IP_ozone_density                         = 4
!   Density of ozone
  INTEGER, Parameter :: IP_water_vapour_density                  = 5
!   Density of water vapour
  INTEGER, Parameter :: IP_lwc                                   = 6
!   Liquid water content
  INTEGER, Parameter :: IP_re                                    = 7
!   (Stratiform) effective radius/characteristic dimension
!   of water cloud
  INTEGER, Parameter :: IP_null_field                            = 8
!   Null field
  INTEGER, Parameter :: IP_temperature_ground                    = 9
!   Surface temperature
  INTEGER, Parameter :: IP_cloud_fraction                        = 10
!   Cloud fraction
  INTEGER, Parameter :: IP_flux_up                               = 11
!   Upward flux
  INTEGER, Parameter :: IP_flux_down                             = 12
!   Downward flux
  INTEGER, Parameter :: IP_flux_direct                           = 13
!   Direct flux
  INTEGER, Parameter :: IP_flux_net                              = 14
!   Net downward flux
  INTEGER, Parameter :: IP_heating_rate                          = 15
!   Heating rate
  INTEGER, Parameter :: IP_flux_total_down                       = 16
!   Total downward flux
  INTEGER, Parameter :: IP_lwm                                   = 17
!   Liquid water mass mixing ratio
  INTEGER, Parameter :: IP_t_dew                                 = 18
!   Dew point temperature
  INTEGER, Parameter :: IP_spec_humidity                         = 19
!   Specific humidity
  INTEGER, Parameter :: IP_hum_mix_ratio                         = 20
!   Humidity mixing ratio
  INTEGER, Parameter :: IP_iwc                                   = 21
!   (Stratiform) ice water content
  INTEGER, Parameter :: IP_iwm                                   = 22
!   (Stratiform) mass mixing ratio of ice
  INTEGER, Parameter :: IP_ire                                   = 23
!   (Stratiform) ice effective radius/characteristic dimension
  INTEGER, Parameter :: IP_tau                                   = 24
!   Optical depth
  INTEGER, Parameter :: IP_ssa                                   = 25
!   Single scattering albedo
  INTEGER, Parameter :: IP_gsc                                   = 26
!   Asymmetry
  INTEGER, Parameter :: IP_fsc                                   = 27
!   Forward scattering
  INTEGER, Parameter :: IP_planck                                = 28
!   Planck function
  INTEGER, Parameter :: IP_cloud_fraction_cv                     = 29
!   Convective cloud fraction
  INTEGER, Parameter :: IP_lwccv                                 = 30
!   Convective LWC
  INTEGER, Parameter :: IP_lwmcv                                 = 31
!   Mass mixing ratio of convective liquid water
  INTEGER, Parameter :: IP_recv                                  = 32
!   (Convective) effective radius/characteristic dimension
!   of water cloud
  INTEGER, Parameter :: IP_iwccv                                 = 33
!   IWC of convective cloud
  INTEGER, Parameter :: IP_iwmcv                                 = 34
!   Mass mixing ratio of ice in convective cloud
  INTEGER, Parameter :: IP_irecv                                 = 35
!   Convective ice effective radius/characteristic dimension
  INTEGER, Parameter :: IP_pressure_ground                       = 36
!   Surface pressure
  INTEGER, Parameter :: IP_t_level                               = 37
!   Temperatures on levels
  INTEGER, Parameter :: IP_rel_hum                               = 38
!   relative humidity
  INTEGER, Parameter :: IP_discard                               = 39
!   Discard data
  INTEGER, Parameter :: IP_cloud_fraction_w                      = 40
!   Water cloud fraction
  INTEGER, Parameter :: IP_cloud_fraction_i                      = 41
!   Ice cloud fraction
  INTEGER, Parameter :: IP_cloud_fraction_w_cv                   = 42
!   Convective water cloud fraction
  INTEGER, Parameter :: IP_cloud_fraction_i_cv                   = 43
!   Convective ice cloud fraction
  INTEGER, Parameter :: IP_solar_zenith_angle                    = 44
!   Solar zenith angle
  INTEGER, Parameter :: IP_solar_azimuth                         = 45
!   Solar azimuthal angle
  INTEGER, Parameter :: IP_solar_toa                             = 46
!   Solar toa irradiance
  INTEGER, Parameter :: IP_surface_type                          = 47
!   Type of surface
  INTEGER, Parameter :: IP_view_polar                            = 48
!   Polar viewing angle
  INTEGER, Parameter :: IP_view_azim                             = 49
!   Azimuthal viewing angle
  INTEGER, Parameter :: IP_radiance                              = 50
!   Radiances
  INTEGER, Parameter :: IP_surface_char                          = 51
!   Surface chracteristics
  INTEGER, Parameter :: IP_optical_water                         = 52
!   Optical properties of water clouds
  INTEGER, Parameter :: IP_optical_ice                           = 53
!   Optical properties of ice clouds
  INTEGER, Parameter :: IP_optical_ss                            = 54
!   Optical properties of an atmospheric column
  INTEGER, Parameter :: IP_iso_inc                               = 55
!   Isotropic radiance at the top of the atmosphere
  INTEGER, Parameter :: IP_view_geom                             = 56
!   Viewing geometry
  INTEGER, Parameter :: IP_photolysis                            = 57
!   Rate of photolysis
  INTEGER, Parameter :: IP_p_level                               = 58
!   Pressure on levels (layer boundaries)
  INTEGER, Parameter :: IP_contrib_funci                         = 59
!   Contribution function (intensity)
  INTEGER, Parameter :: IP_contrib_funcf                         = 60
!   Contribution function (flux)
  INTEGER, Parameter :: IP_actinic_flux                          = 61
!   Actinic flux
  INTEGER, Parameter :: IP_photolysis_rate                       = 62
!   Photolysis rate

!
  CHARACTER  (LEN=len_col_header), Parameter, &
    Dimension(NPD_phys_type) :: header_phys = (/ &
    'PRESS     ', 'HGT       ', 'TEMP      ', 'O3DEN     ', &
    'H2ODEN    ', 'LWC       ', 'RE        ', '          ', &
    'TSTAR     ', 'CLFRAC    ', 'UFLX      ', 'DFLX      ', &
    'SFLX      ', 'NFLX      ', 'HRTS      ', 'VFLX      ', &
    'LWM       ', 'TDEW      ', 'SPH       ', 'HMR       ', &
    'IWC       ', 'IWM       ', 'IRE       ', 'TAU       ', &
    'SSA       ', 'ASYM      ', 'FRWSC     ', 'PLANCK    ', &
    'CVFRAC    ', 'LWCCV     ', 'LWMCV     ', 'RECV      ', &
    'IWCCV     ', 'IWMCV     ', 'IRECV     ', 'PSTAR     ', &
    'TLEV      ', 'RH        ', 'DISCRD    ', 'WCLFRC    ', &
    'ICLFRC    ', 'CWFRAC    ', 'CIFRAC    ', 'SZEN      ', &
    'SAZIM     ', 'STOA      ', 'SURF      ', 'POLAR     ', &
    'AZIM      ', 'RADN      ', 'SRFCHR    ', 'OPWT      ', &
    'OPICE     ', 'OPSS      ', 'ISOS      ', 'GEOM      ', &
    'PHOTOL    ', 'PLEV      ', 'CFI       ', 'CFF       ', &
    'AFLX      ', 'PHRATE    ' /) 
!   Headers for physical data
!
  CHARACTER  (LEN=len_file_suffix), Parameter, &
    Dimension(NPD_phys_type) :: phys_suffix = (/ &
    'p           ', 'hgt         ', 't           ', 'o3d         ', &
    'qd          ', 'lwc         ', 're          ', 'null        ', &
    'tstar       ', 'clfr        ', 'uflx        ', 'dflx        ', &
    'sflx        ', 'nflx        ', 'hrts        ', 'vflx        ', &
    'lwm         ', 'tdw         ', 'q           ', 'hmr         ', &
    'iwc         ', 'iwm         ', 'ire         ', 'tau         ', &
    'ssa         ', 'gsc         ', 'fsc         ', 'plk         ', &
    'ccfr        ', 'lwccv       ', 'lwmcv       ', 'recv        ', &
    'iwccv       ', 'iwmcv       ', 'irecv       ', 'pstar       ', &
    'tl          ', 'rh          ', '            ', 'wclfr       ', &
    'iclfr       ', 'wccfr       ', 'iccfr       ', 'szen        ', &
    'sazim       ', 'stoa        ', 'surf        ', 'vwpol       ', &
    'vwazim      ', 'radn        ', 'surf        ', 'op_water    ', &
    'op_ice      ', 'ss          ', 'isos        ', 'view        ', &
    'photol      ', 'pl          ', 'cfi         ', 'cff         ', &
    'aflx        ', 'ph_rate     '  /) 
!   File suffixes for physical data
!
!
  CHARACTER  (LEN=len_long_title), Parameter, &
    Dimension(NPD_phys_type) :: phys_title = (/ &
    ' Pressure                     ', ' Height                       ', &
    ' Temperature                  ', ' Ozone Mass Density           ', &
    ' Water Vapour Mass Density    ', ' Liquid Water Content         ', &
    ' Effective Radius             ', ' Null Grid                    ', &
    ' Temperature of surface       ', ' Cloud Fraction               ', &
    ' Upward Flux                  ', ' Downward Flux (diffuse)      ', &
    ' Direct Flux                  ', ' Net Downward Flux            ', &
    ' Heating Rates (K, day)       ', ' Total Downward Flux          ', &
    ' Liquid Water Mass Fraction   ', ' Dew point Temperature        ', &
    ' Specific Humidity            ', ' Humidity Mixing Ratio        ', &
    ' Ice Content                  ', ' Ice Mass Fraction            ', &
    ' Ice Effective Radius         ', ' Optical Depth                ', &
    ' Albedo of Single Scattering  ', ' Asymmetry Factor             ', &
    ' Forward Scattering Factor    ', ' Planck Function              ', &
    ' Convective Cloud Fraction    ', ' Conv. Liquid Water Content   ', &
    ' Conv. Liquid Water Mass Frac.', ' Conv. Effective Radius       ', &
    ' Convective Ice Content       ', ' Convective Ice Mass Fraction ', &
    ' Conv. Ice Effective Radius   ', ' Pressure at surface          ', &
    ' Temperature on Levels        ', ' Relative Humidity            ', &
    '                              ', ' Water Cloud Fraction         ', &
    ' Ice Cloud Fraction           ', ' Conv. Water Cloud Fraction   ', &
    ' Conv. Ice Cloud Fraction     ', ' Solar zenith angle           ', &
    ' Solar azimuthal angle        ', ' Solar Irradiance at TOA      ', &
    ' Type of surface              ', ' Polar viewing angle          ', &
    ' Azimuthal viewing angle      ', ' Radiance                     ', &
    ' Surface characteristics      ', ' Optical data for droplets    ', &
    ' Optical data for ice crystals', ' Single scattering properties ', &
    ' Isotropic source             ', ' Viewing Geometry             ', &
    ' Rate of photolysis           ', ' Pressure on Levels           ', &
    ' Contribution function (inty) ', ' Contribution function (flux) ', &
    ' Actinic flux                 ', ' Photolysis rate              '/)
!   Long titles for physical data
!
!
! Aerosol Data:
!
  CHARACTER  (LEN=len_col_header), Parameter, &
    Dimension(NPD_aerosol_component) :: header_aerosol = (/ &
    'WTSOL     ', 'DUST      ', 'OCN       ', 'SOOT      ', &
    'ASH       ', 'SULPH     ', 'NH4SO4    ', 'AUNCH     ', &
    'SAHARA    ', 'ACCUM     ', 'AITKEN    ', 'FRSOOT    ', &
    'AGSOOT    ', 'NACL      ', 'NACLFLM   ', 'NACLJET   ', &
    'DUSTDIV1  ', 'DUSTDIV2  ', 'DUSTDIV3  ', 'DUSTDIV4  ', &
    'DUSTDIV5  ', 'DUSTDIV6  ', 'BIOMS1    ', 'BIOMS2    ', &
    'BIOGENIC  ', 'FROCFF    ', 'AGOCFF    ', 'DELTA     ', &
    '          ', 'NITRATE   ', 'DUST2BIN1 ', 'DUST2BIN2 ' /) 
!   Headers for aerosol data
!
!
  CHARACTER  (LEN=len_file_suffix), Parameter, &
    Dimension(NPD_aerosol_component) :: aerosol_suffix = (/ &
    'wtsol       ', 'dust        ', 'ocn         ', 'soot        ', &
    'ash         ', 'sulph       ', 'nh4so4      ', 'aunch       ', &
    'sahara      ', 'accum       ', 'aitken      ', 'frsoot      ', &
    'agsoot      ', 'nacl        ', 'naclflm     ', 'nacljet     ', &
    'dustdiv1    ', 'dustdiv2    ', 'dustdiv3    ', 'dustdiv4    ', &
    'dustdiv5    ', 'dustdiv6    ', 'bioms1      ', 'bioms2      ', &
    'biogenic    ', 'frocff      ', 'agocff      ', 'delta       ', &
    '            ', 'nitrate     ', 'dust2bin1   ', 'dust2bin2   ' /) 
!   File suffixes for aerosol data
!
!
  CHARACTER  (LEN=len_file_suffix), Parameter, &
    Dimension(NPD_aerosol_component) :: aerosol_opt_suffix = (/ &
    'op_wtsol    ', 'op_dust     ', 'op_ocn      ', 'op_soot     ', &
    'op_ash      ', 'op_sulph    ', 'op_nh4so4   ', 'op_aunch    ', &
    'op_sahara   ', 'op_accum    ', 'op_aitken   ', 'op_frsoot   ', &
    'op_agsoot   ', 'op_nacl     ', 'op_naclflm  ', 'op_nacljet  ', &
    'op_dustdiv1 ', 'op_dustdiv2 ', 'op_dustdiv3 ', 'op_dustdiv4 ', &
    'op_dustdiv5 ', 'op_dustdiv6 ', 'op_bioms1   ', 'op_bioms2   ', &
    'op_biogenic ', 'op_frocff   ', 'op_agocff   ', 'op_delta    ', &
    '            ', 'op_nitrate  ', 'op_dust2bin1', 'op_dust2bin2' /)
!   File suffixes for optical properties of aerosols
!
!
  CHARACTER  (LEN=len_long_title), Parameter, &
    Dimension(NPD_aerosol_component) :: aerosol_title = (/ &
    ' Water-soluble Aerosol        ', ' Dust-like Aerosol            ', &
    ' Oceanic Aerosol              ', ' Soot Aerosol                 ', &
    ' Ash Aerosol                  ', ' 75% Sulphuric acid Aerosol   ', &
    ' Ammonium Sulphate Aerosol    ', ' Uncharacterized Aerosol      ', &
    ' Saharan Dust Aerosol         ', ' Accumulation SO4 Aerosol     ', &
    ' Aitken-mode SO4 Aerosol      ', ' Fresh Soot Aerosol           ', &
    ' Aged Soot Aerosol            ', ' Generic Sodium Chloride      ', &
    ' Sodium Chloride (Film mode)  ', ' Sodium Chloride (Jet mode)   ', &
    ' Dust aerosol (division 1)    ', ' Dust aerosol (division 2)    ', &
    ' Dust aerosol (division 3)    ', ' Dust aerosol (division 4)    ', &
    ' Dust aerosol (division 5)    ', ' Dust aerosol (division 6)    ', &
    ' Biomass aerosol (division 1) ', ' Biomass aerosol (division 2) ', &
    ' Biogenic aerosol             ', ' Fresh fossil-fuel org. carbon', &
    ' Aged fossil-fuel org. carbon ', ' Unspecified (delta) aerosol  ', &
    '                              ', ' Ammonium nitrate aerosol     ', &
    ' Two-bin Dust aerosol (div 1) ', ' Two-bin Dust aerosol (div 2) ' /)
!   Long titles for aerosol data
!
!
!
! Units for data:
!
  INTEGER, Parameter :: IP_unit_pa                = 1
!   Units of pascals
  INTEGER, Parameter :: IP_unit_mb                = 2
!   Units of millibars
  INTEGER, Parameter :: IP_unit_k                 = 3
!   Units of kelvins
  INTEGER, Parameter :: IP_unit_m                 = 4
!   Units of metres
  INTEGER, Parameter :: IP_unit_km                = 5
!   Units of kilometres
  INTEGER, Parameter :: IP_unit_um                = 6
!   Units of microns
  INTEGER, Parameter :: IP_unit_kgm_3             = 7
!   Units of kilograms per cubic metre
  INTEGER, Parameter :: IP_unit_gm_3              = 8
!   Units of grams per cubic metre
  INTEGER, Parameter :: IP_unit_kgm_2             = 9
!   Units of kilograms per square metre
  INTEGER, Parameter :: IP_unit_gm_2              = 10
!   Units of grams per square metre
  INTEGER, Parameter :: IP_unit_gkg               = 11
!   Units of grams per kilogram
  INTEGER, Parameter :: IP_unit_g_g               = 12
!   Units of grams per gram
  INTEGER, Parameter :: IP_unit_none              = 13
!   Dummy unit for dimensionless input
  INTEGER, Parameter :: IP_unit_m_3               = 14
!   Unit number per cubic metre
  INTEGER, Parameter :: IP_unit_cm_3              = 15
!   Unit number per cubic centimetre
  INTEGER, Parameter :: IP_unit_wm_2              = 16
!   Units of watts per cubic metre
  INTEGER, Parameter :: IP_unit_c                 = 17
!   Units of degrees celsius
!
  CHARACTER  (LEN=len_col_header), Parameter, &
    Dimension(NPD_unit) :: name_unit = (/ &
    'PA        ', 'MB        ', 'K         ', 'M         ', &
    'KM        ', 'UM        ', 'KGM-3     ', 'GM-3      ', &
    'KGM-2     ', 'GM-2      ', 'GKG-1     ', 'GG-1      ', &
    'NONE      ', 'M-3       ', 'CM-3      ', 'WM-2      ', &
    'C         ', 'PPMV      ', 'VOL-MARS  ', 'VOL-VENUS ' /) 
!   Names of units
!
!
  REAL  (RealK), Parameter, Dimension(NPD_unit) :: factor_unit = (/ &
    1.0_RealK,     1.0E+02_RealK, 1.0_RealK,     1.0_RealK, &
    1.0E+03_RealK, 1.0E-06_RealK, 1.0_RealK,     1.0E-03_RealK, &
    1.0_RealK,     1.0E-03_RealK, 1.0E-03_RealK, 1.0_RealK, &
    1.0_RealK,     1.0_RealK,     1.0E+06_RealK, 1.0_RealK, &
    1.0_RealK,     -1.0E-06_RealK/28.966_RealK, &
    -1.0_RealK/43.34_RealK, -1.0_RealK/43.45_RealK /)
!   Factors for conversion to S.I. units
  REAL  (RealK), Parameter, Dimension(NPD_unit) :: offset_unit = (/ &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, &
    0.0_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK, &
    2.7315E+02_RealK, 0.0_RealK, 0.0_RealK, 0.0_RealK /)
!   Offsets for conversion to S.I. units
!
!
!
END MODULE input_head_pcf

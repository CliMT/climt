! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Driver using the runes interface to run the two-stream code
! taking input fields from a netCDF file
!
! Method:
!   A netCDF file as it can be produced during an LFRic run with 
!   radiation switched on together with the appropriate configuration 
!   namelist file drive calls to runes to produce desired fluxes.
!   These are added to the netCDF file which may contain the fluxes
!   generated in the LFRic run for possible comparison.
!   
! ---------------------------------------------------------------------

PROGRAM runes_nc

  ! Modules

  USE socrates_set_spectrum, only: set_spectrum, &
                                   get_spectrum, &
                                   set_mcica

  USE socrates_runes, only: runes, &
                            StrDiag, &
                            ip_source_illuminate, &
                            ip_source_thermal, &
                            ip_scatter_full, &
                            ip_scatter_none, &
                            ip_scatter_approx, &
                            ip_scatter_hybrid, &
                            ip_cloud_representation_off, &
                            ip_cloud_representation_ice_water, &
                            ip_cloud_representation_combine_ice_water, &
                            ip_cloud_representation_csiw, &
                            ip_cloud_representation_split_ice_water, &
                            ip_overlap_max_random, &
                            ip_overlap_random, &
                            ip_overlap_exponential_random, &
                            ip_inhom_homogeneous, &
                            ip_inhom_scaling, &
                            ip_inhom_mcica, &
                            ip_inhom_cairns, &
                            ip_inhom_tripleclouds_2019, &
                            ip_cloud_entrapment_max, &
                            ip_cloud_entrapment_zero, &
                            ip_droplet_re_constant, &
                            ip_droplet_re_liu, &
                            ip_droplet_re_default

  USE realtype_rd, only : RealExt

  USE rad_pcf, only: i_normal, &
                     i_err_fatal

  USE errormessagelength_mod, only: errormessagelength

  USE ereport_mod, only: ereport

  USE netcdf

  ! Defines the relevant sections of the configuration namelist
  ! (aerosol, radiation, radiative gases)
  USE runes_nc_nml_mod

  IMPLICIT NONE

  ! Interface blocks for the read and write netcdf subroutines

  INTERFACE

    SUBROUTINE read_nc(input_file,var_name,i_time,dataout)
    USE realtype_rd, ONLY : RealExt
    CHARACTER(LEN=*), INTENT(in) :: input_file
       ! Name of the input file (netCDF)
    CHARACTER(LEN=*), INTENT(in) :: var_name
       ! Name of the field to be read in
    INTEGER, INTENT(in) :: i_time
       ! Time step
    REAL(RealExt), ALLOCATABLE, INTENT(out) :: dataout(:)
       ! On exit holds field as a 1D array
    INTEGER :: ierr
       ! Error flag
    INTEGER :: ncid
       ! netCDF file ID
    INTEGER :: varid
       ! netCDF variable ID
    INTEGER :: ndims
       ! netCDF number of dimensions
    INTEGER, ALLOCATABLE :: dimids(:)
       ! netCDF dimension IDs
    CHARACTER(LEN=200), ALLOCATABLE :: dimension_name(:)
       ! netCDF dimension name
    INTEGER, ALLOCATABLE :: dimension_size(:)
       ! netCDF dimension size
    INTEGER :: vartype
       ! netCDF variable type
    INTEGER :: i,j
       ! Loop indices
    END SUBROUTINE read_nc

    SUBROUTINE write_nc(output_file,var_name,i_time,dataout)
    USE realtype_rd, ONLY : RealExt
    CHARACTER(LEN=*), INTENT(in) :: output_file
       ! Name of the output file (netCDF)
    CHARACTER(LEN=*), INTENT(in) :: var_name
       ! Name of the field to be written out
    INTEGER, INTENT(in) :: i_time
       ! Time step
    REAL(RealExt), ALLOCATABLE, INTENT(in) :: dataout(:)
       ! On entry holds field as a 1D array
    INTEGER :: ierr
       ! Error flag
    INTEGER :: ncid
       ! netCDF file ID
    INTEGER :: varid
       ! netCDF variable ID
    INTEGER :: ndims
       ! netCDF number of dimensions
    INTEGER, ALLOCATABLE :: dimids(:)
       ! netCDF dimension IDs
    CHARACTER(LEN=200), ALLOCATABLE :: dimension_name(:)
       ! netCDF dimension name
    INTEGER, ALLOCATABLE :: dimension_size(:)
       ! netCDF dimension size
    INTEGER :: vartype
       ! netCDF variable type
    INTEGER :: i,j
       ! Loop indices
    END SUBROUTINE write_nc

  END INTERFACE

  REAL(RealExt), PARAMETER :: seconds_per_day    = 8.6400E+04_RealExt
     ! Conversion factor to obtain heating rates per day
  REAL(RealExt), PARAMETER :: cloud_fraction_min = 0.001_RealExt
     ! Threshold fraction for clear/cloudy layer

  ! Dimensions of fields in the netCDF file

  INTEGER :: n_profile
     ! Number of columns
  INTEGER :: nlayers
     ! Number of layers
  INTEGER :: n_full_levels
     ! Number of full levels
  INTEGER :: n_time
     ! Number of time steps
  INTEGER :: n_surf_tile
     ! Number of surface tiles
  INTEGER :: n_lw_band
     ! Number of LW bands
  INTEGER :: n_sw_band
     ! Number of SW bands
  INTEGER :: n_lw_band_tile
     ! Number of LW bands times tiles
  INTEGER :: n_sw_band_tile
     ! Number of SW bands times tiles
  INTEGER :: n_lw_band_aero
     ! Number of LW bands times MODE aerosol
  INTEGER :: n_sw_band_aero
     ! Number of SW bands times MODE aerosol
  INTEGER :: n_aero
     ! Number of MODE aerosols
  INTEGER :: n_aer_mode_sw
     ! Number of SW MODE aerosols
  INTEGER :: n_aer_mode_lw
     ! Number of LW MODE aerosols

  ! Cloud settings
  INTEGER :: i_cloud_representation
     ! Cloud representation
  INTEGER :: i_overlap
     ! Cloud overlap treatment
  INTEGER :: i_inhom
     ! Cloud inhomogenity treatment
  INTEGER :: i_drop_re
     ! Cloud particle effective radius treatment
  INTEGER :: i_scatter_method_lw
     ! Scattering treatment
  INTEGER :: i_cloud_entrapment
     ! Treatment of entrapment

  ! Switches for radiatively active gases
  LOGICAL :: l_cfc11
     ! Include/exclude CFC11
  LOGICAL :: l_cfc113
     ! Include/exclude CFC113
  LOGICAL :: l_cfc12
     ! Include/exclude CFC12
  LOGICAL :: l_ch4
     ! Include/exclude Methane
  LOGICAL :: l_co
     ! Include/exclude Carbon monoxide
  LOGICAL :: l_co2
     ! Include/exclude Carbon dioxide
  LOGICAL :: l_cs
     ! Include/exclude Cesium
  LOGICAL :: l_h2
     ! Include/exclude Hydrogen
  LOGICAL :: l_h2o
     ! Include/exclude Water vapour
  LOGICAL :: l_hcfc22
     ! Include/exclude HCFC22
  LOGICAL :: l_hcn
     ! Include/exclude Hydrogen cyanide
  LOGICAL :: l_he
     ! Include/exclude Helium
  LOGICAL :: l_hfc134a
     ! Include/exclude HFC134a
  LOGICAL :: l_k
     ! Include/exclude Potassium
  LOGICAL :: l_li
     ! Include/exclude Lithium
  LOGICAL :: l_n2
     ! Include/exclude Nitrogen
  LOGICAL :: l_n2o
     ! Include/exclude Dinitrogen oxide
  LOGICAL :: l_na
     ! Include/exclude Sodium
  LOGICAL :: l_nh3
     ! Include/exclude Ammonia
  LOGICAL :: l_o2
     ! Include/exclude Oxygen
  LOGICAL :: l_o3
     ! Include/exclude Ozone
  LOGICAL :: l_rb
     ! Include/exclude Rubidium
  LOGICAL :: l_so2
     ! Include/exclude Sulphur dioxide
  LOGICAL :: l_tio
     ! Include/exclude Titanium oxide
  LOGICAL :: l_vo
     ! Include/exclude Vanadium oxide

  ! Mass mixing ratios of each radiative gas
  ! is constant (.TRUE.) or 
  ! is variable across columns/layers (.FALSE.)
  LOGICAL :: cfc113_well_mixed     
  LOGICAL :: cfc11_well_mixed
  LOGICAL :: cfc12_well_mixed
  LOGICAL :: hcfc22_well_mixed
  LOGICAL :: hfc134a_well_mixed
  LOGICAL :: ch4_well_mixed
  LOGICAL :: co_well_mixed
  LOGICAL :: co2_well_mixed
  LOGICAL :: cs_well_mixed
  LOGICAL :: h2_well_mixed
  LOGICAL :: h2o_well_mixed
  LOGICAL :: hcn_well_mixed
  LOGICAL :: he_well_mixed
  LOGICAL :: k_well_mixed
  LOGICAL :: li_well_mixed
  LOGICAL :: n2_well_mixed
  LOGICAL :: n2o_well_mixed
  LOGICAL :: na_well_mixed
  LOGICAL :: nh3_well_mixed
  LOGICAL :: o2_well_mixed
  LOGICAL :: o3_well_mixed
  LOGICAL :: rb_well_mixed
  LOGICAL :: so2_well_mixed
  LOGICAL :: tio_well_mixed
  LOGICAL :: vo_well_mixed

  LOGICAL :: l_orog
     ! Consider orography (or not)
  LOGICAL :: l_aerosol_mode_sw
     ! MODE aerosols in SW (or not)
  LOGICAL :: l_aerosol_mode_lw
     ! MODE aerosols in LW (or not)

  TYPE(StrDiag) :: sw_diag
     ! Structure holding SW diagnostics in call to runes
  TYPE(StrDiag) :: lw_diag
     ! Structure holding LW diagnostics in call to runes

  ! Arrays to hold the input fields from the LFRic run
  ! (read from netCDF file)

  REAL(RealExt), target, ALLOCATABLE :: sw_heating_rate_rts_nc(:)
     ! SW heating rate
  REAL(RealExt), target, ALLOCATABLE :: lw_heating_rate_rts_nc(:)
     ! LW heating rate
  REAL(RealExt), target, ALLOCATABLE :: sw_up_rts_nc(:)
     ! SW upwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: sw_down_rts_nc(:)
     ! SW downwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: lw_up_rts_nc(:)
     ! LW upwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: lw_down_rts_nc(:)
     ! LW downwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: sw_up_clear_rts_nc(:)
     ! Clear-sky SW upwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: sw_down_clear_rts_nc(:)
     ! Clear-sky SW downwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: lw_up_clear_rts_nc(:)
     ! Clear-sky LW upwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: lw_down_clear_rts_nc(:)
     ! Clear-sky LW downwards flux on radiation levels
  REAL(RealExt), target, ALLOCATABLE :: sw_down_clear_clean_surf_rts_nc(:)
     ! Clear-clean SW downwards surface flux   
  REAL(RealExt), target, ALLOCATABLE :: sw_up_clear_clean_surf_rts_nc(:)
     ! Clear-clean SW upwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: sw_up_clear_clean_toa_rts_nc(:)
     ! Clear-clean SW upwards top-of-atmosphere flux
  REAL(RealExt), target, ALLOCATABLE :: lw_down_clear_clean_surf_rts_nc(:)
     ! Clear-clean LW downwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: lw_up_clear_clean_surf_rts_nc(:)
     ! Clear-clean LW upwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: lw_up_clear_clean_toa_rts_nc(:)
     ! Clear-clean LW upwards top-of-atmosphere flux

  REAL(RealExt), target, ALLOCATABLE :: pressure_in_wth_nc(:)
     ! Pressure in wth space
  REAL(RealExt), target, ALLOCATABLE :: temperature_in_wth_nc(:)
     ! Temperature in wth space
  REAL(RealExt), target, ALLOCATABLE :: d_mass_nc(:)
     ! Mass per square meter of radiation layers
  REAL(RealExt), target, ALLOCATABLE :: rho_in_wth_nc(:)
     ! Density in potential temperature space
  REAL(RealExt), target, ALLOCATABLE :: t_layer_boundaries_nc(:)
     ! Temperature on radiation levels

  REAL(RealExt), target, ALLOCATABLE :: ch4_nc(:)
     ! Mass mixing ratio of Methane
  REAL(RealExt), target, ALLOCATABLE :: co_nc(:)
     ! Mass mixing ratio of Carbon monoxide
  REAL(RealExt), target, ALLOCATABLE :: co2_nc(:)
     ! Mass mixing ratio of Carbon dioxide
  REAL(RealExt), target, ALLOCATABLE :: cs_nc(:)
     ! Mass mixing ratio of Cesium
  REAL(RealExt), target, ALLOCATABLE :: h2_nc(:)
     ! Mass mixing ratio of Hydrogen
  REAL(RealExt), target, ALLOCATABLE :: h2o_nc(:)
     ! Mass mixing ratio of Water vapour
  REAL(RealExt), target, ALLOCATABLE :: hcn_nc(:)
     ! Mass mixing ratio of Hydrogen cyanide
  REAL(RealExt), target, ALLOCATABLE :: he_nc(:)
     ! Mass mixing ratio of Helium
  REAL(RealExt), target, ALLOCATABLE :: potassium_nc(:)
     ! Mass mixing ratio of Potassium
  REAL(RealExt), target, ALLOCATABLE :: li_nc(:)
     ! Mass mixing ratio of Lithium
  REAL(RealExt), target, ALLOCATABLE :: n2_nc(:)
     ! Mass mixing ratio of Nitrogen
  REAL(RealExt), target, ALLOCATABLE :: n2o_nc(:)
     ! Mass mixing ratio of Dinitrogen oxide
  REAL(RealExt), target, ALLOCATABLE :: na_nc(:)
     ! Mass mixing ratio of Sodium
  REAL(RealExt), target, ALLOCATABLE :: nh3_nc(:)
     ! Mass mixing ratio of Ammonia
  REAL(RealExt), target, ALLOCATABLE :: o2_nc(:)
     ! Mass mixing ratio of Oxygen
  REAL(RealExt), target, ALLOCATABLE :: o3_nc(:)
     ! Mass mixing ratio of Ozone
  REAL(RealExt), target, ALLOCATABLE :: rb_nc(:)
     ! Mass mixing ratio of Rubidium
  REAL(RealExt), target, ALLOCATABLE :: so2_nc(:)
     ! Mass mixing ratio field of Sulphur dioxide
  REAL(RealExt), target, ALLOCATABLE :: tio_nc(:)
     ! Mass mixing ratio field of Titanium oxide
  REAL(RealExt), target, ALLOCATABLE :: vo_nc(:)
     ! Mass mixing ratio field of Vanadium oxide

  REAL(RealExt), target, ALLOCATABLE :: tile_fraction_nc(:)
     ! Surface tile fractions
  REAL(RealExt), target, ALLOCATABLE :: tile_temperature_nc(:)
     ! Surface tile temperature
  REAL(RealExt), target, ALLOCATABLE :: tile_lw_albedo_nc(:)
     ! LW tile albedo
  REAL(RealExt), target, ALLOCATABLE :: tile_sw_diffuse_albedo_nc(:)
     ! SW diffuse tile albedo
  REAL(RealExt), target, ALLOCATABLE :: tile_sw_direct_albedo_nc(:)
     ! SW direct tile albedo
  REAL(RealExt), target, ALLOCATABLE :: radiative_cloud_fraction_nc(:)
     ! Large scale cloud fraction
  REAL(RealExt), target, ALLOCATABLE :: liquid_fraction_nc(:)
     ! Liquid cloud fraction field
  REAL(RealExt), target, ALLOCATABLE :: frozen_fraction_nc(:)
     ! Frozen cloud fraction field
  REAL(RealExt), target, ALLOCATABLE :: mcl_nc(:)
     ! Cloud liquid field
  REAL(RealExt), target, ALLOCATABLE :: mci_nc(:)
     ! Cloud ice field
  REAL(RealExt), target, ALLOCATABLE :: n_ice_nc(:)
     ! Ice number concentration
  REAL(RealExt), target, ALLOCATABLE :: cloud_drop_no_conc_nc(:)
     ! Cloud drop number concentration
  REAL(RealExt), target, ALLOCATABLE :: radiative_conv_fraction_nc(:)
     ! Convective cloud fraction
  REAL(RealExt), target, ALLOCATABLE :: conv_liquid_fraction_nc(:)
     ! Convective liquid cloud fraction
  REAL(RealExt), target, ALLOCATABLE :: conv_frozen_fraction_nc(:)
     ! Convective ice cloud fraction
  REAL(RealExt), target, ALLOCATABLE :: conv_liquid_mmr_nc(:)
     ! Convective liquid gridbox MMR
  REAL(RealExt), target, ALLOCATABLE :: conv_frozen_mmr_nc(:)
     ! Convective frozen gridbox MMR
  REAL(RealExt), target, ALLOCATABLE :: sigma_mc_nc(:)
     ! Fractional standard deviation of condensate
  REAL(RealExt), target, ALLOCATABLE :: layer_heat_capacity_nc(:)
     ! Heat capacity of radiation layers
  REAL(RealExt), target, ALLOCATABLE :: sulphuric_nc(:)
     ! Sulphuric acid aerosol
  REAL(RealExt), target, ALLOCATABLE :: cos_zenith_angle_rts_nc(:)
     ! Cosine of the stellar zenith angle
  REAL(RealExt), target, ALLOCATABLE :: lit_fraction_rts_nc(:)
     ! Lit fraction of the time step
  REAL(RealExt), target, ALLOCATABLE :: stellar_irradiance_rts_nc(:)
     ! Stellar irradiance at the planet
  REAL(RealExt), target, ALLOCATABLE :: orographic_correction_rts_nc(:)
     ! Orographic correction
  REAL(RealExt), target, ALLOCATABLE :: rand_seed_real_radinput_nc(:)
     ! Random seed field for cloud generator 
     ! (saved as real in netCDF) 
  REAL(RealExt), target, ALLOCATABLE :: aer_mix_ratio_nc(:)
     ! MODE aerosol mixing ratios
  REAL(RealExt), target, ALLOCATABLE :: aer_sw_absorption_nc(:)
     ! MODE aerosol SW absorption
  REAL(RealExt), target, ALLOCATABLE :: aer_sw_scattering_nc(:)
     ! MODE aerosol SW scattering
  REAL(RealExt), target, ALLOCATABLE :: aer_sw_asymmetry_nc(:)
     ! MODE aerosol SW asymmetry
  REAL(RealExt), target, ALLOCATABLE :: aer_lw_absorption_nc(:)
     ! MODE aerosol LW absorption
  REAL(RealExt), target, ALLOCATABLE :: aer_lw_scattering_nc(:)
     ! MODE aerosol LW scattering
  REAL(RealExt), target, ALLOCATABLE :: aer_lw_asymmetry_nc(:)
     ! MODE aerosol LW asymmetry

  REAL(RealExt), target, ALLOCATABLE :: dummy_nc(:)
     ! Dummy array to point to
     ! for any field that is missing in the netCDF file

  ! Pointers to the output fluxes/heating rates from the LFRic run
  ! (read from netCDF file)

  REAL(RealExt), pointer :: sw_heating_rate_rts(:)
     ! SW heating rate
  REAL(RealExt), pointer :: lw_heating_rate_rts(:)
     ! LW heating rate
  REAL(RealExt), pointer :: sw_up_rts(:)
     ! SW upwards flux on radiation levels
  REAL(RealExt), pointer :: sw_down_rts(:)
     ! SW downwards flux on radiation levels
  REAL(RealExt), pointer :: lw_up_rts(:)
     ! LW upwards flux on radiation levels
  REAL(RealExt), pointer :: lw_down_rts(:)
     ! LW downwards flux on radiation levels
  REAL(RealExt), pointer :: sw_up_clear_rts(:)
     ! Clear-sky SW upwards flux on radiation levels
  REAL(RealExt), pointer :: sw_down_clear_rts(:)
     ! Clear-sky SW downwards flux on radiation levels
  REAL(RealExt), pointer :: lw_up_clear_rts(:)
     ! Clear-sky LW upwards flux on radiation levels
  REAL(RealExt), pointer :: lw_down_clear_rts(:)
     ! Clear-sky LW downwards flux on radiation levels
  REAL(RealExt), pointer :: sw_down_clear_clean_surf_rts(:)
     ! Clear-clean SW downwards surface flux
  REAL(RealExt), pointer :: sw_up_clear_clean_surf_rts(:)
     ! Clear-clean SW upwards surface flux
  REAL(RealExt), pointer :: sw_up_clear_clean_toa_rts(:)
     ! Clear-clean SW upwards top-of-atmosphere flux
  REAL(RealExt), pointer :: lw_down_clear_clean_surf_rts(:)
     ! Clear-clean LW downwards surface flux
  REAL(RealExt), pointer :: lw_up_clear_clean_surf_rts(:)
     ! Clear-clean LW upwards surface flux
  REAL(RealExt), pointer :: lw_up_clear_clean_toa_rts(:)
     ! Clear-clean LW upwards top-of-atmosphere flux

  ! Pointers to the input fields from the LFRic run
  ! (read from netCDF file)
    
  REAL(RealExt), pointer :: pressure_in_wth(:)
     ! Pressure in wth space
  REAL(RealExt), pointer :: temperature_in_wth(:)
     ! Temperature in wth space
  REAL(RealExt), pointer :: d_mass(:)
     ! Mass per square meter of radiation layers
  REAL(RealExt), pointer :: rho_in_wth(:)
     ! Density in potential temperature space
  REAL(RealExt), pointer :: t_layer_boundaries(:)
     ! Temperature on radiation levels

  REAL(RealExt), pointer :: ch4(:)
     ! Mass mixing ratio of Methane
  REAL(RealExt), pointer :: co(:)
     ! Mass mixing ratio of Carbon monoxide
  REAL(RealExt), pointer :: co2(:)
     ! Mass mixing ratio of Carbon dioxide
  REAL(RealExt), pointer :: cs(:)
     ! Mass mixing ratio of Cesium
  REAL(RealExt), pointer :: h2(:)
     ! Mass mixing ratio of Hydrogen
  REAL(RealExt), pointer :: h2o(:)
     ! Mass mixing ratio of Water vapour
  REAL(RealExt), pointer :: hcn(:)
     ! Mass mixing ratio of Hydrogen cyanide
  REAL(RealExt), pointer :: he(:)
     ! Mass mixing ratio of Helium
  REAL(RealExt), pointer :: potassium(:)
     ! Mass mixing ratio of Potassium
  REAL(RealExt), pointer :: li(:)
     ! Mass mixing ratio of Lithium
  REAL(RealExt), pointer :: n2(:)
     ! Mass mixing ratio of Nitrogen
  REAL(RealExt), pointer :: n2o(:)
     ! Mass mixing ratio of Dinitrogen oxide
  REAL(RealExt), pointer :: na(:)
     ! Mass mixing ratio of Sodium
  REAL(RealExt), pointer :: nh3(:)
     ! Mass mixing ratio of Ammonia
  REAL(RealExt), pointer :: o2(:)
     ! Mass mixing ratio of Oxygen
  REAL(RealExt), pointer :: o3(:)
     ! Mass mixing ratio of Ozone
  REAL(RealExt), pointer :: rb(:)
     ! Mass mixing ratio of Rubidium
  REAL(RealExt), pointer :: so2(:)
     ! Mass mixing ratio of Sulphur dioxide
  REAL(RealExt), pointer :: tio(:)
     ! Mass mixing ratio of Titanium oxide
  REAL(RealExt), pointer :: vo(:)
     ! Mass mixing ratio of Vanadium oxide

  REAL(RealExt), pointer :: tile_fraction(:)
     ! Surface tile fractions
  REAL(RealExt), pointer :: tile_temperature(:)
     ! Surface tile temperature
  REAL(RealExt), pointer :: tile_lw_albedo(:)
     ! LW tile albedos
  REAL(RealExt), pointer :: tile_sw_diffuse_albedo(:)
     ! Diffuse SW tile albedo
  REAL(RealExt), pointer :: tile_sw_direct_albedo(:)
     ! Direct SW tile albedo
  REAL(RealExt), pointer :: radiative_cloud_fraction(:)
     ! Large scale cloud fraction
  REAL(RealExt), pointer :: liquid_fraction(:)
     ! Liquid cloud fraction field
  REAL(RealExt), pointer :: frozen_fraction(:)
     ! Frozen cloud fraction field
  REAL(RealExt), pointer :: mcl(:)
     ! Cloud liquid field
  REAL(RealExt), pointer :: mci(:)
     ! Cloud ice field
  REAL(RealExt), pointer :: n_ice(:)
     ! Ice number concentration
  REAL(RealExt), pointer :: cloud_drop_no_conc(:)
     ! Cloud drop number concentration
  REAL(RealExt), pointer :: radiative_conv_fraction(:)
     ! Convective cloud fraction
  REAL(RealExt), pointer :: conv_liquid_fraction(:)
     ! Convective liquid cloud fraction
  REAL(RealExt), pointer :: conv_frozen_fraction(:)
     ! Convective frozen cloud fraction
  REAL(RealExt), pointer :: conv_liquid_mmr(:)
     ! Convective liquid gridbox MMR
  REAL(RealExt), pointer :: conv_frozen_mmr(:)
     ! Convective frozen gridbox MMR
  REAL(RealExt), pointer :: sigma_mc(:)
     ! Fractional standard deviation of condensate
  REAL(RealExt), pointer :: layer_heat_capacity(:)
     ! Heat capacity of radiation layers
  REAL(RealExt), pointer :: sulphuric(:)
     ! Sulphuric acid aerosol
  REAL(RealExt), pointer :: cos_zenith_angle_rts(:)
     ! Cosine of the stellar zenith angle
  REAL(RealExt), pointer :: lit_fraction_rts(:)
     ! Lit fraction of the time step
  REAL(RealExt), pointer :: stellar_irradiance_rts(:)
     ! Stellar irradiance at the planet
  REAL(RealExt), pointer :: orographic_correction_rts(:)
     ! Orographic correction
  REAL(RealExt), pointer :: rand_seed_real_radinput(:)
     ! Random seed field for cloud generator 
     ! (saved as real in netCDF) 
  REAL(RealExt), pointer :: aer_mix_ratio(:)
     ! MODE aerosol mixing ratios
  REAL(RealExt), pointer :: aer_sw_absorption(:)
     ! MODE aerosol SW absorption
  REAL(RealExt), pointer :: aer_sw_scattering(:)
     ! MODE aerosol SW scattering
  REAL(RealExt), pointer :: aer_sw_asymmetry(:)
     ! MODE aerosol SW asymmetry
  REAL(RealExt), pointer :: aer_lw_absorption(:)
     ! MODE aerosol LW absorption
  REAL(RealExt), pointer :: aer_lw_scattering(:)
     ! MODE aerosol LW scattering
  REAL(RealExt), pointer :: aer_lw_asymmetry(:)
     ! MODE aerosol LW asymmetry

  ! Arrays to hold the output flux/heating rate fields
  ! from the offline run

  REAL(RealExt), target, ALLOCATABLE :: sw_heating_rate(:,:)
     ! SW heating rate
  REAL(RealExt), target, ALLOCATABLE :: lw_heating_rate(:,:)
     ! LW heating rate
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_up(:,:)
     ! SW upwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_down(:,:)
     ! SW downwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_up_clear(:,:)
     ! Clear-sky SW upwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_down_clear(:,:)
     ! Clear-sky SW downwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_down_clear_clean_surf(:)
     ! Clear-clean SW downwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_up_clear_clean_surf(:)
     ! Clear-clean SW upwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: sw_flux_up_clear_clean_toa(:)
     ! Clear-clean SW upwards top-of-atmosphere flux
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_up(:,:)
     ! LW upwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_down(:,:)
     ! LW downwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_up_clear(:,:)
     ! Clear-sky LW upwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_down_clear(:,:)
     ! Clear-sky LW downwards flux on levels
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_down_clear_clean_surf(:)
     ! Clear-clean LW downwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_up_clear_clean_surf(:)
     ! Clear-clean LW upwards surface flux
  REAL(RealExt), target, ALLOCATABLE :: lw_flux_up_clear_clean_toa(:)
     ! Clear-clean LW upwards top-of-atmosphere flux

  ! Local variables
  
  INTEGER, ALLOCATABLE :: rand_seed(:)
     ! Random seed field for cloud generator
     ! (required as integer in runes interface)
  INTEGER :: n_profile_list
     ! Number of columns with given equal cloud top
  INTEGER, ALLOCATABLE :: profile_list(:)
     ! List to hold indices of columns 
     ! with given equal cloud top
  INTEGER, ALLOCATABLE :: n_cloud_layer(:)
     ! Number of cloud layers
       
  REAL(RealExt), target, ALLOCATABLE :: dummy_arr_out(:)
     ! Generic array for fields to be written to netCDF
     ! after packing them in 1D
  
  LOGICAL :: l_debug
     ! Flag to turn on/off write statements for debugging
  LOGICAL :: l_ascii
     ! Flag to turn on/off flux/heating rates output in ascii format
     ! in addition to netCDF format

  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'runes_nc'

  CHARACTER(LEN=200) :: conf_file
     ! Configuration file name
  CHARACTER(LEN=200) :: input_file
     ! Input (netCDF) file
  CHARACTER(LEN=200) :: output_file
     ! Output file
  CHARACTER(LEN=200) :: spectral_dir
     ! Spectral file directory

  ! Namelist unit
  INTEGER :: fn_conf
  
  ! Units for ascii output files
  INTEGER :: iu_out_sw_lfric
  INTEGER :: iu_out_sw_offline
  INTEGER :: iu_out_lw_lfric
  INTEGER :: iu_out_lw_offline

  ! Local variables for netCDF support
  INTEGER :: ncid
  INTEGER :: ierr
  INTEGER :: ios

  INTEGER :: n_dimension
  INTEGER :: dimension_size
  INTEGER :: dimids
  INTEGER :: dimout(3)

  CHARACTER(LEN=200) :: varname
  CHARACTER(LEN=200) :: dimension_name
  CHARACTER(LEN=200) :: arg
  
  ! Loop variables
  
  INTEGER :: i, j, k, l, ic

  INTEGER :: i_time

  INTEGER :: nn, wth_start, wth_end, kk

  ! Local variables analog to those for the LFRic kernels
  INTEGER :: twod_1
  INTEGER :: twod_last
  INTEGER :: wth_1
  INTEGER :: wth_last
  INTEGER :: flux_0
  INTEGER :: flux_last
  INTEGER :: tile_1
  INTEGER :: tile_last
  INTEGER :: rtile_sw_1
  INTEGER :: rtile_sw_last
  INTEGER :: rtile_lw_1
  INTEGER :: rtile_lw_last
  INTEGER :: mode_1
  INTEGER :: mode_last
  INTEGER :: rmode_sw_1
  INTEGER :: rmode_sw_last
  INTEGER :: rmode_lw_1
  INTEGER :: rmode_lw_last

  ! Default settings for output

  l_ascii = .false.
  l_debug = .false.

  IF (l_debug) THEN
    WRITE(*,*) 'PROGRAM runes_nc start: '
  END IF ! l_debug


  ! Read command line arguments

  IF (l_debug) THEN
    WRITE(*,*) 'command_argument_count: ', command_argument_count()
  END IF ! l_debug
  ic=0
  DO
    ic=ic+1
    IF (l_debug) THEN
      WRITE(*,*) 'ic: ', ic
    END IF ! l_debug
    IF (ic > command_argument_count()) EXIT
    CALL get_command_argument(ic, arg)
      SELECT CASE(arg)
        CASE ('-conf_file')
          CALL get_command_argument(ic+1, conf_file)
          IF (l_debug) THEN
            WRITE(*,*) 'conf_file: ', TRIM(conf_file)
          END IF ! l_debug
          ic=ic+1
        CASE ('-input_file')
          CALL get_command_argument(ic+1, input_file)
          IF (l_debug) THEN
            WRITE(*,*) 'input_file: ', TRIM(input_file)
          END IF ! l_debug
          ic=ic+1
        CASE ('-spectral_dir')
          CALL get_command_argument(ic+1, spectral_dir)
          IF (l_debug) THEN
            WRITE(*,*) 'spectral_file_dir: ', TRIM(spectral_dir)
          END IF ! l_debug
          ic=ic+1
        CASE ('-ascii')
          l_ascii=.true.
          IF (l_debug) THEN
            WRITE(*,*) 'write ascii output'
          END IF ! l_debug
          ic=ic+1
        CASE ('-debug')
          l_debug=.true.
          IF (l_debug) THEN
            WRITE(*,*) 'write debug output'
          END IF ! l_debug
          ic=ic+1
      END SELECT
  END DO
  
  ! The fluxes/heating rates calculated offline 
  ! are added to the netCDF input file

  output_file=input_file

  
  ! Read namelist element defaults

  CALL set_nml_default

  ! Update with the actual namelist entries

  OPEN(newunit=fn_conf, file=TRIM(conf_file), &
    status='old', action="read", iostat=ios)
  IF (ios.EQ.0) THEN
    READ(fn_conf, nml=aerosol, iostat=ierr)
    REWIND(fn_conf)
    READ(fn_conf, nml=radiation, iostat=ierr)
    IF (ierr.NE.0) THEN
      WRITE(cmessage,'(A)') 'no radiation section in the conf namelist file'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage) 
    END IF
    REWIND(fn_conf)
    READ(fn_conf, nml=radiative_gases, iostat=ierr)
  ELSE
    WRITE(cmessage,'(A)') 'problem reading configuration namelist file'
    ierr = i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage) 
  END IF
  CLOSE(fn_conf)


  ! Settings of radiative gases

  l_cfc113 = .true.
  cfc113_well_mixed = .true.
  SELECT CASE ( TRIM(cfc113_rad_opt) )
    CASE ('off')
      l_cfc113 = .false.
      cfc113_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (cfc113_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'cfc113_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      cfc113_well_mixed = .false.
    CASE ('prognostic')
      cfc113_well_mixed = .false.
      WRITE(cmessage,'(A)') 'cfc113 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_cfc11 = .true.
  cfc11_well_mixed = .true.
  SELECT CASE ( TRIM(cfc11_rad_opt) )
    CASE ('off')
      l_cfc11 = .false.
      cfc11_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (cfc11_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'cfc11_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      cfc11_well_mixed = .false.
    CASE ('prognostic')
      cfc11_well_mixed = .false.
      WRITE(cmessage,'(A)') 'cfc11 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_cfc12 = .true.
  cfc12_well_mixed = .true.
  SELECT CASE ( TRIM(cfc12_rad_opt) )
    CASE ('off')
      l_cfc12 = .false.
      cfc12_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (cfc12_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'cfc12_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      cfc12_well_mixed = .false.
    CASE ('prognostic')
      cfc12_well_mixed = .false.
      WRITE(cmessage,'(A)') 'cfc12 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_ch4 = .true.
  ch4_well_mixed = .true.
  SELECT CASE ( TRIM(ch4_rad_opt) )
    CASE ('off')
      l_ch4 = .false.
      ch4_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (ch4_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'ch4_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      ch4_well_mixed = .false.
    CASE ('prognostic')
      ch4_well_mixed = .false.
      WRITE(cmessage,'(A)') 'ch4 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_co2 = .true.
  co2_well_mixed = .true.
  SELECT CASE ( TRIM(co2_rad_opt) )
    CASE ('off')
      l_co2 = .false.
      co2_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (co2_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'co2_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      co2_well_mixed = .false.
    CASE ('prognostic')
      co2_well_mixed = .false.
      WRITE(cmessage,'(A)') 'co2 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_co = .true.
  co_well_mixed = .true.
  SELECT CASE ( TRIM(co_rad_opt) )
    CASE ('off')
      l_co = .false.
      co_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (co_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'co_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      co_well_mixed = .false.
    CASE ('prognostic')
      co_well_mixed = .false.
      WRITE(cmessage,'(A)') 'co prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_h2 = .true.
  h2_well_mixed = .true.
  SELECT CASE ( TRIM(h2_rad_opt) )
    CASE ('off')
      l_h2 = .false.
      h2_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (h2_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'h2_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      h2_well_mixed = .false.
    CASE ('prognostic')
      h2_well_mixed = .false.
      WRITE(cmessage,'(A)') 'h2 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_h2o = .true.
  h2o_well_mixed = .true.
  SELECT CASE ( TRIM(h2o_rad_opt) )
    CASE ('off')
      l_h2o = .false.
      h2o_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (h2o_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'h2o_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      h2o_well_mixed = .false.
    CASE ('prognostic')
      h2o_well_mixed = .false.
    CASE DEFAULT
  END SELECT

  l_hcfc22 = .true.
  hcfc22_well_mixed = .true.
  SELECT CASE ( TRIM(hcfc22_rad_opt) )
    CASE ('off')
      l_hcfc22 = .false. 
      hcfc22_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (hcfc22_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'hcfc22_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      hcfc22_well_mixed = .false.
    CASE ('prognostic')
      hcfc22_well_mixed = .false.
      WRITE(cmessage,'(A)') 'hcfc22 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_hcn = .true.
  hcn_well_mixed = .true.
  SELECT CASE ( TRIM(hcn_rad_opt) )
    CASE ('off')
      l_hcn = .false.
      hcn_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (hcn_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'hcn_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      hcn_well_mixed = .false.
    CASE ('prognostic')
      hcn_well_mixed = .false.
      WRITE(cmessage,'(A)') 'hcn prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_he = .true.
  he_well_mixed = .true.
  SELECT CASE ( TRIM(he_rad_opt) )
    CASE ('off')
      l_he = .false.
      he_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (he_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'he_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      he_well_mixed = .false.
    CASE ('prognostic')
      he_well_mixed = .false.
      WRITE(cmessage,'(A)') 'he prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_hfc134a = .true.
  hfc134a_well_mixed = .true.
  SELECT CASE ( TRIM(hfc134a_rad_opt) )
    CASE ('off')
      l_hfc134a = .false.
      hfc134a_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (hfc134a_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'hfc134a_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      hfc134a_well_mixed = .false.
    CASE ('prognostic')
      hfc134a_well_mixed = .false.
      WRITE(cmessage,'(A)') 'hfc134a prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_n2 = .true.
  n2_well_mixed = .true.
  SELECT CASE ( TRIM(n2_rad_opt) )
    CASE ('off')
      l_n2 = .false.
      n2_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (n2_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'n2_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      n2_well_mixed = .false.
    CASE ('prognostic')
      n2_well_mixed = .false.
      WRITE(cmessage,'(A)') 'n2 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_n2o = .true.
  n2o_well_mixed = .true.
  SELECT CASE ( TRIM(n2o_rad_opt) )
    CASE ('off')
      l_n2o = .false.
      n2o_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (n2o_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'n2o_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      n2o_well_mixed = .false.
    CASE ('prognostic')
      n2o_well_mixed = .false.
      WRITE(cmessage,'(A)') 'n2o prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_nh3 = .true.
  nh3_well_mixed = .true.
  SELECT CASE ( TRIM(nh3_rad_opt) )
    CASE ('off')
      l_nh3 = .false.
      nh3_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (nh3_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'nh3_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      nh3_well_mixed = .false.
    CASE ('prognostic')
      nh3_well_mixed = .false.
      WRITE(cmessage,'(A)') 'nh3 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_o2 = .true.
  o2_well_mixed = .true.
  SELECT CASE ( TRIM(o2_rad_opt) )
    CASE ('off')
      l_o2 = .false.
      o2_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (o2_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'o2_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      o2_well_mixed = .false.
    CASE ('prognostic')
      o2_well_mixed = .false.
      WRITE(cmessage,'(A)') 'o2 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  l_o3 = .true.
  o3_well_mixed = .true.
  SELECT CASE ( TRIM(o3_rad_opt) )
    CASE ('off')
      l_o3 = .false.
      o3_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (o3_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'o3_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      o3_well_mixed = .false.
    CASE ('prognostic')
      o3_well_mixed = .false.
    CASE DEFAULT
  END SELECT

  l_so2 = .true.
  so2_well_mixed = .true.
  SELECT CASE ( TRIM(so2_rad_opt) )
    CASE ('off')
      l_so2 = .false.
      so2_mix_ratio = 0.0_RealExt
    CASE ('constant')
      IF (so2_mix_ratio <= 0.0_RealExt) THEN
        WRITE(cmessage,'(A)') 'so2_mix_ratio not set'
        ierr = i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    CASE ('ancil')
      so2_well_mixed = .false.
    CASE ('prognostic')
      so2_well_mixed = .false.
      WRITE(cmessage,'(A)') 'so2 prognostic not implemented'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    CASE DEFAULT
  END SELECT

  ! Set logical for the orographic correction
  if (TRIM(topography) == 'flat') then
    l_orog = .false.
  else
    l_orog = .true.
  end if

  ! Settings of cloud properties

  select case (TRIM(scatter_method_lw))
  case ('full')
    i_scatter_method_lw = ip_scatter_full
  case ('none')
    i_scatter_method_lw = ip_scatter_none
  case ('approx')
    i_scatter_method_lw = ip_scatter_approx
  case ('hybrid')
    i_scatter_method_lw = ip_scatter_hybrid
  case default
    i_scatter_method_lw = ip_scatter_full
  end select

  select case (TRIM(cloud_representation))
  case ('no_cloud')
    i_cloud_representation = ip_cloud_representation_off
  case ('liquid_and_ice')
    i_cloud_representation = ip_cloud_representation_ice_water
  case ('combined')
    i_cloud_representation = ip_cloud_representation_combine_ice_water
  case ('conv_strat_liq_ice')
    i_cloud_representation = ip_cloud_representation_csiw
  case ('split')
    i_cloud_representation = ip_cloud_representation_split_ice_water
  case default
    i_cloud_representation = ip_cloud_representation_off
  end select

  select case (TRIM(cloud_overlap))
  case ('maximum_random')
    i_overlap = ip_overlap_max_random
  case ('random')
    i_overlap = ip_overlap_random
  case ('exponential_random')
    i_overlap = ip_overlap_exponential_random
  case default
    i_overlap = ip_overlap_max_random
  end select

  select case (TRIM(cloud_inhomogeneity))
  case ('homogeneous')
    i_inhom = ip_inhom_homogeneous
  case ('scaling')
    i_inhom = ip_inhom_scaling
  case ('mcica')
    i_inhom = ip_inhom_mcica
  case ('cairns')
    i_inhom = ip_inhom_cairns
  case ('tripleclouds')
    i_inhom = ip_inhom_tripleclouds_2019
  case default
    i_inhom = ip_inhom_homogeneous
  end select

  select case (TRIM(cloud_entrapment))
  case ('max')
    i_cloud_entrapment = ip_cloud_entrapment_max
  case ('zero')
    i_cloud_entrapment = ip_cloud_entrapment_zero
  case default
    i_cloud_entrapment = ip_cloud_entrapment_zero
  end select

  select case (TRIM(droplet_effective_radius))
  case ('constant')
    i_drop_re = ip_droplet_re_constant
  case ('liu')
    i_drop_re = ip_droplet_re_liu
  case default
    i_drop_re = ip_droplet_re_default
  end select

  ! Settings of aerosols

  IF (l_radaer .or. easyaerosol_sw) THEN
    l_aerosol_mode_sw = .true.
  ELSE
    l_aerosol_mode_sw = .false.
  END IF

  IF (l_radaer .or. easyaerosol_lw) THEN
    l_aerosol_mode_lw = .true.
  ELSE
    l_aerosol_mode_lw = .false.
  END IF
  

  ! READ netCDF file (LFRic output)

  IF (l_debug) THEN
    WRITE(*,*) 'Reading from netCDF file: ', TRIM(input_file)
  END IF ! l_debug
  
  ierr=NF90_OPEN(Trim(input_file),NF90_NOWRITE,ncid)
  IF (l_debug) THEN
    WRITE(*,*) 'NF90_OPEN: '
    WRITE(*,*) 'ierr: ', ierr
    WRITE(*,*) 'ncid: ', ncid
  END IF ! l_debug
  IF (ierr.NE.0) THEN
    WRITE(*,*) 'input file: ',Trim(input_file),' does not exist'
    STOP
  END IF

  ierr=NF90_INQUIRE(ncid,ndimensions=n_dimension)
  IF (l_debug) THEN
    WRITE(*,*) 'NF90_INQUIRE: '
    WRITE(*,*) 'ierr: ', ierr
    WRITE(*,*) 'ncid: ', ncid
    WRITE(*,*) 'n_dimension: ', n_dimension
  END IF ! l_debug

  DO i=1,n_dimension

    IF (l_debug) THEN
    WRITE(*,*) 'i (of n_dimension): '
    END IF ! l_debug
    
    ierr=NF90_INQUIRE_DIMENSION(ncid,i,name=dimension_name,len=dimension_size)
    IF (l_debug) THEN
      WRITE(*,*) 'NF90_INQUIRE_DIMENSION: '
      WRITE(*,*) 'ierr: ', ierr
      WRITE(*,*) 'ncid: ', ncid
      WRITE(*,*) 'dimension_name: ', Trim(dimension_name)
      WRITE(*,*) 'dimension_size: ', dimension_size
    END IF ! l_debug

    IF (Trim(dimension_name)=='nMesh2d_face') THEN
       n_profile=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_profile: ', n_profile
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='full_levels') THEN
       n_full_levels=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_full_levels: ', n_full_levels
       END IF ! l_debug
       nlayers = n_full_levels - 1
    END IF
    IF (Trim(dimension_name)=='surface_tiles') THEN
       n_surf_tile=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_surf_tile: ', n_surf_tile
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='lw_bands_surface_tiles') THEN
       n_lw_band_tile=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_lw_band_tile: ', n_lw_band_tile
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='sw_bands_surface_tiles') THEN
       n_sw_band_tile=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_sw_band_tile: ', n_sw_band_tile
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='aero_modes') THEN
       n_aero=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_aero: ', n_aero
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='lw_bands_aero_modes') THEN
       n_lw_band_aero=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_lw_band_aero: ', n_lw_band_aero
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='sw_bands_aero_modes') THEN
       n_sw_band_aero=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_sw_band_aero: ', n_sw_band_aero
       END IF ! l_debug
    END IF
    IF (Trim(dimension_name)=='time') THEN
       n_time=dimension_size
       IF (l_debug) THEN
         WRITE(*,*) 'Found: ', Trim(dimension_name)
         WRITE(*,*) 'Assigned to n_time: ', n_time
       END IF ! l_debug
    END IF

  END DO

  ierr=NF90_CLOSE(ncid)
  IF (l_debug) THEN
    WRITE(*,*) 'NF90_CLOSE: '
    WRITE(*,*) 'ierr: ', ierr
    WRITE(*,*) 'ncid: ', ncid
  END IF ! l_debug

  IF (l_debug) THEN
    WRITE(*,*) 'nlayers: ', nlayers
  END IF ! l_debug


  ! Work out appropriate indices of the fields 
  ! in the runes interface

  twod_1        = 1
  twod_last     = n_profile
  wth_1         = 2
  wth_last      = n_profile*(nlayers+1)
  flux_0        = 1
  flux_last     = n_profile*(nlayers+1)
  tile_1        = 1
  tile_last     = n_profile*n_surf_tile
  rtile_sw_1    = 1
  rtile_sw_last = n_profile*n_sw_band_tile
  rtile_lw_1    = 1
  rtile_lw_last = n_profile*n_lw_band_tile
  mode_1        = 1
  mode_last     = n_profile*(nlayers+1)*n_aero
  rmode_sw_1    = 1
  rmode_sw_last = n_profile*(nlayers+1)*n_sw_band_aero
  rmode_lw_1    = 1
  rmode_lw_last = n_profile*(nlayers+1)*n_lw_band_aero

  n_aer_mode_sw = n_aero
  n_aer_mode_lw = n_aero

  ! Allocate the fields to hold the fluxes/heating rates
  ! calculated by the offline driver

  ALLOCATE(sw_heating_rate(0:nlayers, 1:n_profile))
  ALLOCATE(lw_heating_rate(0:nlayers, 1:n_profile))
  ALLOCATE(sw_flux_up(0:nlayers, 1:n_profile))
  ALLOCATE(sw_flux_down(0:nlayers, 1:n_profile))
  ALLOCATE(sw_flux_up_clear(0:nlayers, 1:n_profile))
  ALLOCATE(sw_flux_down_clear(0:nlayers, 1:n_profile))
  ALLOCATE(lw_flux_up(0:nlayers, 1:n_profile))
  ALLOCATE(lw_flux_down(0:nlayers, 1:n_profile))
  ALLOCATE(lw_flux_up_clear(0:nlayers, 1:n_profile))
  ALLOCATE(lw_flux_down_clear(0:nlayers, 1:n_profile))

  ALLOCATE(sw_flux_down_clear_clean_surf(1:n_profile))
  ALLOCATE(sw_flux_up_clear_clean_surf(1:n_profile))
  ALLOCATE(sw_flux_up_clear_clean_toa(1:n_profile))
  ALLOCATE(lw_flux_down_clear_clean_surf(1:n_profile))
  ALLOCATE(lw_flux_up_clear_clean_surf(1:n_profile))
  ALLOCATE(lw_flux_up_clear_clean_toa(1:n_profile))

  ! Allocate additional fields required for the offline driver call

  ALLOCATE(rand_seed(1:n_profile))
  ALLOCATE(n_cloud_layer(1:n_profile))
  
  ! Allocate the dummy array to a sufficient size

  ALLOCATE(dummy_nc(n_profile*(nlayers+1) &
    *max(1,max(n_sw_band_aero,n_lw_band_aero))))
  
  ALLOCATE(dummy_arr_out(n_profile*(nlayers+1)))

  ! Read in spectral files once at the beginning of a run
  CALL set_spectrum(                        &
    spectrum_name = 'sw',                   &
    spectral_file = TRIM(spectral_dir)//TRIM(spectral_file_sw), &
    l_cfc11       = l_cfc11,                &
    l_cfc113      = l_cfc113,               &
    l_cfc12       = l_cfc12,                &
    l_ch4         = l_ch4,                  &
    l_co          = l_co,                   &
    l_co2         = l_co2,                  &
    l_cs          = l_cs,                   &
    l_h2          = l_h2,                   &
    l_h2o         = l_h2o,                  &
    l_hcfc22      = l_hcfc22,               &
    l_hcn         = l_hcn,                  &
    l_he          = l_he,                   &
    l_hfc134a     = l_hfc134a,              &
    l_k           = l_k,                    &
    l_li          = l_li,                   &
    l_n2          = l_n2,                   &
    l_n2o         = l_n2o,                  &
    l_na          = l_na,                   &
    l_nh3         = l_nh3,                  &
    l_o2          = l_o2,                   &
    l_o3          = l_o3,                   &
    l_rb          = l_rb,                   &
    l_so2         = l_so2,                  &
    l_tio         = l_tio,                  &
    l_vo          = l_vo )

  CALL set_spectrum(                        &
    spectrum_name = 'lw',                   &
    spectral_file = TRIM(spectral_dir)//TRIM(spectral_file_lw), &
    l_cfc11       = l_cfc11,                &
    l_cfc113      = l_cfc113,               &
    l_cfc12       = l_cfc12,                &
    l_ch4         = l_ch4,                  &
    l_co          = l_co,                   &
    l_co2         = l_co2,                  &
    l_cs          = l_cs,                   &
    l_h2          = l_h2,                   &
    l_h2o         = l_h2o,                  &
    l_hcfc22      = l_hcfc22,               &
    l_hcn         = l_hcn,                  &
    l_he          = l_he,                   &
    l_hfc134a     = l_hfc134a,              &
    l_k           = l_k,                    &
    l_li          = l_li,                   &
    l_n2          = l_n2,                   &
    l_n2o         = l_n2o,                  &
    l_na          = l_na,                   &
    l_nh3         = l_nh3,                  &
    l_o2          = l_o2,                   &
    l_o3          = l_o3,                   &
    l_rb          = l_rb,                   &
    l_so2         = l_so2,                  &
    l_tio         = l_tio,                  &
    l_vo          = l_vo )    

  CALL get_spectrum( &
    spectrum_name = 'sw', &
    n_band = n_sw_band )

  CALL get_spectrum( &
    spectrum_name = 'lw', &
    n_band = n_lw_band )

  CALL set_mcica(TRIM(spectral_dir)//mcica_data_file , 'sw', 'lw')


  ! Open files for ascii output
  
  IF (l_ascii) THEN
    OPEN(newunit=iu_out_sw_lfric, file='sw_fluxes_lfric.txt', &
      status='unknown', iostat=ios)
    IF (ios.NE.0) THEN
      WRITE(cmessage,'(A)') 'problem opening ascii file for sw lfric fluxes'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF
    OPEN(newunit=iu_out_sw_offline, file='sw_fluxes_offline.txt', &
      status='unknown', iostat=ios)
    IF (ios.NE.0) THEN
      WRITE(cmessage,'(A)') 'problem opening ascii file for sw offline fluxes'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF
    OPEN(newunit=iu_out_lw_lfric, file='lw_fluxes_lfric.txt', &
      status='unknown', iostat=ios)
    IF (ios.NE.0) THEN
      WRITE(cmessage,'(A)') 'problem opening ascii file for lw lfric fluxes'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF
    OPEN(newunit=iu_out_lw_offline, file='lw_fluxes_offline.txt', &
      status='unknown', iostat=ios)
    IF (ios.NE.0) THEN
      WRITE(cmessage,'(A)') 'problem opening ascii file for lw offline fluxes'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF
  END IF ! l_ascii 


  ! The outer loop is over time steps

  DO i_time=1, n_time

    IF (l_debug) THEN
      WRITE(*,*) 'time_counter: ', i_time
    END IF ! l_debug

    ! Read the fields from the netCDF file
    ! (one per call to read routine)

    CALL read_nc(input_file,'pressure_in_layers_radinput',i_time,pressure_in_wth_nc)
    IF (ALLOCATED(pressure_in_wth_nc)) THEN
      pressure_in_wth(wth_1:wth_last) => pressure_in_wth_nc(wth_1:wth_last)
    ELSE
      pressure_in_wth(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'temperature_in_layers_radinput',i_time,temperature_in_wth_nc)
    IF (ALLOCATED(temperature_in_wth_nc)) THEN
      temperature_in_wth(wth_1:wth_last) => temperature_in_wth_nc(wth_1:wth_last)
    ELSE
      temperature_in_wth(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'dry_mass_in_layers_radinput',i_time,d_mass_nc)
    IF (ALLOCATED(d_mass_nc)) THEN
      d_mass(wth_1:wth_last) => d_mass_nc(wth_1:wth_last)
    ELSE
      d_mass(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'dry_density_in_layers_radinput',i_time,rho_in_wth_nc)
    IF (ALLOCATED(rho_in_wth_nc)) THEN
      rho_in_wth(wth_1:wth_last) => rho_in_wth_nc(wth_1:wth_last)
    ELSE
      rho_in_wth(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'t_layer_boundaries_radinput',i_time,t_layer_boundaries_nc)
    IF (ALLOCATED(rho_in_wth_nc)) THEN
      t_layer_boundaries(flux_0:flux_last) => t_layer_boundaries_nc(flux_0:flux_last)
    ELSE
      t_layer_boundaries(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'ch4_mmr_radinput',i_time,ch4_nc)
    IF (ALLOCATED(ch4_nc)) THEN
      ch4(wth_1:wth_last) => ch4_nc(wth_1:wth_last)
    ELSE
      ch4(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'co_mmr_radinput',i_time,co_nc)
    IF (ALLOCATED(co_nc)) THEN
      co(wth_1:wth_last) => co_nc(wth_1:wth_last)
    ELSE
      co(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'co2_mmr_radinput',i_time,co2_nc)
    IF (ALLOCATED(co2_nc)) THEN
      co2(wth_1:wth_last) => co2_nc(wth_1:wth_last)
    ELSE
      co2(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'cs_mmr_radinput',i_time,cs_nc)
    IF (ALLOCATED(cs_nc)) THEN
      cs(wth_1:wth_last) => cs_nc(wth_1:wth_last)
    ELSE
      cs(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'h2_mmr_radinput',i_time,h2_nc)
    IF (ALLOCATED(h2_nc)) THEN
      h2(wth_1:wth_last) => h2_nc(wth_1:wth_last)
    ELSE
      h2(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'h2o_mmr_radinput',i_time,h2o_nc)
    IF (ALLOCATED(h2o_nc)) THEN
      h2o(wth_1:wth_last) => h2o_nc(wth_1:wth_last)
    ELSE
      h2o(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'hcn_mmr_radinput',i_time,hcn_nc)
    IF (ALLOCATED(hcn_nc)) THEN
      hcn(wth_1:wth_last) => hcn_nc(wth_1:wth_last)
    ELSE
      hcn(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'he_mmr_radinput',i_time,he_nc)
    IF (ALLOCATED(he_nc)) THEN
      he(wth_1:wth_last) => he_nc(wth_1:wth_last)
    ELSE
      he(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'potassium_mmr_radinput',i_time,potassium_nc)
    IF (ALLOCATED(potassium_nc)) THEN
      potassium(wth_1:wth_last) => potassium_nc(wth_1:wth_last)
    ELSE
      potassium(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'li_mmr_radinput',i_time,li_nc)
    IF (ALLOCATED(li_nc)) THEN
      li(wth_1:wth_last) => li_nc(wth_1:wth_last)
    ELSE
      li(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'n2_mmr_radinput',i_time,n2_nc)
    IF (ALLOCATED(n2_nc)) THEN
      n2(wth_1:wth_last) => n2_nc(wth_1:wth_last)
    ELSE
      n2(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'n2o_mmr_radinput',i_time,n2o_nc)
    IF (ALLOCATED(n2o_nc)) THEN
      n2o(wth_1:wth_last) => n2o_nc(wth_1:wth_last)
    ELSE
      n2o(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'na_mmr_radinput',i_time,na_nc)
    IF (ALLOCATED(na_nc)) THEN
      na(wth_1:wth_last) => na_nc(wth_1:wth_last)
    ELSE
      na(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'nh3_mmr_radinput',i_time,nh3_nc)
    IF (ALLOCATED(nh3_nc)) THEN
      nh3(wth_1:wth_last) => nh3_nc(wth_1:wth_last)
    ELSE
      nh3(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'o2_mmr_radinput',i_time,o2_nc)
    IF (ALLOCATED(o2_nc)) THEN
      o2(wth_1:wth_last) => o2_nc(wth_1:wth_last)
    ELSE
      o2(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'o3_mmr_radinput',i_time,o3_nc)
    IF (ALLOCATED(o3_nc)) THEN
      o3(wth_1:wth_last) => o3_nc(wth_1:wth_last)
    ELSE
      o3(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'rb_mmr_radinput',i_time,rb_nc)
    IF (ALLOCATED(rb_nc)) THEN
      rb(wth_1:wth_last) => rb_nc(wth_1:wth_last)
    ELSE
      rb(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF 
    CALL read_nc(input_file,'so2_mmr_radinput',i_time,so2_nc)
    IF (ALLOCATED(so2_nc)) THEN
      so2(wth_1:wth_last) => so2_nc(wth_1:wth_last)
    ELSE
      so2(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'tio_mmr_radinput',i_time,tio_nc)
    IF (ALLOCATED(tio_nc)) THEN
      tio(wth_1:wth_last) => tio_nc(wth_1:wth_last)
    ELSE
      tio(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF 

    CALL read_nc(input_file,'vo_mmr_radinput',i_time,vo_nc)
    IF (ALLOCATED(vo_nc)) THEN
      vo(wth_1:wth_last) => vo_nc(wth_1:wth_last)
    ELSE
      vo(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'tile_fraction_radinput',i_time,tile_fraction_nc)
    IF (ALLOCATED(tile_fraction_nc)) THEN
      tile_fraction(tile_1:tile_last) => tile_fraction_nc(tile_1:tile_last)
    ELSE
      tile_fraction(tile_1:tile_last) => dummy_nc(tile_1:tile_last)
    END IF

    CALL read_nc(input_file,'tile_temperature_radinput',i_time,tile_temperature_nc)
    IF (ALLOCATED(tile_temperature_nc)) THEN
      tile_temperature(tile_1:tile_last) => tile_temperature_nc(tile_1:tile_last)
    ELSE
      tile_temperature(tile_1:tile_last) => dummy_nc(tile_1:tile_last)
    END IF
    
    CALL read_nc(input_file,'tile_lw_albedo_radinput',i_time,tile_lw_albedo_nc)
    IF (ALLOCATED(tile_lw_albedo_nc)) THEN
      tile_lw_albedo(rtile_lw_1:rtile_lw_last) => tile_lw_albedo_nc(rtile_lw_1:rtile_lw_last)
    ELSE
      tile_lw_albedo(rtile_lw_1:rtile_lw_last) => dummy_nc(rtile_lw_1:rtile_lw_last)
    END IF

    CALL read_nc(input_file,'tile_sw_diffuse_albedo_radinput',i_time,tile_sw_diffuse_albedo_nc)
    IF (ALLOCATED(tile_sw_diffuse_albedo_nc)) THEN
      tile_sw_diffuse_albedo(rtile_sw_1:rtile_sw_last) => tile_sw_diffuse_albedo_nc(rtile_sw_1:rtile_sw_last)
    ELSE
      tile_sw_diffuse_albedo(rtile_sw_1:rtile_sw_last) => dummy_nc(rtile_sw_1:rtile_sw_last)
    END IF

    CALL read_nc(input_file,'tile_sw_direct_albedo_radinput',i_time,tile_sw_direct_albedo_nc)
    IF (ALLOCATED(tile_sw_direct_albedo_nc)) THEN
      tile_sw_direct_albedo(rtile_sw_1:rtile_sw_last) => tile_sw_direct_albedo_nc(rtile_sw_1:rtile_sw_last)
    ELSE
      tile_sw_direct_albedo(rtile_sw_1:rtile_sw_last) => dummy_nc(rtile_sw_1:rtile_sw_last)
    END IF

    CALL read_nc(input_file,'radiative_cloud_fraction_radinput',i_time,radiative_cloud_fraction_nc)
    IF (ALLOCATED(radiative_cloud_fraction_nc)) THEN
      radiative_cloud_fraction(wth_1:wth_last) => radiative_cloud_fraction_nc(wth_1:wth_last)
    ELSE
      radiative_cloud_fraction(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'liquid_fraction_radinput',i_time,liquid_fraction_nc)
    IF (ALLOCATED(liquid_fraction_nc)) THEN
      liquid_fraction(wth_1:wth_last) => liquid_fraction_nc(wth_1:wth_last)
    ELSE
      liquid_fraction(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'frozen_fraction_radinput',i_time,frozen_fraction_nc)
    IF (ALLOCATED(frozen_fraction_nc)) THEN
      frozen_fraction(wth_1:wth_last) => frozen_fraction_nc(wth_1:wth_last)
    ELSE
      frozen_fraction(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'liquid_mmr_radinput',i_time,mcl_nc)
    IF (ALLOCATED(mcl_nc)) THEN
      mcl(wth_1:wth_last) => mcl_nc(wth_1:wth_last)
    ELSE
      mcl(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'frozen_mmr_radinput',i_time,mci_nc)
    IF (ALLOCATED(mci_nc)) THEN
      mci(wth_1:wth_last) => mci_nc(wth_1:wth_last)
    ELSE
      mci(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'ice_no_conc_radinput',i_time,n_ice_nc)
    IF (ALLOCATED(n_ice_nc)) THEN
      n_ice(wth_1:wth_last) => n_ice_nc(wth_1:wth_last)
    ELSE
      n_ice(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'cloud_drop_no_conc_radinput',i_time,cloud_drop_no_conc_nc)
    IF (ALLOCATED(cloud_drop_no_conc_nc)) THEN
      cloud_drop_no_conc(wth_1:wth_last) => cloud_drop_no_conc_nc(wth_1:wth_last)
    ELSE
      cloud_drop_no_conc(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'radiative_conv_fraction_radinput',i_time,radiative_conv_fraction_nc)
    IF (ALLOCATED(radiative_conv_fraction_nc)) THEN
      radiative_conv_fraction(wth_1:wth_last) => radiative_conv_fraction_nc(wth_1:wth_last)
    ELSE
      radiative_conv_fraction(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'conv_liquid_fraction_radinput',i_time,conv_liquid_fraction_nc)
    IF (ALLOCATED(conv_liquid_fraction_nc)) THEN
      conv_liquid_fraction(wth_1:wth_last) => conv_liquid_fraction_nc(wth_1:wth_last)
    ELSE
      conv_liquid_fraction(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'conv_frozen_fraction_radinput',i_time,conv_frozen_fraction_nc)
    IF (ALLOCATED(conv_frozen_fraction_nc)) THEN
      conv_frozen_fraction(wth_1:wth_last) => conv_frozen_fraction_nc(wth_1:wth_last)
    ELSE
      conv_frozen_fraction(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'conv_liquid_mmr_radinput',i_time,conv_liquid_mmr_nc)
    IF (ALLOCATED(conv_liquid_mmr_nc)) THEN
      conv_liquid_mmr(wth_1:wth_last) => conv_liquid_mmr_nc(wth_1:wth_last)
    ELSE
      conv_liquid_mmr(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'conv_frozen_mmr_radinput',i_time,conv_frozen_mmr_nc)
    IF (ALLOCATED(conv_frozen_mmr_nc)) THEN
      conv_frozen_mmr(wth_1:wth_last) => conv_frozen_mmr_nc(wth_1:wth_last)
    ELSE
      conv_frozen_mmr(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'sigma_condensate_radinput',i_time,sigma_mc_nc)
    IF (ALLOCATED(sigma_mc_nc)) THEN
      sigma_mc(wth_1:wth_last) => sigma_mc_nc(wth_1:wth_last)
    ELSE
      sigma_mc(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF
    
    CALL read_nc(input_file,'layer_heat_capacity_radinput',i_time,layer_heat_capacity_nc)
    IF (ALLOCATED(layer_heat_capacity_nc)) THEN
      layer_heat_capacity(wth_1:wth_last) => layer_heat_capacity_nc(wth_1:wth_last)
    ELSE
      layer_heat_capacity(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'sulphuric_aerosol_radinput',i_time,sulphuric_nc)
    IF (ALLOCATED(sulphuric_nc)) THEN
      sulphuric(wth_1:wth_last) => sulphuric_nc(wth_1:wth_last)
    ELSE
      sulphuric(wth_1:wth_last) => dummy_nc(wth_1:wth_last)
    END IF

    CALL read_nc(input_file,'cos_zenith_angle_rts',i_time,cos_zenith_angle_rts_nc)
    IF (ALLOCATED(cos_zenith_angle_rts_nc)) THEN
      cos_zenith_angle_rts(twod_1:twod_last) => cos_zenith_angle_rts_nc(twod_1:twod_last)
    ELSE
      cos_zenith_angle_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF

    CALL read_nc(input_file,'lit_fraction_rts',i_time,lit_fraction_rts_nc)
    IF (ALLOCATED(cos_zenith_angle_rts_nc)) THEN
      lit_fraction_rts(twod_1:twod_last) => lit_fraction_rts_nc(twod_1:twod_last)
    ELSE
      lit_fraction_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF

    CALL read_nc(input_file,'stellar_irradiance_rts',i_time,stellar_irradiance_rts_nc)
    IF (ALLOCATED(stellar_irradiance_rts_nc)) THEN
      stellar_irradiance_rts(twod_1:twod_last) => stellar_irradiance_rts_nc(twod_1:twod_last)
    ELSE
      stellar_irradiance_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF

    CALL read_nc(input_file,'orographic_correction_rts',i_time,orographic_correction_rts_nc)
    IF (ALLOCATED(orographic_correction_rts_nc)) THEN
      orographic_correction_rts(twod_1:twod_last) => orographic_correction_rts_nc(twod_1:twod_last)
    ELSE
      orographic_correction_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF
    
    CALL read_nc(input_file,'rand_seed_real_radinput',i_time,rand_seed_real_radinput_nc)
    IF (ALLOCATED(rand_seed_real_radinput_nc)) THEN
      rand_seed_real_radinput(twod_1:twod_last) => rand_seed_real_radinput_nc(twod_1:twod_last)
    ELSE
      rand_seed_real_radinput(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF

    CALL read_nc(input_file,'aer_mix_ratio',i_time,aer_mix_ratio_nc)
    IF (ALLOCATED(aer_mix_ratio_nc)) THEN
      aer_mix_ratio(mode_1:mode_last) => aer_mix_ratio_nc(mode_1:mode_last)
    ELSE
      aer_mix_ratio(mode_1:mode_last) => dummy_nc(mode_1:mode_last)
    END IF

    CALL read_nc(input_file,'aer_sw_absorption',i_time,aer_sw_absorption_nc)
    IF (ALLOCATED(aer_sw_absorption_nc)) THEN
      aer_sw_absorption(rmode_sw_1:rmode_sw_last) => aer_sw_absorption_nc(rmode_sw_1:rmode_sw_last)
    ELSE
      aer_sw_absorption(rmode_sw_1:rmode_sw_last) => dummy_nc(rmode_sw_1:rmode_sw_last)
    END IF

    CALL read_nc(input_file,'aer_sw_scattering',i_time,aer_sw_scattering_nc)
    IF (ALLOCATED(aer_sw_scattering_nc)) THEN
      aer_sw_scattering(rmode_sw_1:rmode_sw_last) => aer_sw_scattering_nc(rmode_sw_1:rmode_sw_last)
    ELSE
      aer_sw_scattering(rmode_sw_1:rmode_sw_last) => dummy_nc(rmode_sw_1:rmode_sw_last)
    END IF
    
    CALL read_nc(input_file,'aer_sw_asymmetry',i_time,aer_sw_asymmetry_nc)
    IF (ALLOCATED(aer_sw_asymmetry_nc)) THEN
      aer_sw_asymmetry(rmode_sw_1:rmode_sw_last) => aer_sw_asymmetry_nc(rmode_sw_1:rmode_sw_last)
    ELSE
      aer_sw_asymmetry(rmode_sw_1:rmode_sw_last) => dummy_nc(rmode_sw_1:rmode_sw_last)
    END IF

    CALL read_nc(input_file,'aer_lw_absorption',i_time,aer_lw_absorption_nc)
    IF (ALLOCATED(aer_lw_absorption_nc)) THEN
      aer_lw_absorption(rmode_lw_1:rmode_lw_last) => aer_lw_absorption_nc(rmode_lw_1:rmode_lw_last)
    ELSE
      aer_lw_absorption(rmode_lw_1:rmode_lw_last) => dummy_nc(rmode_lw_1:rmode_lw_last)
    END IF

    CALL read_nc(input_file,'aer_lw_scattering',i_time,aer_lw_scattering_nc)
    IF (ALLOCATED(aer_lw_scattering_nc)) THEN
      aer_lw_scattering(rmode_lw_1:rmode_lw_last) => aer_lw_scattering_nc(rmode_lw_1:rmode_lw_last)
    ELSE
      aer_lw_scattering(rmode_lw_1:rmode_lw_last) => dummy_nc(rmode_lw_1:rmode_lw_last)
    END IF
    
    CALL read_nc(input_file,'aer_lw_absorption',i_time,aer_lw_asymmetry_nc)
    IF (ALLOCATED(aer_lw_asymmetry_nc)) THEN
      aer_lw_asymmetry(rmode_lw_1:rmode_lw_last) => aer_lw_asymmetry_nc(rmode_lw_1:rmode_lw_last)
    ELSE
      aer_lw_asymmetry(rmode_lw_1:rmode_lw_last) => dummy_nc(rmode_lw_1:rmode_lw_last)
    END IF

    CALL read_nc(input_file,'sw_up_rts',i_time,sw_up_rts_nc)
    IF (ALLOCATED(sw_up_rts_nc)) THEN
      sw_up_rts(flux_0:flux_last) => sw_up_rts_nc(flux_0:flux_last)
    ELSE
      sw_up_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF

    CALL read_nc(input_file,'sw_down_rts',i_time,sw_down_rts_nc)
    IF (ALLOCATED(sw_down_rts_nc)) THEN
      sw_down_rts(flux_0:flux_last) => sw_down_rts_nc(flux_0:flux_last)
    ELSE
      sw_down_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF

    CALL read_nc(input_file,'sw_heating_rate_rts',i_time,sw_heating_rate_rts_nc)
    IF (ALLOCATED(sw_heating_rate_rts_nc)) THEN
      sw_heating_rate_rts(flux_0:flux_last) => sw_heating_rate_rts_nc(flux_0:flux_last)
    ELSE
      sw_heating_rate_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'sw_up_clear_rts',i_time,sw_up_clear_rts_nc)
    IF (ALLOCATED(sw_up_clear_rts_nc)) THEN
      sw_up_clear_rts(flux_0:flux_last) => sw_up_clear_rts_nc(flux_0:flux_last)
    ELSE
      sw_up_clear_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'sw_down_clear_rts',i_time,sw_down_clear_rts_nc)
    IF (ALLOCATED(sw_down_clear_rts_nc)) THEN
      sw_down_clear_rts(flux_0:flux_last) => sw_down_clear_rts_nc(flux_0:flux_last)
    ELSE
      sw_down_clear_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'sw_down_clear_clean_surf_rts',i_time,sw_down_clear_clean_surf_rts_nc)
    IF (ALLOCATED(sw_down_clear_clean_surf_rts_nc)) THEN
      sw_down_clear_clean_surf_rts(twod_1:twod_last) => sw_down_clear_clean_surf_rts_nc(twod_1:twod_last)
    ELSE
      sw_down_clear_clean_surf_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF
    
    CALL read_nc(input_file,'sw_up_clear_clean_surf_rts',i_time,sw_up_clear_clean_surf_rts_nc)
    IF (ALLOCATED(sw_up_clear_clean_surf_rts_nc)) THEN
      sw_up_clear_clean_surf_rts(twod_1:twod_last) => sw_up_clear_clean_surf_rts_nc(twod_1:twod_last)
    ELSE
      sw_up_clear_clean_surf_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF
    
    CALL read_nc(input_file,'sw_up_clear_clean_toa_rts',i_time,sw_up_clear_clean_toa_rts_nc)
    IF (ALLOCATED(sw_up_clear_clean_toa_rts_nc)) THEN
      sw_up_clear_clean_toa_rts(twod_1:twod_last) => sw_up_clear_clean_toa_rts_nc(twod_1:twod_last)
    ELSE
      sw_up_clear_clean_toa_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF

    CALL read_nc(input_file,'lw_up_rts',i_time,lw_up_rts_nc)
    IF (ALLOCATED(lw_up_rts_nc)) THEN
      lw_up_rts(flux_0:flux_last) => lw_up_rts_nc(flux_0:flux_last)
    ELSE
      lw_up_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF

    CALL read_nc(input_file,'lw_down_rts',i_time,lw_down_rts_nc)
    IF (ALLOCATED(lw_down_rts_nc)) THEN
      lw_down_rts(flux_0:flux_last) => lw_down_rts_nc(flux_0:flux_last)
    ELSE
      lw_down_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF

    CALL read_nc(input_file,'lw_heating_rate_rts',i_time,lw_heating_rate_rts_nc)
    IF (ALLOCATED(lw_heating_rate_rts_nc)) THEN
      lw_heating_rate_rts(flux_0:flux_last) => lw_heating_rate_rts_nc(flux_0:flux_last)
    ELSE
      lw_heating_rate_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'lw_up_clear_rts',i_time,lw_up_clear_rts_nc)
    IF (ALLOCATED(lw_up_clear_rts_nc)) THEN
      lw_up_clear_rts(flux_0:flux_last) => lw_up_clear_rts_nc(flux_0:flux_last)
    ELSE
      lw_up_clear_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'lw_down_clear_rts',i_time,lw_down_clear_rts_nc)
    IF (ALLOCATED(lw_down_clear_rts_nc)) THEN
      lw_down_clear_rts(flux_0:flux_last) => lw_down_clear_rts_nc(flux_0:flux_last)
    ELSE
      lw_down_clear_rts(flux_0:flux_last) => dummy_nc(flux_0:flux_last)
    END IF
    
    CALL read_nc(input_file,'lw_down_clear_clean_surf_rts',i_time,lw_down_clear_clean_surf_rts_nc)
    IF (ALLOCATED(lw_down_clear_clean_surf_rts_nc)) THEN
      lw_down_clear_clean_surf_rts(twod_1:twod_last) => lw_down_clear_clean_surf_rts_nc(twod_1:twod_last)
    ELSE
      lw_down_clear_clean_surf_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF
    
    CALL read_nc(input_file,'lw_up_clear_clean_surf_rts',i_time,lw_up_clear_clean_surf_rts_nc)
    IF (ALLOCATED(lw_up_clear_clean_surf_rts_nc)) THEN
      lw_up_clear_clean_surf_rts(twod_1:twod_last) => lw_up_clear_clean_surf_rts_nc(twod_1:twod_last)
    ELSE
      lw_up_clear_clean_surf_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF
    
    CALL read_nc(input_file,'lw_up_clear_clean_toa_rts',i_time,lw_up_clear_clean_toa_rts_nc)
    IF (ALLOCATED(lw_up_clear_clean_toa_rts_nc)) THEN
      lw_up_clear_clean_toa_rts(twod_1:twod_last) => lw_up_clear_clean_toa_rts_nc(twod_1:twod_last)
    ELSE
      lw_up_clear_clean_toa_rts(twod_1:twod_last) => dummy_nc(twod_1:twod_last)
    END IF

    rand_seed=NINT(rand_seed_real_radinput,RealExt)


    ! Set the number of cloud layers to the highest cloud layer

    n_cloud_layer(1:n_profile) = 0
    DO j = 1, n_profile
      DO i = nlayers, 1, -1
        k = (j-1)*(nlayers+1)+i+1
        IF ( (radiative_cloud_fraction(k) + radiative_conv_fraction(k)) &
           > cloud_fraction_min ) then
           n_cloud_layer(j) = i
           EXIT
        END IF
      END DO ! nlayers
    END DO ! n_profile


    ! SW call

    sw_diag%heating_rate => sw_heating_rate
    sw_diag%flux_up => sw_flux_up
    sw_diag%flux_down => sw_flux_down
    sw_diag%flux_up_clear => sw_flux_up_clear
    sw_diag%flux_down_clear => sw_flux_down_clear

    sw_diag%flux_down_clear_clean_surf => sw_flux_down_clear_clean_surf
    sw_diag%flux_up_clear_clean_surf => sw_flux_up_clear_clean_surf
    sw_diag%flux_up_clear_clean_toa => sw_flux_up_clear_clean_toa

    DO k = 0, nlayers
      profile_list = pack( [(l, l=1, n_profile)], &
                           lit_fraction_rts(twod_1:twod_last) > 0.0_RealExt &
                           .AND. n_cloud_layer(twod_1:twod_last) == k )
      n_profile_list = size(profile_list)
      IF (n_profile_list > 0) THEN
      CALL runes(n_profile_list, nlayers, sw_diag,                            &
        spectrum_name          = 'sw',                                        &
        i_source               = ip_source_illuminate,                        &
        profile_list           = profile_list,                                &
        n_layer_stride         = nlayers+1,                                   &
        n_cloud_layer          = k,                                           &
        p_layer_1d             = pressure_in_wth(wth_1:wth_last),             &
        t_layer_1d             = temperature_in_wth(wth_1:wth_last),          &
        mass_1d                = d_mass(wth_1:wth_last),                      &
        density_1d             = rho_in_wth(wth_1:wth_last),                  &
        ch4_1d                 = ch4(wth_1:wth_last),                         &
        co_1d                  = co(wth_1:wth_last),                          &
        co2_1d                 = co2(wth_1:wth_last),                         &
        cs_1d                  = cs(wth_1:wth_last),                          &
        h2_1d                  = h2(wth_1:wth_last),                          &
        h2o_1d                 = h2o(wth_1:wth_last),                         &
        hcn_1d                 = hcn(wth_1:wth_last),                         &
        he_1d                  = he(wth_1:wth_last),                          &
        k_1d                   = potassium(wth_1:wth_last),                   &
        li_1d                  = li(wth_1:wth_last),                          &
        n2_1d                  = n2(wth_1:wth_last),                          &
        n2o_1d                 = n2o(wth_1:wth_last),                         &
        na_1d                  = na(wth_1:wth_last),                          &
        nh3_1d                 = nh3(wth_1:wth_last),                         &
        o2_1d                  = o2(wth_1:wth_last),                          &
        o3_1d                  = o3(wth_1:wth_last),                          &
        rb_1d                  = rb(wth_1:wth_last),                          &
        so2_1d                 = so2(wth_1:wth_last),                         &
        tio_1d                 = tio(wth_1:wth_last),                         &
        vo_1d                  = vo(wth_1:wth_last),                          &
        cfc11_mix_ratio        = cfc11_mix_ratio,                             &
        cfc113_mix_ratio       = cfc113_mix_ratio,                            &
        cfc12_mix_ratio        = cfc12_mix_ratio,                             &
        ch4_mix_ratio          = ch4_mix_ratio,                               &
        co_mix_ratio           = co_mix_ratio,                                &
        co2_mix_ratio          = co2_mix_ratio,                               &
        cs_mix_ratio           = cs_mix_ratio,                                &
        h2_mix_ratio           = h2_mix_ratio,                                &
        h2o_mix_ratio          = h2o_mix_ratio,                               &
        hcfc22_mix_ratio       = hcfc22_mix_ratio,                            &
        hcn_mix_ratio          = hcn_mix_ratio,                               &
        he_mix_ratio           = he_mix_ratio,                                &
        hfc134a_mix_ratio      = hfc134a_mix_ratio,                           &
        k_mix_ratio            = k_mix_ratio,                                 &
        li_mix_ratio           = li_mix_ratio,                                &
        n2_mix_ratio           = n2_mix_ratio,                                &
        n2o_mix_ratio          = n2o_mix_ratio,                               &
        na_mix_ratio           = na_mix_ratio,                                &
        nh3_mix_ratio          = nh3_mix_ratio,                               &
        o2_mix_ratio           = o2_mix_ratio,                                &
        o3_mix_ratio           = o3_mix_ratio,                                &
        rb_mix_ratio           = rb_mix_ratio,                                &
        so2_mix_ratio          = so2_mix_ratio,                               &
        tio_mix_ratio          = tio_mix_ratio,                               &
        vo_mix_ratio           = vo_mix_ratio,                                &
        l_ch4_well_mixed       = ch4_well_mixed,                              &
        l_co_well_mixed        = co_well_mixed,                               &
        l_co2_well_mixed       = co2_well_mixed,                              &
        l_cs_well_mixed        = cs_well_mixed,                               &
        l_h2_well_mixed        = h2_well_mixed,                               &
        l_h2o_well_mixed       = h2o_well_mixed,                              &
        l_hcn_well_mixed       = hcn_well_mixed,                              &
        l_he_well_mixed        = he_well_mixed,                               &
        l_k_well_mixed         = k_well_mixed,                                &
        l_li_well_mixed        = li_well_mixed,                               &
        l_n2_well_mixed        = n2_well_mixed,                               &
        l_n2o_well_mixed       = n2o_well_mixed,                              &
        l_na_well_mixed        = na_well_mixed,                               &
        l_nh3_well_mixed       = nh3_well_mixed,                              &
        l_o2_well_mixed        = o2_well_mixed,                               &
        l_o3_well_mixed        = o3_well_mixed,                               &
        l_rb_well_mixed        = rb_well_mixed,                               &
        l_so2_well_mixed       = so2_well_mixed,                              &
        l_tio_well_mixed       = tio_well_mixed,                              &
        l_vo_well_mixed        = vo_well_mixed,                               &
        cos_zenith_angle       = cos_zenith_angle_rts(twod_1:twod_last),      &
        solar_irrad            = stellar_irradiance_rts(twod_1:twod_last),    &
        l_orog                 = l_orog,                                      &
        orog_corr              = orographic_correction_rts(twod_1:twod_last), &
        n_tile                 = n_surf_tile,                                 &
        frac_tile_1d           = tile_fraction(tile_1:tile_last),             &
        albedo_diff_tile_1d    =                                              &
          tile_sw_diffuse_albedo(rtile_sw_1:rtile_sw_last),                   &
        albedo_dir_tile_1d     =                                              &
          tile_sw_direct_albedo(rtile_sw_1:rtile_sw_last),                    &
        cloud_frac_1d          = radiative_cloud_fraction(wth_1:wth_last),    &
        liq_frac_1d            = liquid_fraction(wth_1:wth_last),             &
        ice_frac_1d            = frozen_fraction(wth_1:wth_last),             &
        liq_mmr_1d             = mcl(wth_1:wth_last),                         &
        ice_mmr_1d             = mci(wth_1:wth_last),                         &
        ice_nc_1d              = n_ice(wth_1:wth_last),                       &
        ice_conv_nc_1d         = n_ice(wth_1:wth_last),                       &
        liq_dim_constant       = constant_droplet_effective_radius,           &
        liq_nc_1d              = cloud_drop_no_conc(wth_1:wth_last),          &
        conv_frac_1d           = radiative_conv_fraction(wth_1:wth_last),     &
        liq_conv_frac_1d       = conv_liquid_fraction(wth_1:wth_last),        &
        ice_conv_frac_1d       = conv_frozen_fraction(wth_1:wth_last),        &
        liq_conv_mmr_1d        = conv_liquid_mmr(wth_1:wth_last),             &
        ice_conv_mmr_1d        = conv_frozen_mmr(wth_1:wth_last),             &
        liq_conv_dim_constant  = constant_droplet_effective_radius,           &
        liq_conv_nc_1d         = cloud_drop_no_conc(wth_1:wth_last),          &
        liq_rsd_1d             = sigma_mc(wth_1:wth_last),                    &
        ice_rsd_1d             = sigma_mc(wth_1:wth_last),                    &
        cloud_vertical_decorr  = cloud_vertical_decorr,                       &
        conv_vertical_decorr   = cloud_vertical_decorr,                       &
        liq_dim_aparam         = liu_aparam,                                  &
        liq_dim_bparam         = liu_bparam,                                  &
        rand_seed              = rand_seed(twod_1:twod_last),                 &
        layer_heat_capacity_1d = layer_heat_capacity(wth_1:wth_last),         &
        l_rayleigh             = l_rayleigh_sw,                               &
        l_mixing_ratio         = .true.,                                      &
        i_cloud_representation = i_cloud_representation,                      &
        i_overlap              = i_overlap,                                   &
        i_inhom                = i_inhom,                                     &
        i_cloud_entrapment     = i_cloud_entrapment,                          &
        i_drop_re              = i_drop_re,                                   &
        i_st_water             = i_cloud_liq_type_sw,                         &
        i_st_ice               = i_cloud_ice_type_sw,                         &
        i_cnv_water            = i_cloud_liq_type_sw,                         &
        i_cnv_ice              = i_cloud_ice_type_sw,                         &
        l_sulphuric            = sulphuric_strat_climatology,                 &
        sulphuric_1d           = sulphuric(wth_1:wth_last),                   &
        l_aerosol_mode         = l_aerosol_mode_sw,                           &
        n_aer_mode             = n_aer_mode_sw,                               &
        aer_mix_ratio_1d       = aer_mix_ratio(mode_1:mode_last),             &
        aer_absorption_1d      = aer_sw_absorption(rmode_sw_1:rmode_sw_last), &
        aer_scattering_1d      = aer_sw_scattering(rmode_sw_1:rmode_sw_last), &
        aer_asymmetry_1d       = aer_sw_asymmetry(rmode_sw_1:rmode_sw_last),  &
        l_invert               = .true.,                                      &
        l_profile_last         = .true.)
      END IF
    END DO

    IF (l_ascii) THEN

      WRITE(iu_out_sw_lfric,'(A, I0)') 'time: ', i_time
      DO j=1, n_profile
        IF (lit_fraction_rts(j) > 0.0_RealExt) THEN
          WRITE(iu_out_sw_lfric,'(A, I0)') 'profile: ', j
          WRITE(iu_out_sw_lfric,'(A,A)') &
            'SW: fluxes (Wm-2): down | up | down clear | up clear', &
            ' | heating rate (K/day)'
          DO i=0, nlayers
            WRITE(iu_out_sw_lfric, '(I4, 2X, 5F11.5)') i, &
              sw_down_rts((j-1)*(nlayers+1)+i+1), &
              sw_up_rts((j-1)*(nlayers+1)+i+1), &              
              sw_down_clear_rts((j-1)*(nlayers+1)+i+1), &
              sw_up_clear_rts((j-1)*(nlayers+1)+i+1), &
              sw_heating_rate_rts((j-1)*(nlayers+1)+i+1)*seconds_per_day
          END DO ! nlayers
          WRITE(iu_out_sw_lfric,'(A)') &
            'SW: clear-clean fluxes (Wm-2): down surf | up surf | up toa'
          WRITE(iu_out_sw_lfric, '(3F11.5,/)') &
            sw_down_clear_clean_surf_rts(j), &
            sw_up_clear_clean_surf_rts(j), &
            sw_up_clear_clean_toa_rts(j)
        END IF ! lit_fraction_rts > 0.0_RealExt
      END DO ! n_profile

      WRITE(iu_out_sw_offline,'(A, I0)') 'time: ', i_time
      DO j=1, n_profile
        IF (lit_fraction_rts(j) > 0.0_RealExt) THEN
          WRITE(iu_out_sw_offline,'(A, I0)') 'profile: ', j
          WRITE(iu_out_sw_offline,'(A,A)') &
            'SW: fluxes (Wm-2): down | up | down clear | up clear', &
            ' | heating rate (K/day)'
          DO i=0, nlayers
            WRITE(iu_out_sw_offline, '(I4, 2X, 5F11.5)') i, &
              sw_flux_down(i, j), &
              sw_flux_up(i, j), &
              sw_flux_down_clear(i, j), &
              sw_flux_up_clear(i, j), &
              sw_heating_rate(i, j)*seconds_per_day
          END DO ! nlayers
          WRITE(iu_out_sw_offline,'(A)') &
            'SW: clear-clean fluxes (Wm-2): down surf | up surf | up toa'
          WRITE(iu_out_sw_offline, '(3F11.5,/)') &
            sw_flux_down_clear_clean_surf(j), &
            sw_flux_up_clear_clean_surf(j), &
            sw_flux_up_clear_clean_toa(j)
        END IF ! lit_fraction_rts > 0.0_RealExt
      END DO ! n_profile
    
    END IF ! l_ascii
    

    ! LW call

    lw_diag%heating_rate => lw_heating_rate
    lw_diag%flux_up => lw_flux_up
    lw_diag%flux_down => lw_flux_down
    lw_diag%flux_up_clear => lw_flux_up_clear
    lw_diag%flux_down_clear => lw_flux_down_clear

    lw_diag%flux_down_clear_clean_surf => lw_flux_down_clear_clean_surf
    lw_diag%flux_up_clear_clean_surf => lw_flux_up_clear_clean_surf
    lw_diag%flux_up_clear_clean_toa => lw_flux_up_clear_clean_toa

    DO k=0, nlayers
      profile_list = pack( [(l, l=1, n_profile)], &
                           n_cloud_layer(twod_1:twod_last) == k )
      n_profile_list = size(profile_list)
      IF (n_profile_list > 0) THEN
      CALL runes(n_profile_list, nlayers, lw_diag,                            &
        spectrum_name          = 'lw',                                        &
        i_source               = ip_source_thermal,                           &
        profile_list           = profile_list,                                &
        n_layer_stride         = nlayers+1,                                   &
        n_cloud_layer          = k,                                           &
        p_layer_1d             = pressure_in_wth(wth_1:wth_last),             &
        t_layer_1d             = temperature_in_wth(wth_1:wth_last),          &
        mass_1d                = d_mass(wth_1:wth_last),                      &
        density_1d             = rho_in_wth(wth_1:wth_last),                  &
        t_level_1d             = t_layer_boundaries(flux_0:flux_last),        &
        ch4_1d                 = ch4(wth_1:wth_last),                         &
        co_1d                  = co(wth_1:wth_last),                          &
        co2_1d                 = co2(wth_1:wth_last),                         &
        cs_1d                  = cs(wth_1:wth_last),                          &
        h2_1d                  = h2(wth_1:wth_last),                          &
        h2o_1d                 = h2o(wth_1:wth_last),                         &
        hcn_1d                 = hcn(wth_1:wth_last),                         &
        he_1d                  = he(wth_1:wth_last),                          &
        k_1d                   = potassium(wth_1:wth_last),                   &
        li_1d                  = li(wth_1:wth_last),                          &
        n2_1d                  = n2(wth_1:wth_last),                          &
        n2o_1d                 = n2o(wth_1:wth_last),                         &
        na_1d                  = na(wth_1:wth_last),                          &
        nh3_1d                 = nh3(wth_1:wth_last),                         &
        o2_1d                  = o2(wth_1:wth_last),                          &
        o3_1d                  = o3(wth_1:wth_last),                          &
        rb_1d                  = rb(wth_1:wth_last),                          &
        so2_1d                 = so2(wth_1:wth_last),                         &
        tio_1d                 = tio(wth_1:wth_last),                         &
        vo_1d                  = vo(wth_1:wth_last),                          &
        cfc11_mix_ratio        = cfc11_mix_ratio,                             &
        cfc113_mix_ratio       = cfc113_mix_ratio,                            &
        cfc12_mix_ratio        = cfc12_mix_ratio,                             &
        ch4_mix_ratio          = ch4_mix_ratio,                               &
        co_mix_ratio           = co_mix_ratio,                                &
        co2_mix_ratio          = co2_mix_ratio,                               &
        cs_mix_ratio           = cs_mix_ratio,                                &
        h2_mix_ratio           = h2_mix_ratio,                                &
        h2o_mix_ratio          = h2o_mix_ratio,                               &
        hcfc22_mix_ratio       = hcfc22_mix_ratio,                            &
        hcn_mix_ratio          = hcn_mix_ratio,                               &
        he_mix_ratio           = he_mix_ratio,                                &
        hfc134a_mix_ratio      = hfc134a_mix_ratio,                           &
        k_mix_ratio            = k_mix_ratio,                                 &
        li_mix_ratio           = li_mix_ratio,                                &
        n2_mix_ratio           = n2_mix_ratio,                                &
        n2o_mix_ratio          = n2o_mix_ratio,                               &
        na_mix_ratio           = na_mix_ratio,                                &
        nh3_mix_ratio          = nh3_mix_ratio,                               &
        o2_mix_ratio           = o2_mix_ratio,                                &
        o3_mix_ratio           = o3_mix_ratio,                                &
        rb_mix_ratio           = rb_mix_ratio,                                &
        so2_mix_ratio          = so2_mix_ratio,                               &
        tio_mix_ratio          = tio_mix_ratio,                               &
        vo_mix_ratio           = vo_mix_ratio,                                &
        l_ch4_well_mixed       = ch4_well_mixed,                              &
        l_co_well_mixed        = co_well_mixed,                               &
        l_co2_well_mixed       = co2_well_mixed,                              &
        l_cs_well_mixed        = cs_well_mixed,                               &
        l_h2_well_mixed        = h2_well_mixed,                               &
        l_h2o_well_mixed       = h2o_well_mixed,                              &
        l_hcn_well_mixed       = hcn_well_mixed,                              &
        l_he_well_mixed        = he_well_mixed,                               &
        l_k_well_mixed         = k_well_mixed,                                &
        l_li_well_mixed        = li_well_mixed,                               &
        l_n2_well_mixed        = n2_well_mixed,                               &
        l_n2o_well_mixed       = n2o_well_mixed,                              &
        l_na_well_mixed        = na_well_mixed,                               &
        l_nh3_well_mixed       = nh3_well_mixed,                              &
        l_o2_well_mixed        = o2_well_mixed,                               &
        l_o3_well_mixed        = o3_well_mixed,                               &
        l_rb_well_mixed        = rb_well_mixed,                               &
        l_so2_well_mixed       = so2_well_mixed,                              &
        l_tio_well_mixed       = tio_well_mixed,                              &
        l_vo_well_mixed        = vo_well_mixed,                               &
        n_tile                 = n_surf_tile,                                 &
        frac_tile_1d           = tile_fraction(tile_1:tile_last),             &
        t_tile_1d              = tile_temperature(tile_1:tile_last),          &
        albedo_diff_tile_1d    = tile_lw_albedo(rtile_lw_1:rtile_lw_last),    &
        cloud_frac_1d          = radiative_cloud_fraction(wth_1:wth_last),    &
        liq_frac_1d            = liquid_fraction(wth_1:wth_last),             &
        ice_frac_1d            = frozen_fraction(wth_1:wth_last),             &
        liq_mmr_1d             = mcl(wth_1:wth_last),                         &
        ice_mmr_1d             = mci(wth_1:wth_last),                         &
        ice_nc_1d              = n_ice(wth_1:wth_last),                       &
        ice_conv_nc_1d         = n_ice(wth_1:wth_last),                       &
        liq_dim_constant       = constant_droplet_effective_radius,           &
        liq_nc_1d              = cloud_drop_no_conc(wth_1:wth_last),          &
        conv_frac_1d           = radiative_conv_fraction(wth_1:wth_last),     &
        liq_conv_frac_1d       = conv_liquid_fraction(wth_1:wth_last),        &
        ice_conv_frac_1d       = conv_frozen_fraction(wth_1:wth_last),        &
        liq_conv_mmr_1d        = conv_liquid_mmr(wth_1:wth_last),             &
        ice_conv_mmr_1d        = conv_frozen_mmr(wth_1:wth_last),             &
        liq_conv_dim_constant  = constant_droplet_effective_radius,           &
        liq_conv_nc_1d         = cloud_drop_no_conc(wth_1:wth_last),          &
        liq_rsd_1d             = sigma_mc(wth_1:wth_last),                    &
        ice_rsd_1d             = sigma_mc(wth_1:wth_last),                    &
        cloud_vertical_decorr  = cloud_vertical_decorr,                       &
        conv_vertical_decorr   = cloud_vertical_decorr,                       &
        liq_dim_aparam         = liu_aparam,                                  &
        liq_dim_bparam         = liu_bparam,                                  &
        rand_seed              = rand_seed(twod_1:twod_last),                 &
        layer_heat_capacity_1d = layer_heat_capacity(wth_1:wth_last),         &
        l_mixing_ratio         = .true.,                                      &
        i_scatter_method       = i_scatter_method_lw,                         &
        i_cloud_representation = i_cloud_representation,                      &
        i_overlap              = i_overlap,                                   &
        i_inhom                = i_inhom,                                     &
        i_cloud_entrapment     = i_cloud_entrapment,                          &
        i_drop_re              = i_drop_re,                                   &
        i_st_water             = i_cloud_liq_type_lw,                         &
        i_st_ice               = i_cloud_ice_type_lw,                         &
        i_cnv_water            = i_cloud_liq_type_lw,                         &
        i_cnv_ice              = i_cloud_ice_type_lw,                         &
        l_sulphuric            = sulphuric_strat_climatology,                 &
        sulphuric_1d           = sulphuric(wth_1:wth_last),                   &
        l_aerosol_mode         = l_aerosol_mode_lw,                           &
        n_aer_mode             = n_aer_mode_lw,                               &
        aer_mix_ratio_1d       = aer_mix_ratio(mode_1:mode_last),             &
        aer_absorption_1d      = aer_lw_absorption(rmode_lw_1:rmode_lw_last), &
        aer_scattering_1d      = aer_lw_scattering(rmode_lw_1:rmode_lw_last), &
        aer_asymmetry_1d       = aer_lw_asymmetry(rmode_lw_1:rmode_lw_last),  &
        l_invert               = .true.,                                      &
        l_profile_last         = .true.)
      END IF
    END DO


    IF (l_ascii) THEN

      WRITE(iu_out_lw_lfric,'(A, I0)') 'time: ', i_time
      DO j=1, n_profile
        WRITE(iu_out_lw_lfric,'(A, I0)') 'profile: ', j
        WRITE(iu_out_lw_lfric,'(A,A)') &
          'LW: fluxes (Wm-2): down | up | down clear | up clear', &
          ' | heating rate (K/day)'
        DO i=0, nlayers
          WRITE(iu_out_lw_lfric, '(I4, 2X, 5F11.5)') i, &
            lw_down_rts((j-1)*(nlayers+1)+i+1), &
            lw_up_rts((j-1)*(nlayers+1)+i+1), &
            lw_down_clear_rts((j-1)*(nlayers+1)+i+1), &
            lw_up_clear_rts((j-1)*(nlayers+1)+i+1), &
            lw_heating_rate_rts((j-1)*(nlayers+1)+i+1)*seconds_per_day
        END DO ! nlayers
        WRITE(iu_out_lw_lfric,'(A)') &
          'LW: clear-clean fluxes (Wm-2): down surf | up surf | up toa'
        WRITE(iu_out_lw_lfric, '(3F11.5,/)') &
          lw_down_clear_clean_surf_rts(j), &
          lw_up_clear_clean_surf_rts(j), &
          lw_up_clear_clean_toa_rts(j)
      END DO ! n_profile

      WRITE(iu_out_lw_offline,'(A, I0)') 'time: ', i_time
      DO j=1, n_profile
        WRITE(iu_out_lw_offline,'(A, I0)') 'profile: ', j
        WRITE(iu_out_lw_offline,'(A,A)') &
          'LW: fluxes (Wm-2): down | up | down clear | up clear', &
          ' | heating rate (K/day)'
        DO i=0, nlayers
          WRITE(iu_out_lw_offline, '(I4, 2X, 5F11.5)') i, &
            lw_flux_down(i, j), &
            lw_flux_up(i, j), &
            lw_flux_down_clear(i, j), &
            lw_flux_up_clear(i, j), &
            lw_heating_rate(i, j)*seconds_per_day
        END DO ! nlayers
        WRITE(iu_out_lw_offline,'(A)') &
          'LW: clear-clean fluxes (Wm-2): down surf | up surf | up toa'
        WRITE(iu_out_lw_offline, '(3F11.5,/)') &
          lw_flux_down_clear_clean_surf(j), &
          lw_flux_up_clear_clean_surf(j), &
          lw_flux_up_clear_clean_toa(j)
      END DO ! n_profile

    END IF ! l_ascii


    ! Add flux/heating rate fields calculated by the offline driver
    ! to the netCDF file
    ! (The write routine expects 1D fields)

    dummy_arr_out=pack(sw_heating_rate,.TRUE.)
    CALL write_nc(output_file,'sw_heating_rate_rts',i_time,dummy_arr_out)

    dummy_arr_out=pack(sw_flux_up,.TRUE.)
    CALL write_nc(output_file,'sw_up_rts',i_time,dummy_arr_out)

    dummy_arr_out=pack(sw_flux_down,.TRUE.)
    CALL write_nc(output_file,'sw_down_rts',i_time,dummy_arr_out)
  
    dummy_arr_out=pack(sw_flux_up_clear,.TRUE.)
    CALL write_nc(output_file,'sw_up_clear_rts',i_time,dummy_arr_out)

    dummy_arr_out=pack(sw_flux_down_clear,.TRUE.)
    CALL write_nc(output_file,'sw_down_clear_rts',i_time,dummy_arr_out)
    
    CALL write_nc(output_file,'sw_down_clear_clean_surf_rts',i_time,sw_flux_down_clear_clean_surf)
    CALL write_nc(output_file,'sw_up_clear_clean_surf_rts',i_time,sw_flux_up_clear_clean_surf)
    CALL write_nc(output_file,'sw_up_clear_clean_toa_rts',i_time,sw_flux_up_clear_clean_toa)

    dummy_arr_out=pack(lw_heating_rate,.TRUE.)
    CALL write_nc(output_file,'lw_heating_rate_rts',i_time,dummy_arr_out)
            
    dummy_arr_out=pack(lw_flux_up,.TRUE.)
    CALL write_nc(output_file,'lw_up_rts',i_time,dummy_arr_out)

    dummy_arr_out=pack(lw_flux_down,.TRUE.)
    CALL write_nc(output_file,'lw_down_rts',i_time,dummy_arr_out)

    dummy_arr_out=pack(lw_flux_up_clear,.TRUE.)
    CALL write_nc(output_file,'lw_up_clear_rts',i_time,dummy_arr_out)

    dummy_arr_out=pack(lw_flux_down_clear,.TRUE.)
    CALL write_nc(output_file,'lw_down_clear_rts',i_time,dummy_arr_out)

    CALL write_nc(output_file,'lw_down_clear_clean_surf_rts',i_time,lw_flux_down_clear_clean_surf)
    CALL write_nc(output_file,'lw_up_clear_clean_surf_rts',i_time,lw_flux_up_clear_clean_surf)    
    CALL write_nc(output_file,'lw_up_clear_clean_toa_rts',i_time,lw_flux_up_clear_clean_toa)
        
  END DO ! i_time


  ! Close ascii files
  IF (l_ascii) THEN
    CLOSE(iu_out_sw_lfric)
    CLOSE(iu_out_sw_offline)
    CLOSE(iu_out_lw_lfric)
    CLOSE(iu_out_lw_offline)
  END IF ! l_ascii  


  ! Deallocate everything

  IF (ALLOCATED(sw_heating_rate)) DEALLOCATE(sw_heating_rate)
  IF (ALLOCATED(lw_heating_rate)) DEALLOCATE(lw_heating_rate)
  IF (ALLOCATED(sw_flux_up)) DEALLOCATE(sw_flux_up)
  IF (ALLOCATED(sw_flux_down)) DEALLOCATE(sw_flux_down)
  IF (ALLOCATED(sw_flux_up_clear)) DEALLOCATE(sw_flux_up_clear)
  IF (ALLOCATED(sw_flux_down_clear)) DEALLOCATE(sw_flux_down_clear)
  IF (ALLOCATED(sw_flux_down_clear_clean_surf)) DEALLOCATE(sw_flux_down_clear_clean_surf)
  IF (ALLOCATED(sw_flux_up_clear_clean_surf)) DEALLOCATE(sw_flux_up_clear_clean_surf)
  IF (ALLOCATED(sw_flux_up_clear_clean_toa)) DEALLOCATE(sw_flux_up_clear_clean_toa)
  IF (ALLOCATED(lw_flux_up)) DEALLOCATE(lw_flux_up)
  IF (ALLOCATED(lw_flux_down)) DEALLOCATE(lw_flux_down)
  IF (ALLOCATED(lw_flux_up_clear)) DEALLOCATE(lw_flux_up_clear)
  IF (ALLOCATED(lw_flux_down_clear)) DEALLOCATE(lw_flux_down_clear)
  IF (ALLOCATED(lw_flux_down_clear_clean_surf)) DEALLOCATE(lw_flux_down_clear_clean_surf)
  IF (ALLOCATED(lw_flux_up_clear_clean_surf)) DEALLOCATE(lw_flux_up_clear_clean_surf)
  IF (ALLOCATED(lw_flux_up_clear_clean_toa)) DEALLOCATE(lw_flux_up_clear_clean_toa)

  IF (ALLOCATED(sw_heating_rate_rts_nc)) DEALLOCATE(sw_heating_rate_rts_nc)
  IF (ALLOCATED(lw_heating_rate_rts_nc)) DEALLOCATE(lw_heating_rate_rts_nc)
  IF (ALLOCATED(sw_up_rts_nc)) DEALLOCATE(sw_up_rts_nc)
  IF (ALLOCATED(sw_down_rts_nc)) DEALLOCATE(sw_down_rts_nc)
  IF (ALLOCATED(lw_up_rts_nc)) DEALLOCATE(lw_up_rts_nc)
  IF (ALLOCATED(lw_down_rts_nc)) DEALLOCATE(lw_down_rts_nc)
  IF (ALLOCATED(sw_up_clear_rts_nc)) DEALLOCATE(sw_up_clear_rts_nc)
  IF (ALLOCATED(sw_down_clear_rts_nc)) DEALLOCATE(sw_down_clear_rts_nc)
  IF (ALLOCATED(lw_up_clear_rts_nc)) DEALLOCATE(lw_up_clear_rts_nc)
  IF (ALLOCATED(lw_down_clear_rts_nc)) DEALLOCATE(lw_down_clear_rts_nc)
  IF (ALLOCATED(sw_down_clear_clean_surf_rts_nc)) DEALLOCATE(sw_down_clear_clean_surf_rts_nc)
  IF (ALLOCATED(sw_up_clear_clean_surf_rts_nc)) DEALLOCATE(sw_up_clear_clean_surf_rts_nc)
  IF (ALLOCATED(sw_up_clear_clean_toa_rts_nc)) DEALLOCATE(sw_up_clear_clean_toa_rts_nc)
  IF (ALLOCATED(lw_down_clear_clean_surf_rts_nc)) DEALLOCATE(lw_down_clear_clean_surf_rts_nc)
  IF (ALLOCATED(lw_up_clear_clean_surf_rts_nc)) DEALLOCATE(lw_up_clear_clean_surf_rts_nc)
  IF (ALLOCATED(lw_up_clear_clean_toa_rts_nc)) DEALLOCATE(lw_up_clear_clean_toa_rts_nc)

  IF (ALLOCATED(pressure_in_wth_nc)) DEALLOCATE(pressure_in_wth_nc)
  IF (ALLOCATED(temperature_in_wth_nc)) DEALLOCATE(temperature_in_wth_nc)
  IF (ALLOCATED(d_mass_nc)) DEALLOCATE(d_mass_nc)
  IF (ALLOCATED(rho_in_wth_nc)) DEALLOCATE(rho_in_wth_nc)
  IF (ALLOCATED(t_layer_boundaries_nc)) DEALLOCATE(t_layer_boundaries_nc)

  IF (ALLOCATED(ch4_nc)) DEALLOCATE(ch4_nc)
  IF (ALLOCATED(co_nc)) DEALLOCATE(co_nc)
  IF (ALLOCATED(co2_nc)) DEALLOCATE(co2_nc)
  IF (ALLOCATED(cs_nc)) DEALLOCATE(cs_nc)
  IF (ALLOCATED(h2_nc)) DEALLOCATE(h2_nc)
  IF (ALLOCATED(h2o_nc)) DEALLOCATE(h2o_nc)
  IF (ALLOCATED(hcn_nc)) DEALLOCATE(hcn_nc)
  IF (ALLOCATED(he_nc)) DEALLOCATE(he_nc)
  IF (ALLOCATED(potassium_nc)) DEALLOCATE(potassium_nc)
  IF (ALLOCATED(li_nc)) DEALLOCATE(li_nc)
  IF (ALLOCATED(n2_nc)) DEALLOCATE(n2_nc)
  IF (ALLOCATED(n2o_nc)) DEALLOCATE(n2o_nc)
  IF (ALLOCATED(na_nc)) DEALLOCATE(na_nc)
  IF (ALLOCATED(nh3_nc)) DEALLOCATE(nh3_nc)
  IF (ALLOCATED(o2_nc)) DEALLOCATE(o2_nc)
  IF (ALLOCATED(o3_nc)) DEALLOCATE(o3_nc)
  IF (ALLOCATED(rb_nc)) DEALLOCATE(rb_nc)
  IF (ALLOCATED(so2_nc)) DEALLOCATE(so2_nc)
  IF (ALLOCATED(tio_nc)) DEALLOCATE(tio_nc)
  IF (ALLOCATED(vo_nc)) DEALLOCATE(vo_nc)

  IF (ALLOCATED(tile_fraction_nc)) DEALLOCATE(tile_fraction_nc)
  IF (ALLOCATED(tile_temperature_nc)) DEALLOCATE(tile_temperature_nc)
  IF (ALLOCATED(tile_lw_albedo_nc)) DEALLOCATE(tile_lw_albedo_nc)
  IF (ALLOCATED(tile_sw_diffuse_albedo_nc)) DEALLOCATE(tile_sw_diffuse_albedo_nc)
  IF (ALLOCATED(tile_sw_direct_albedo_nc)) DEALLOCATE(tile_sw_direct_albedo_nc)
  IF (ALLOCATED(radiative_cloud_fraction_nc)) DEALLOCATE(radiative_cloud_fraction_nc)
  IF (ALLOCATED(liquid_fraction_nc)) DEALLOCATE(liquid_fraction_nc)
  IF (ALLOCATED(frozen_fraction_nc)) DEALLOCATE(frozen_fraction_nc)
  IF (ALLOCATED(mcl_nc)) DEALLOCATE(mcl_nc)
  IF (ALLOCATED(mci_nc)) DEALLOCATE(mci_nc)
  IF (ALLOCATED(n_ice_nc)) DEALLOCATE(n_ice_nc)
  IF (ALLOCATED(cloud_drop_no_conc_nc)) DEALLOCATE(cloud_drop_no_conc_nc)
  IF (ALLOCATED(radiative_conv_fraction_nc)) DEALLOCATE(radiative_conv_fraction_nc)
  IF (ALLOCATED(conv_liquid_fraction_nc)) DEALLOCATE(conv_liquid_fraction_nc)
  IF (ALLOCATED(conv_frozen_fraction_nc)) DEALLOCATE(conv_frozen_fraction_nc)
  IF (ALLOCATED(conv_liquid_mmr_nc)) DEALLOCATE(conv_liquid_mmr_nc)
  IF (ALLOCATED(conv_frozen_mmr_nc)) DEALLOCATE(conv_frozen_mmr_nc)
  IF (ALLOCATED(sigma_mc_nc)) DEALLOCATE(sigma_mc_nc)
  IF (ALLOCATED(layer_heat_capacity_nc)) DEALLOCATE(layer_heat_capacity_nc)
  IF (ALLOCATED(sulphuric_nc)) DEALLOCATE(sulphuric_nc)
  IF (ALLOCATED(cos_zenith_angle_rts_nc)) DEALLOCATE(cos_zenith_angle_rts_nc)
  IF (ALLOCATED(lit_fraction_rts_nc)) DEALLOCATE(lit_fraction_rts_nc)
  IF (ALLOCATED(stellar_irradiance_rts_nc)) DEALLOCATE(stellar_irradiance_rts_nc)
  IF (ALLOCATED(orographic_correction_rts_nc)) DEALLOCATE(orographic_correction_rts_nc)
  IF (ALLOCATED(rand_seed_real_radinput_nc)) DEALLOCATE(rand_seed_real_radinput_nc)

  IF (ALLOCATED(aer_mix_ratio_nc)) DEALLOCATE(aer_mix_ratio_nc)
  IF (ALLOCATED(aer_sw_absorption_nc)) DEALLOCATE(aer_sw_absorption_nc)
  IF (ALLOCATED(aer_sw_scattering_nc)) DEALLOCATE(aer_sw_scattering_nc)
  IF (ALLOCATED(aer_sw_asymmetry_nc)) DEALLOCATE(aer_sw_asymmetry_nc)
  IF (ALLOCATED(aer_lw_absorption_nc)) DEALLOCATE(aer_lw_absorption_nc)
  IF (ALLOCATED(aer_lw_scattering_nc)) DEALLOCATE(aer_lw_scattering_nc)
  IF (ALLOCATED(aer_lw_asymmetry_nc)) DEALLOCATE(aer_lw_asymmetry_nc)

  IF (ALLOCATED(rand_seed)) DEALLOCATE(rand_seed)
  IF (ALLOCATED(profile_list)) DEALLOCATE(profile_list)
  IF (ALLOCATED(n_cloud_layer)) DEALLOCATE(n_cloud_layer)
  
  IF (ALLOCATED(dummy_nc)) DEALLOCATE(dummy_nc)
  IF (ALLOCATED(dummy_arr_out)) DEALLOCATE(dummy_arr_out)
    
  ! Nullify pointers
    
  NULLIFY(sw_heating_rate_rts)
  NULLIFY(lw_heating_rate_rts)
  NULLIFY(sw_up_rts)
  NULLIFY(sw_down_rts)
  NULLIFY(lw_up_rts)
  NULLIFY(lw_down_rts)
  NULLIFY(sw_up_clear_rts)
  NULLIFY(sw_down_clear_rts)
  NULLIFY(lw_up_clear_rts)
  NULLIFY(lw_down_clear_rts)
  NULLIFY(sw_down_clear_clean_surf_rts)
  NULLIFY(sw_up_clear_clean_surf_rts)
  NULLIFY(sw_up_clear_clean_toa_rts)
  NULLIFY(lw_down_clear_clean_surf_rts)
  NULLIFY(lw_up_clear_clean_surf_rts)
  NULLIFY(lw_up_clear_clean_toa_rts)
  
  NULLIFY(pressure_in_wth)
  NULLIFY(temperature_in_wth)
  NULLIFY(d_mass)
  NULLIFY(rho_in_wth)
  NULLIFY(t_layer_boundaries)
  
  NULLIFY(ch4)
  NULLIFY(co)
  NULLIFY(co2)
  NULLIFY(cs)
  NULLIFY(h2)
  NULLIFY(h2o)
  NULLIFY(hcn)
  NULLIFY(he)
  NULLIFY(potassium)
  NULLIFY(li)
  NULLIFY(n2)
  NULLIFY(n2o)
  NULLIFY(na)
  NULLIFY(nh3)
  NULLIFY(o2)
  NULLIFY(o3)
  NULLIFY(rb)
  NULLIFY(so2)
  NULLIFY(tio)
  NULLIFY(vo)

  NULLIFY(tile_fraction)
  NULLIFY(tile_temperature)
  NULLIFY(tile_lw_albedo)
  NULLIFY(tile_sw_diffuse_albedo)
  NULLIFY(tile_sw_direct_albedo)
  NULLIFY(radiative_cloud_fraction)
  NULLIFY(liquid_fraction)
  NULLIFY(frozen_fraction)
  NULLIFY(mcl)
  NULLIFY(mci)
  NULLIFY(n_ice)
  NULLIFY(cloud_drop_no_conc)
  NULLIFY(radiative_conv_fraction)
  NULLIFY(conv_liquid_fraction)
  NULLIFY(conv_frozen_fraction)
  NULLIFY(conv_liquid_mmr)
  NULLIFY(conv_frozen_mmr)
  NULLIFY(sigma_mc)
  NULLIFY(layer_heat_capacity)
  NULLIFY(sulphuric)
  NULLIFY(cos_zenith_angle_rts)
  NULLIFY(lit_fraction_rts)
  NULLIFY(stellar_irradiance_rts)
  NULLIFY(orographic_correction_rts)
  NULLIFY(rand_seed_real_radinput)

  NULLIFY(aer_mix_ratio)
  NULLIFY(aer_sw_absorption)
  NULLIFY(aer_sw_scattering)
  NULLIFY(aer_sw_asymmetry)
  NULLIFY(aer_lw_absorption)
  NULLIFY(aer_lw_scattering)
  NULLIFY(aer_lw_asymmetry)

  IF (l_debug) THEN
    WRITE(*,*) 'PROGRAM runes_nc end: '
  END IF

END PROGRAM runes_nc

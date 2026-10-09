! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate grey optical properties.
!
! Method:
!   For each activated optical process, excluding gaseous
!   absorption, increments are calculated for the total and
!   scattering extinctions, and the products of the asymmetry
!   factor and the forward scattering factor in clear and
!   cloudy regions. These increments are summed, and the grey
!   total and scattering extinctions and the asymmetry and forward
!   scattering factors are thus calculated.
!
!- ---------------------------------------------------------------------
MODULE grey_opt_prop_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'GREY_OPT_PROP_MOD'
CONTAINS
SUBROUTINE grey_opt_prop(ierr, control, radout, i_band                  &
    , n_profile, n_layer, p, t, density                                 &
    , n_order_phase, l_solar_phf, n_direction, cos_sol_view             &
    , i_rayleigh_scheme, n_gas_rayleigh, index_rayleigh                 &
    , rayleigh_coeff_tot, rayleigh_coeff_gas                            &
    , gas_mix_ratio                                                     &
    , l_continuum, n_continuum, k_continuum, amount_continuum           &
    , n_aerosol_mr, aerosol_mix_ratio                                   &
    , aerosol_mr_source, aerosol_mr_type_index                          &
    , i_aerosol_parametrization                                         &
    , i_humidity_pointer, humidities, delta_humidity                    &
    , mean_rel_humidity                                                 &
    , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc         &
    , n_ukca_mode, ukca_modal_mixr                                      &
    , ukca_absorption, ukca_scattering, ukca_asymmetry                  &
    , n_opt_level_aerosol_prsc, aerosol_pressure_prsc                   &
    , aerosol_absorption_prsc, aerosol_scattering_prsc                  &
    , aerosol_phase_fnc_prsc                                            &
    , n_cloud_profile, i_cloud_profile                                  &
    , n_cloud_top, n_condensed, l_cloud_cmp, i_phase_cmp                &
    , i_condensed_param, condensed_param_list                           &
    , condensed_mix_ratio, condensed_dim_char                           &
    , condensed_rel_var_dens                                            &
    , n_cloud_type, i_cloud_type                                        &
    , n_opt_level_drop_prsc                                             &
    , drop_pressure_prsc, drop_absorption_prsc                          &
    , drop_scattering_prsc, drop_phase_fnc_prsc                         &
    , n_opt_level_ice_prsc, ice_pressure_prsc                           &
    , ice_absorption_prsc, ice_scattering_prsc, ice_phase_fnc_prsc      &
    , ss_prop                                                           &
    , frac_cloud                                                        &
    , cloud_extinction, cloud_absorptivity                              &
    , ls_cloud_extinction, ls_cloud_absorptivity                        &
    , cnv_cloud_extinction, cnv_cloud_absorptivity                      &
    , nd_profile, nd_radiance_profile, nd_layer                         &
    , nd_layer_clr, id_ct                                               &
    , nd_continuum, nd_species                                          &
    , nd_aerosol_species, nd_aerosol_mixratio                           &
    , nd_humidities, nd_cloud_parameter, nd_cloud_component             &
    , nd_cloud_type, nd_phase_term                                      &
    , nd_phf_term_aerosol_prsc, nd_phf_term_cloud_prsc                  &
    , nd_max_order, nd_direction                                        &
    , nd_ukca_mode, nd_profile_aerosol_prsc, nd_profile_cloud_prsc      &
    , nd_opt_level_aerosol_prsc, nd_opt_level_cloud_prsc                &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_out,     ONLY: StrOut
  USE def_ss_prop, ONLY: str_ss_prop
  USE rad_pcf, ONLY: i_err_fatal, ip_rayleigh_total, ip_rayleigh_custom,&
                     ip_cairns, ip_cloud_mcica, ip_cloud_type_si,       &
                     ip_cloud_type_sw, ip_direct_csr_scaling,           &
                     ip_phase_ice, ip_phase_water
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE circumsolar_fraction_mod, ONLY: circumsolar_fraction
  USE opt_prop_aerosol_mod, ONLY: opt_prop_aerosol
  USE opt_prop_ice_cloud_mod, ONLY: opt_prop_ice_cloud
  USE opt_prop_inhom_corr_cairns_mod, ONLY: opt_prop_inhom_corr_cairns
  USE opt_prop_ukca_aerosol_mod, ONLY: opt_prop_ukca_aerosol
  USE opt_prop_water_cloud_mod, ONLY: opt_prop_water_cloud

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Output fields:
  TYPE(StrOut),       INTENT(INOUT) :: radout

! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_radiance_profile                                               &
!       Size allocated for profiles of quantities specifically
!       used in calulating radiances
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_layer_clr                                                      &
!       Size allocated for completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_species                                                        &
!       Size allocated for nd_species
    , nd_aerosol_species                                                &
!       Size allocated for aerosols in spectral information
    , nd_aerosol_mixratio                                               &
!       Size allocated for aerosols in aerosol_mix_ratio array
    , nd_humidities                                                     &
!       Size allocated for humidities
    , nd_continuum                                                      &
!       Size allocated for continua
    , nd_phase_term                                                     &
!       Size allocated for terms in the phase function
    , nd_phf_term_aerosol_prsc                                          &
!       Size allocated for terms in the prescribed aerosol phase function
    , nd_phf_term_cloud_prsc                                            &
!       Size allocated for terms in the prescribed cloud phase function
    , nd_max_order                                                      &
!       Size allocated for the order of the calculation
    , nd_cloud_parameter                                                &
!       Size allocated for cloud parameters
    , nd_cloud_component                                                &
!       Size allocated for components of clouds
    , nd_cloud_type                                                     &
!       Size allocated for types of clouds
    , nd_ukca_mode                                                      &
!       Size allocated for UKCA aerosols
    , nd_profile_aerosol_prsc                                           &
!       Size allocated for profiles of prescribed
!       cloudy optical properties
    , nd_profile_cloud_prsc                                             &
!       Size allocated for profiles of prescribed
!       aerosol optical properties
    , nd_opt_level_aerosol_prsc                                         &
!       Size allocated for levels of prescribed
!       cloudy optical properties
    , nd_opt_level_cloud_prsc
!       Size allocated for levels of prescribed
!       aerosol optical properties


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

  INTEGER, INTENT(IN) ::                                                &
      i_band
!       Spectral band

! Basic atmospheric properties:
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers

  REAL (RealK), INTENT(IN) ::                                           &
      p(nd_profile, nd_layer)                                           &
!       Pressure
    , t(nd_profile, nd_layer)                                           &
!       Temperature
    , density(nd_profile, nd_layer)
!       Density at levels


! Optical switches:
  LOGICAL, INTENT(IN) ::                                                &
      l_solar_phf
!       Flag to use an extended phase function for solar radiation
  INTEGER, INTENT(IN) ::                                                &
      n_order_phase
!       Order of terms in the phase function

! Directional information
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      cos_sol_view(nd_radiance_profile, nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing direction


! Rayleigh scattering:
  INTEGER, INTENT(IN) ::                                                &
      i_rayleigh_scheme                                                 &
!       Rayleigh scattering scheme
    , n_gas_rayleigh                                                    &
!       Number of gases to include in calculation of Rayleigh scattering
!       coefficient
    , index_rayleigh(nd_species)
!       Index of gases to include in calculation of Rayleigh scattering
!       coefficient
  REAL (RealK), INTENT(IN) ::                                           &
      rayleigh_coeff_tot                                                &
!       Rayleigh scattering coefficient for total gas
    , rayleigh_coeff_gas(nd_species)                                    &
!       Rayleigh scattering coefficient for each gas
    , gas_mix_ratio(nd_profile, nd_layer, nd_species)
!       Gaseous mass mixing ratios

! Continuum processes:
  LOGICAL, INTENT(IN) ::                                                &
      l_continuum
!       Continuum absorption activated

  INTEGER, INTENT(IN) ::                                                &
      n_continuum
!       Number of continua

  REAL (RealK), INTENT(IN) ::                                           &
      k_continuum(nd_continuum)                                         &
!       Continuum extinction
    , amount_continuum(nd_profile, nd_layer, nd_continuum)
!       Amounts for continua


! Properties of aerosols:
  INTEGER, INTENT(IN) ::                                                &
      n_aerosol_mr                                                      &
!       Number of aerosol species in aerosol_mix_ratio array
    , aerosol_mr_type_index(nd_aerosol_mixratio)                        &
!       Index relating aerosol_mix_ratio aerosols to aerosols in
!       the spectral information
    , aerosol_mr_source(nd_aerosol_mixratio)                            &
!       Scheme/source of the aerosol data, to determine use in
!       changing radiative fluxes and use in diagnostics
    , i_aerosol_parametrization(nd_aerosol_species)                     &
!       Parametrizations of aerosols
    , i_humidity_pointer(nd_profile,  nd_layer)
!       Pointer to aerosol look-up table

  REAL (RealK), INTENT(IN) ::                                           &
      aerosol_mix_ratio(nd_profile, nd_layer, nd_aerosol_mixratio)      &
!       Number densty of aerosols
    , aerosol_absorption(nd_humidities, nd_aerosol_species)             &
!       Aerosol absorption in band for a mixing ratio of unity
    , aerosol_scattering(nd_humidities, nd_aerosol_species)             &
!       Aerosol scattering in band for a mixing ratio of unity
    , aerosol_phase_fnc(nd_humidities                                   &
        , nd_phase_term, nd_aerosol_species)                            &
!       Aerosol phase function in band
    , humidities(nd_humidities, nd_aerosol_species)                     &
!       Array of humidities
    , delta_humidity                                                    &
!       Increment in humidity
    , mean_rel_humidity(nd_profile, nd_layer)
!       Mixing ratio of water vapour


! Properties of UKCA aerosols:
  INTEGER, INTENT(IN) ::                                                &
      n_ukca_mode
!       Actual number of UKCA aerosol modes

  REAL (RealK), INTENT(IN) ::                                           &
      ukca_modal_mixr(nd_profile, nd_layer, MAX(1,nd_ukca_mode))        &
!       Modal mass mixing ratio
    , ukca_absorption(nd_profile, nd_layer, MAX(1,nd_ukca_mode))        &
!       Waveband-averaged specific coefficient for absorption
    , ukca_scattering(nd_profile, nd_layer, MAX(1,nd_ukca_mode))        &
!       Waveband-averaged specific coefficient for scattering
    , ukca_asymmetry(nd_profile, nd_layer, MAX(1,nd_ukca_mode))
!       Waveband-averaged asymmetry factor

! Observational properties of aerosols:
  INTEGER, INTENT(IN) ::                                                &
      n_opt_level_aerosol_prsc(nd_aerosol_species)
!       Number of levels of prescribed optical properties
!       of aerosols

  REAL  (RealK), INTENT(IN) ::                                          &
      aerosol_pressure_prsc(nd_profile_aerosol_prsc                     &
        , nd_opt_level_aerosol_prsc, nd_aerosol_species)                &
!       Pressures at which optical properties of aerosols
!       are prescribed
    , aerosol_absorption_prsc(nd_profile_aerosol_prsc                   &
        , nd_opt_level_aerosol_prsc, nd_aerosol_species)                &
!       Prescribed absorption by aerosols
    , aerosol_scattering_prsc(nd_profile_aerosol_prsc                   &
        , nd_opt_level_aerosol_prsc, nd_aerosol_species)                &
!       Prescribed scattering by aerosols
    , aerosol_phase_fnc_prsc(nd_profile_aerosol_prsc                    &
        , nd_opt_level_aerosol_prsc                                     &
        , nd_phf_term_aerosol_prsc, nd_aerosol_species)
!       Prescribed phase functions of aerosols


! Properties of clouds:

! Geometry of clouds:
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_cloud_type                                                      &
!       Number of types of clouds
    , n_cloud_profile(id_ct: nd_layer)                                  &
!       Number of cloudy profiles in each layer
    , i_cloud_profile(nd_profile, id_ct: nd_layer)                      &
!       Profiles containing clouds
    , i_cloud_type(nd_cloud_component)
!       Types of cloud to which each component contributes

! Microphysical quantities:
  INTEGER, INTENT(IN) ::                                                &
      n_condensed                                                       &
!       Number of condensed components
    , i_phase_cmp(nd_cloud_component)                                   &
!       Phases of cloudy components
    , i_condensed_param(nd_cloud_component)
!       Parametrization schemes for cloudy components

  LOGICAL, INTENT(IN) ::                                                &
      l_cloud_cmp(nd_cloud_component)
!       Flags to activate cloudy components

  REAL (RealK), INTENT(IN) ::                                           &
      condensed_param_list(nd_cloud_parameter                           &
        , nd_cloud_component)                                           &
!       Coefficients in parametrization schemes
    , condensed_mix_ratio(nd_profile, id_ct: nd_layer                   &
        , nd_cloud_component)                                           &
!       Mixing ratios of cloudy components
    , condensed_dim_char(nd_profile, id_ct: nd_layer                    &
        , nd_cloud_component)                                           &
!       Characteristic dimensions of cloudy components
    , condensed_rel_var_dens(nd_profile, id_ct: nd_layer                &
        , nd_cloud_component)
!       Relative variance of cloud density

! Prescribed cloudy optical properties:
  INTEGER, INTENT(IN) ::                                                &
      n_opt_level_drop_prsc                                             &
!       Number of levels of prescribed
!       optical properties of droplets
    , n_opt_level_ice_prsc
!       Number of levels of prescribed
!       optical properties of ice crystals

  REAL  (RealK), INTENT(IN) ::                                          &
      drop_pressure_prsc(nd_profile_cloud_prsc                          &
        , nd_opt_level_cloud_prsc)                                      &
!       Pressures at which optical properties of
!       droplets are prescribed
    , drop_absorption_prsc(nd_profile_cloud_prsc                        &
        , nd_opt_level_cloud_prsc)                                      &
!       Prescribed absorption by droplets
    , drop_scattering_prsc(nd_profile_cloud_prsc                        &
        , nd_opt_level_cloud_prsc)                                      &
!       Prescribed scattering by droplets
    , drop_phase_fnc_prsc(nd_profile_cloud_prsc                         &
        , nd_opt_level_cloud_prsc, nd_phf_term_cloud_prsc)              &
!       Prescribed phase function of droplets
    , ice_pressure_prsc(nd_profile_cloud_prsc                           &
        , nd_opt_level_cloud_prsc)                                      &
!       Pressures at which optical properties of
!       ice crystals are prescribed
    , ice_absorption_prsc(nd_profile_cloud_prsc                         &
        , nd_opt_level_cloud_prsc)                                      &
!       Prescribed absorption by ice crystals
    , ice_scattering_prsc(nd_profile_cloud_prsc                         &
        , nd_opt_level_cloud_prsc)                                      &
!       Prescribed scattering by ice crystals
    , ice_phase_fnc_prsc(nd_profile_cloud_prsc                          &
        , nd_opt_level_cloud_prsc, nd_phf_term_cloud_prsc)
!       Prescribed phase functions of ice crystals

  REAL (RealK), INTENT(IN) ::                                           &
      frac_cloud(nd_profile, id_ct:nd_layer, nd_cloud_type)
!       Fractions of each type of cloud


! Optical properties
  TYPE(str_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

  REAL (RealK), INTENT(OUT) ::                                          &
      cloud_extinction(nd_profile, nd_layer)                            &
!       Mean cloud extinction
    , cloud_absorptivity(nd_profile, nd_layer)                          &
!       Mean cloud absorptivity
    , ls_cloud_extinction(nd_profile, nd_layer)                         &
!       Mean large-scale cloud extinction
    , ls_cloud_absorptivity(nd_profile, nd_layer)                       &
!       Mean large-scale cloud absorptivity
    , cnv_cloud_extinction(nd_profile, nd_layer)                        &
!       Mean convective cloud extinction
    , cnv_cloud_absorptivity(nd_profile, nd_layer)
!       Mean convective cloud absorptivity



! Local variables.
  INTEGER ::                                                            &
      l                                                                 &
!       Loop variable
    , ll                                                                &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , ls                                                                &
!       Loop variable
    , n_index                                                           &
!       Number of indices satisfying the test
    , indx(nd_profile)
!       Indices satifying the test
  REAL (RealK) ::                                                       &
      rayleigh_coeff(nd_profile, nd_layer)
!       Calculated total Rayleigh coefficient

! Temporary variable for the divisions
  REAL (RealK) :: tmp_inv(nd_profile)

  REAL (RealK), PARAMETER :: tiny_k=TINY(ss_prop%k_ext_scat)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER(LEN=*), PARAMETER :: RoutineName='GREY_OPT_PROP'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If using a separate solar phase function that must be initialized.
  IF (l_solar_phf) THEN
    DO id=1, n_direction
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%phase_fnc_solar_clr(l, i, id)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%phase_fnc_solar(l, i, id, 0)=0.0_RealK
        END DO
      END DO
    END DO
    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        ss_prop%forward_solar_clr(l, i)=0.0_RealK
      END DO
    END DO
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        ss_prop%forward_solar(l, i, 0)=0.0_RealK
      END DO
    END DO
  END IF



! Consider each optical process in turn.

! Rayleigh scattering:

  IF (control%l_rayleigh) THEN

!   Calculate Rayleigh scattering coefficient for current band
    SELECT CASE (i_rayleigh_scheme)

    CASE (ip_rayleigh_total)
!     Rayleigh scattering coefficients are tabulated for total gas
      rayleigh_coeff(1:n_profile, 1:n_layer) = rayleigh_coeff_tot
    
    CASE (ip_rayleigh_custom)
!     Compute Rayleigh scattering coefficient of total gas from individual
!     coefficients for each gas ignoring non-ideal gas effects due to water
!     vapour
      DO l=1, n_profile
        DO i=1, n_layer
          rayleigh_coeff(l, i) &
            =SUM(rayleigh_coeff_gas(1:n_gas_rayleigh) &
            *gas_mix_ratio(l, i, index_rayleigh(1:n_gas_rayleigh)))
        END DO
      END DO

    CASE DEFAULT
      WRITE(cmessage, '(a)')                                            &
        'i_rayleigh_scheme has been set to an illegal value.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    
    END SELECT

!   Forward scattering is required only when delta-rescaling
!   is performed.
    IF (control%l_rescale) THEN
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=rayleigh_coeff(l, i)
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
          ss_prop%forward_scatter_clr(l, i)=0.0_RealK
          ss_prop%forward_scatter_clr_csr(l, i)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=rayleigh_coeff(l, i)
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
          ss_prop%forward_scatter(l, i, 0)=0.0_RealK
          ss_prop%forward_scatter_csr(l, i, 0)=0.0_RealK
        END DO
      END DO
    ELSE
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=rayleigh_coeff(l, i)
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=rayleigh_coeff(l, i)
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
        END DO
      END DO
    END IF

!   Only the second Lengendre polynomial contributes.
    IF (n_order_phase >= 2) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%phase_fnc_clr(l, i, 2)=rayleigh_coeff(l, i)*1.0e-01_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%phase_fnc(l, i, 2, 0)=rayleigh_coeff(l, i)*1.0e-01_RealK
        END DO
      END DO
      DO ls=3, n_order_phase
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            ss_prop%phase_fnc_clr(l, i, ls)=0.0_RealK
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%phase_fnc(l, i, ls, 0)=0.0_RealK
          END DO
        END DO
      END DO
    END IF

!   No formal rescaling is applied to the phase function for
!   Rayleigh scattering, as only g_2 is non-zero.

    IF (l_solar_phf) THEN

      DO id=1, n_direction
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            ss_prop%phase_fnc_solar_clr(l, i, id)                       &
              =ss_prop%phase_fnc_solar_clr(l, i, id)                    &
              +rayleigh_coeff(l, i)                                     &
              *0.75_RealK*(1.0_RealK+cos_sol_view(l, id)**2)
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%phase_fnc_solar(l, i, id, 0)                        &
              =ss_prop%phase_fnc_solar(l, i, id, 0)                     &
              +rayleigh_coeff(l, i)                                     &
              *0.75_RealK*(1.0_RealK+cos_sol_view(l, id)**2)
          END DO
        END DO
      END DO

    END IF

  ELSE

    IF (control%l_rescale) THEN
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=0.0_RealK
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
          ss_prop%forward_scatter_clr(l, i)=0.0_RealK
          ss_prop%forward_scatter_clr_csr(l, i)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=0.0_RealK
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
          ss_prop%forward_scatter(l, i, 0)=0.0_RealK
          ss_prop%forward_scatter_csr(l, i, 0)=0.0_RealK
        END DO
      END DO
    ELSE
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=0.0_RealK
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=0.0_RealK
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
        END DO
      END DO
    END IF
    DO ls=2, n_order_phase
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%phase_fnc_clr(l, i, ls)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%phase_fnc(l, i, ls, 0)=0.0_RealK
        END DO
      END DO
    END DO

  END IF

  IF (control%l_aerosol) THEN
!   Include the effects of aerosol.
!   Above clouds.
    CALL opt_prop_aerosol(ierr                                          &
      , n_profile, 1, n_cloud_top-1                                     &
      , n_order_phase, control%l_rescale, control%n_order_forward       &
      , control%l_henyey_greenstein_pf                                  &
      , n_aerosol_mr, aerosol_mix_ratio                                 &
      , aerosol_mr_source, aerosol_mr_type_index                        &
      , i_aerosol_parametrization                                       &
      , i_humidity_pointer, humidities, delta_humidity                  &
      , mean_rel_humidity                                               &
      , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc       &
      , l_solar_phf, control%n_order_phase_solar, n_direction           &
      , cos_sol_view                                                    &
      , p, density                                                      &
      , n_opt_level_aerosol_prsc, aerosol_pressure_prsc                 &
      , aerosol_absorption_prsc, aerosol_scattering_prsc                &
      , aerosol_phase_fnc_prsc                                          &
      , ss_prop%k_grey_tot_clr, ss_prop%k_ext_scat_clr                  &
      , ss_prop%phase_fnc_clr, ss_prop%forward_scatter_clr              &
      , ss_prop%forward_solar_clr, ss_prop%phase_fnc_solar_clr          &
      , nd_profile, nd_radiance_profile, nd_layer                       &
      , 1, nd_layer_clr                                                 &
      , nd_aerosol_species, nd_aerosol_mixratio, nd_humidities          &
      , nd_phase_term, nd_max_order, nd_direction                       &
      , nd_profile_aerosol_prsc, nd_opt_level_aerosol_prsc              &
      , nd_phf_term_aerosol_prsc                                        &
      )
!   Within clouds:
    CALL opt_prop_aerosol(ierr                                          &
      , n_profile, n_cloud_top, n_layer                                 &
      , n_order_phase, control%l_rescale, control%n_order_forward       &
      , control%l_henyey_greenstein_pf                                  &
      , n_aerosol_mr, aerosol_mix_ratio                                 &
      , aerosol_mr_source, aerosol_mr_type_index                        &
      , i_aerosol_parametrization                                       &
      , i_humidity_pointer, humidities, delta_humidity                  &
      , mean_rel_humidity                                               &
      , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc       &
      , l_solar_phf, control%n_order_phase_solar, n_direction           &
      , cos_sol_view                                                    &
      , p, density                                                      &
      , n_opt_level_aerosol_prsc, aerosol_pressure_prsc                 &
      , aerosol_absorption_prsc, aerosol_scattering_prsc                &
      , aerosol_phase_fnc_prsc                                          &
      , ss_prop%k_grey_tot(:, :, 0)                                     &
      , ss_prop%k_ext_scat(:, :, 0)                                     &
      , ss_prop%phase_fnc(:, :, :, 0)                                   &
      , ss_prop%forward_scatter(:, :, 0)                                &
      , ss_prop%forward_solar(:, :, 0)                                  &
      , ss_prop%phase_fnc_solar(:, :, :, 0)                             &
      , nd_profile, nd_radiance_profile, nd_layer                       &
      , id_ct, nd_layer                                                 &
      , nd_aerosol_species, nd_aerosol_mixratio, nd_humidities          &
      , nd_phase_term, nd_max_order, nd_direction                       &
      , nd_profile_aerosol_prsc, nd_opt_level_aerosol_prsc              &
      , nd_phf_term_aerosol_prsc                                        &
      )
  END IF

  IF (control%l_aerosol_mode) THEN
!   Include the effects of UKCA aerosols.
!   Above clouds.
    CALL opt_prop_ukca_aerosol(                                         &
        n_profile, 1, n_cloud_top-1                                     &
      , n_order_phase, control%l_rescale, control%n_order_forward       &
      , n_ukca_mode, ukca_modal_mixr                                    &
      , ukca_absorption, ukca_scattering, ukca_asymmetry                &
      , ss_prop%k_grey_tot_clr, ss_prop%k_ext_scat_clr                  &
      , ss_prop%phase_fnc_clr, ss_prop%forward_scatter_clr              &
      , nd_profile, nd_layer, 1, nd_layer_clr, nd_ukca_mode             &
      , nd_max_order                                                    &
      )
!   Within clouds:
    CALL opt_prop_ukca_aerosol(                                         &
        n_profile, n_cloud_top, n_layer                                 &
      , n_order_phase, control%l_rescale, control%n_order_forward       &
      , n_ukca_mode, ukca_modal_mixr                                    &
      , ukca_absorption, ukca_scattering, ukca_asymmetry                &
      , ss_prop%k_grey_tot(:, :, 0)                                     &
      , ss_prop%k_ext_scat(:, :, 0)                                     &
      , ss_prop%phase_fnc(:, :, :, 0)                                   &
      , ss_prop%forward_scatter(:, :, 0)                                &
      , nd_profile, nd_layer, id_ct, nd_layer, nd_ukca_mode             &
      , nd_max_order                                                    &
      )
  END IF

! Output aerosol optical property diagnostics
  IF (control%l_aerosol_absorption_band) THEN
    DO i=1, n_cloud_top-1
      DO l=1, n_profile
!       At this point k_grey_tot includes the absorption only
        radout%aerosol_absorption_band(l, i, i_band)                    &
          = ss_prop%k_grey_tot_clr(l, i)
      END DO
    END DO
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        radout%aerosol_absorption_band(l, i, i_band)                    &
          = ss_prop%k_grey_tot(l, i, 0)
      END DO
    END DO
  END IF
  IF (control%l_aerosol_scattering_band) THEN
    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        radout%aerosol_scattering_band(l, i, i_band)                    &
          = ss_prop%k_ext_scat_clr(l, i)
      END DO
    END DO
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        radout%aerosol_scattering_band(l, i, i_band)                    &
          = ss_prop%k_ext_scat(l, i, 0)
      END DO
    END DO
    IF (control%l_rayleigh) THEN
      DO i=1, n_layer
        DO l=1, n_profile
!         Remove the Rayleigh scattering contribution
          radout%aerosol_scattering_band(l,i,i_band) =                  &
            radout%aerosol_scattering_band(l,i,i_band) -                &
            rayleigh_coeff(l, i)
        END DO
      END DO
    END IF
  END IF
  IF (control%l_aerosol_asymmetry_band) THEN
!   The first moment of the phase function is weighted by the scattering
!   at this point: the diagnostic will be passed out as a weighted value
    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        radout%aerosol_asymmetry_band(l, i, i_band)                    &
          = ss_prop%phase_fnc_clr(l, i, 1)
      END DO
    END DO
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        radout%aerosol_asymmetry_band(l, i, i_band)                    &
          = ss_prop%phase_fnc(l, i, 1, 0)
      END DO
    END DO
  END IF

  IF (l_continuum) THEN
!   Include continuum absorption.
    DO j=1, n_continuum
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=ss_prop%k_grey_tot_clr(l, i)     &
            +k_continuum(j)*amount_continuum(l, i, j)
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=ss_prop%k_grey_tot(l, i, 0)       &
            +k_continuum(j)*amount_continuum(l, i, j)
        END DO
      END DO
    END DO
  END IF


! Add the scattering on to the total extinction. The final clear-sky
! phase function not calculated here since the product of the phase
! function and scattering is also needed to calculate the cloudy
! phase function.
  DO i=1, n_cloud_top-1
    DO l=1, n_profile
      ss_prop%k_grey_tot_clr(l, i)=ss_prop%k_grey_tot_clr(l, i)         &
        +ss_prop%k_ext_scat_clr(l, i)
    END DO
  END DO
  DO i=n_cloud_top, n_layer
    DO l=1, n_profile
      ss_prop%k_grey_tot(l, i, 0)=ss_prop%k_grey_tot(l, i, 0)           &
        +ss_prop%k_ext_scat(l, i, 0)
    END DO
  END DO


! If there are no clouds calculate the final optical properties
! and return to the calling routine.
  IF (.NOT.control%l_cloud) THEN

    IF (control%l_rescale) THEN
      DO i=1, n_cloud_top-1
        n_index=0
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
            n_index=n_index+1
            indx(n_index)=l
          END IF
        END DO
        DO l=1, n_index 
          tmp_inv(l)=1.0_RealK/ss_prop%k_ext_scat_clr(indx(l), i)
          ss_prop%forward_scatter_clr(indx(l), i)                       &
            =ss_prop%forward_scatter_clr(indx(l), i)*tmp_inv(l)
          DO ls=1, n_order_phase
            ss_prop%phase_fnc_clr(indx(l), i, ls)                       &
              =ss_prop%phase_fnc_clr(indx(l), i, ls)*tmp_inv(l)
          END DO
        END DO
!----------------------------------------------------------------------
! Calculate forward scattering fraction due to the CircumSolar Radiation
! (CSR) within a FOV of pyrheliometer. This fraction replaces g^2 
! in the Delta-Eddington scaling in the direct flux calculation.
!---------------------------------------------------------------------- 
        IF (control%i_direct_tau == ip_direct_csr_scaling ) THEN
! Above cloud top.
          CALL circumsolar_fraction(n_index                             &
           , indx, control%half_angle                                   &
           , ss_prop%phase_fnc_clr(:, i, 1)                             &
           , ss_prop%forward_scatter_clr_csr(:, i)                      &
           , nd_profile                                                 &
           )
        END IF
      END DO
      DO i=n_cloud_top, n_layer
        n_index=0
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
            n_index=n_index+1
            indx(n_index)=l
          END IF
        END DO
        DO l=1, n_index
          tmp_inv(l)=1.0_RealK/ss_prop%k_ext_scat(indx(l), i, 0)
          ss_prop%forward_scatter(indx(l), i, 0)                         &
            =ss_prop%forward_scatter(indx(l), i, 0)*tmp_inv(l)
          DO ls=1, n_order_phase
            ss_prop%phase_fnc(indx(l), i, ls, 0)                         &
              =ss_prop%phase_fnc(indx(l), i, ls, 0)*tmp_inv(l)
          END DO
        END DO
        IF (control%i_direct_tau == ip_direct_csr_scaling ) THEN
  ! Below cloud top.
          CALL circumsolar_fraction(n_index                               &
           , indx, control%half_angle                                     &
           , ss_prop%phase_fnc(:, i, 1, 0)                                &
           , ss_prop%forward_scatter_csr(:, i, 0)                         &
           , nd_profile                                                   &
           )
        END IF 
      END DO
    ELSE
      DO ls=1, n_order_phase
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
              ss_prop%phase_fnc_clr(l, i, ls)                           &
                =ss_prop%phase_fnc_clr(l, i, ls)                        &
                /ss_prop%k_ext_scat_clr(l, i)
            END IF
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
              ss_prop%phase_fnc(l, i, ls, 0)                            &
                =ss_prop%phase_fnc(l, i, ls, 0)                         &
                /ss_prop%k_ext_scat(l, i, 0)
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (l_solar_phf) THEN

      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
            ss_prop%forward_solar_clr(l, i)                             &
              =ss_prop%forward_solar_clr(l, i)                          &
              /ss_prop%k_ext_scat_clr(l, i)
          END IF
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
            ss_prop%forward_solar(l, i, 0)                              &
              =ss_prop%forward_solar(l, i, 0)                           &
              /ss_prop%k_ext_scat(l, i, 0)
          END IF
        END DO
      END DO

      DO id=1, n_direction
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
              ss_prop%phase_fnc_solar_clr(l, i, id)                     &
                =ss_prop%phase_fnc_solar_clr(l, i, id)                  &
                /ss_prop%k_ext_scat_clr(l, i)
            END IF
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
              ss_prop%phase_fnc_solar(l, i, id, 0)                      &
                =ss_prop%phase_fnc_solar(l, i, id, 0)                   &
                /ss_prop%k_ext_scat(l, i, 0)
            END IF
          END DO
        END DO
      END DO

    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN

  END IF



! Addition of cloudy properties:


! Add in background contibutions:


! All the processes occurring outside clouds also occur within them.
  DO k=1, n_cloud_type
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        ss_prop%k_grey_tot(l, i, k)=ss_prop%k_grey_tot(l, i, 0)
        ss_prop%k_ext_scat(l, i, k)=ss_prop%k_ext_scat(l, i, 0)
        ss_prop%forward_scatter(l, i, k)                                &
          =ss_prop%forward_scatter(l, i, 0)
        ss_prop%forward_scatter_csr(l, i, k)                            &
          =ss_prop%forward_scatter_csr(l, i, 0)
      END DO
      DO ls=1, n_order_phase
        DO l=1, n_profile
          ss_prop%phase_fnc(l, i, ls, k)                                &
            =ss_prop%phase_fnc(l, i, ls, 0)
        END DO
      END DO
    END DO
!   If using a separate solar phase function that must be initialized.
    IF (l_solar_phf) THEN
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%forward_solar(l, i, k)                                &
            =ss_prop%forward_solar(l, i, 0)
        END DO
      END DO
      DO id=1, n_direction
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%phase_fnc_solar(l, i, id, k)                        &
              =ss_prop%phase_fnc_solar(l, i, id, 0)
          END DO
        END DO
      END DO
    END IF
  END DO

! For use with McICA, save the clear-sky phase function (actually the
! product of the phase function and the scattering) before it is
! divided by the mean scatering.
  IF (control%i_cloud == ip_cloud_mcica .AND.                           &
      .NOT. control%l_avg_phase_fnc) THEN
    DO i=n_cloud_top, n_layer
      DO ls=1, n_order_phase
        DO l=1, n_profile
          ss_prop%phase_fnc_no_cloud(l, i, ls)                          &
            =ss_prop%phase_fnc(l, i, ls, 0)
        END DO
      END DO
    END DO

    IF (control%l_rescale) THEN
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%forward_scatter_no_cloud(l, i)                        &
            =ss_prop%forward_scatter(l, i, 0)
        END DO
      END DO
    END IF
  END IF


! Initialize arrays for diagnostic use.
  IF (control%l_cloud_extinction) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cloud_extinction(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (control%l_cloud_absorptivity) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cloud_absorptivity(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (control%l_ls_cloud_extinction) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           ls_cloud_extinction(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (control%l_ls_cloud_absorptivity) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           ls_cloud_absorptivity(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (control%l_cnv_cloud_extinction) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cnv_cloud_extinction(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (control%l_cnv_cloud_absorptivity) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cnv_cloud_absorptivity(l, i)=0.0_RealK
        END DO
     END DO
  END IF



! Add on the terms representing processes within clouds.

! Loop over the condensed components, calculating their optical
! properties and then assign them to the arrays for the types of
! cloud.

  DO k=1, n_condensed

!   Flags for dealing with components were set in the subroutine
!   set_cloud_pointer. we now determine whether the component is
!   to be included and calculate its optical properties according
!   to the phase of the component. these contributions are added
!   to the arrays for the selected type of cloud.

    IF (l_cloud_cmp(k)) THEN

      IF (i_phase_cmp(k) == ip_phase_water) THEN

!       Include scattering by water droplets.

        CALL opt_prop_water_cloud(ierr                                  &
          , n_profile, n_layer, n_cloud_top                             &
          , n_cloud_profile, i_cloud_profile                            &
          , n_order_phase, control%l_rescale, control%n_order_forward   &
          , control%l_henyey_greenstein_pf, l_solar_phf                 &
          , control%l_lanczos                                           &
          , control%n_order_phase_solar, n_direction, cos_sol_view      &
          , i_condensed_param(k)                                        &
          , condensed_param_list(:, k)                                  &
          , condensed_mix_ratio(:, :, k)                                &
          , condensed_dim_char(:, :, k)                                 &
          , p, density                                                  &
          , n_opt_level_drop_prsc, drop_pressure_prsc                   &
          , drop_absorption_prsc, drop_scattering_prsc                  &
          , drop_phase_fnc_prsc                                         &
          , ss_prop%k_ext_tot_cloud_comp(:, :, k)                       &
          , ss_prop%k_ext_scat_cloud_comp(:, :, k)                      &
          , ss_prop%phase_fnc_cloud_comp(:, :, :, k)                    &
          , ss_prop%forward_scatter_cloud_comp(:, :, k)                 &
          , ss_prop%forward_solar_cloud_comp(:, :, k)                   &
          , ss_prop%phase_fnc_solar_cloud_comp(:, :, :, k)              &
          , nd_profile, nd_radiance_profile, nd_layer, id_ct            &
          , nd_direction, nd_phf_term_cloud_prsc, nd_max_order          &
          , nd_cloud_parameter                                          &
          , nd_profile_cloud_prsc, nd_opt_level_cloud_prsc              &
          )

      ELSE IF (i_phase_cmp(k) == ip_phase_ice) THEN

!       Include scattering by ice crystals.

        CALL opt_prop_ice_cloud(ierr                                    &
          , n_profile, n_layer, n_cloud_top                             &
          , n_cloud_profile, i_cloud_profile                            &
          , n_order_phase, control%l_rescale, control%n_order_forward   &
          , control%l_henyey_greenstein_pf, l_solar_phf                 &
          , control%l_lanczos                                           &
          , control%n_order_phase_solar, n_direction, cos_sol_view      &
          , i_condensed_param(k)                                        &
          , condensed_param_list(:, k)                                  &
          , condensed_mix_ratio(:, :, k)                                &
          , condensed_dim_char(:, :, k)                                 &
          , p, t, density                                               &
          , n_opt_level_ice_prsc, ice_pressure_prsc                     &
          , ice_absorption_prsc, ice_scattering_prsc                    &
          , ice_phase_fnc_prsc                                          &
          , ss_prop%k_ext_tot_cloud_comp(:, :, k)                       &
          , ss_prop%k_ext_scat_cloud_comp(:, :, k)                      &
          , ss_prop%phase_fnc_cloud_comp(:, :, :, k)                    &
          , ss_prop%forward_scatter_cloud_comp(:, :, k)                 &
          , ss_prop%forward_solar_cloud_comp(:, :, k)                   &
          , ss_prop%phase_fnc_solar_cloud_comp(:, :, :, k)              &
          , nd_profile, nd_radiance_profile, nd_layer, id_ct            &
          , nd_direction                                                &
          , nd_phf_term_cloud_prsc, nd_max_order, nd_cloud_parameter    &
          , nd_profile_cloud_prsc, nd_opt_level_cloud_prsc              &
          )

      END IF

!     Apply cloud inhomogeneity correction

      IF (control%i_inhom == ip_cairns) THEN
        CALL opt_prop_inhom_corr_cairns(                                &
            n_layer, n_cloud_top                                        &
          , n_cloud_profile, i_cloud_profile                            &
          , control%l_rescale, control%n_order_forward                  &
          , condensed_rel_var_dens(:, :, k)                             &
          , ss_prop%k_ext_tot_cloud_comp(:, :, k)                       &
          , ss_prop%k_ext_scat_cloud_comp(:, :, k)                      &
          , ss_prop%phase_fnc_cloud_comp(:, :, :, k)                    &
          , ss_prop%forward_scatter_cloud_comp(:, :, k)                 &
          , nd_profile, nd_layer, id_ct                                 &
          )
      END IF

!     Increment the arrays of optical properties.

      IF (control%i_cloud /= ip_cloud_mcica .OR. control%l_avg_phase_fnc) THEN
      IF (control%l_rescale) THEN
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            ss_prop%k_grey_tot(l, i, i_cloud_type(k))                   &
              =ss_prop%k_grey_tot(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_tot_cloud_comp(l, i, k)
            ss_prop%k_ext_scat(l, i, i_cloud_type(k))                   &
              =ss_prop%k_ext_scat(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_scat_cloud_comp(l, i, k)
            DO ls=1, n_order_phase
              ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))              &
                =ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))           &
                +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)
            END DO
            ss_prop%forward_scatter(l, i, i_cloud_type(k))              &
              =ss_prop%forward_scatter(l, i, i_cloud_type(k))           &
              +ss_prop%forward_scatter_cloud_comp(l, i, k)
          END DO
        END DO
      ELSE
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            ss_prop%k_grey_tot(l, i, i_cloud_type(k))                   &
              =ss_prop%k_grey_tot(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_tot_cloud_comp(l, i, k)
            ss_prop%k_ext_scat(l, i, i_cloud_type(k))                   &
              =ss_prop%k_ext_scat(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_scat_cloud_comp(l, i, k)
            DO ls=1, n_order_phase
              ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))              &
                =ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))           &
                +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)
            END DO
          END DO
        END DO
      END IF
      IF (l_solar_phf) THEN
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            ss_prop%forward_solar(l, i, i_cloud_type(k))                &
              =ss_prop%forward_solar(l, i, i_cloud_type(k))             &
              +ss_prop%forward_solar_cloud_comp(l, i, k)
          END DO
        END DO
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO id=1, n_direction
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              ss_prop%phase_fnc_solar(l, i, id, i_cloud_type(k))        &
                =ss_prop%phase_fnc_solar(l, i, id, i_cloud_type(k))     &
                +ss_prop%phase_fnc_solar_cloud_comp(l, i, id, k)
            END DO
          END DO
        END DO
      END IF
      END IF


!     Extra calculations for diagnostics.

      IF (control%l_cloud_extinction) THEN
         DO i=n_cloud_top, n_layer
!CDIR NODEP
            DO ll=1, n_cloud_profile(i)
               l=i_cloud_profile(ll, i)
               cloud_extinction(l, i)                                   &
                  =cloud_extinction(l, i)                               &
                  +ss_prop%k_ext_tot_cloud_comp(l, i, k)                &
                  *frac_cloud(l, i, i_cloud_type(k))
            END DO
         END DO
      END IF


      IF (control%l_cloud_absorptivity) THEN
         DO i=n_cloud_top, n_layer
!CDIR NODEP
            DO ll=1, n_cloud_profile(i)
               l=i_cloud_profile(ll, i)
               cloud_absorptivity(l, i)                                 &
                  =cloud_absorptivity(l, i)                             &
                  +(ss_prop%k_ext_tot_cloud_comp(l, i, k)               &
                  -ss_prop%k_ext_scat_cloud_comp(l, i, k))              &
                  *frac_cloud(l, i, i_cloud_type(k))
            END DO
         END DO
      END IF

      IF ((i_cloud_type(k) == ip_cloud_type_sw).OR.                     &
          (i_cloud_type(k) == ip_cloud_type_si)) THEN

        IF (control%l_ls_cloud_extinction) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 ls_cloud_extinction(l, i)                              &
                    =ls_cloud_extinction(l, i)                          &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF


        IF (control%l_ls_cloud_absorptivity) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 ls_cloud_absorptivity(l, i)                            &
                    =ls_cloud_absorptivity(l, i)                        &
                    +(ss_prop%k_ext_tot_cloud_comp(l, i, k)             &
                    -ss_prop%k_ext_scat_cloud_comp(l, i, k))            &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF

      ELSE  !  Cloud is of convective type

        IF (control%l_cnv_cloud_extinction) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 cnv_cloud_extinction(l, i)                             &
                    =cnv_cloud_extinction(l, i)                         &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF


        IF (control%l_cnv_cloud_absorptivity) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 cnv_cloud_absorptivity(l, i)                           &
                    =cnv_cloud_absorptivity(l, i)                       &
                    +(ss_prop%k_ext_tot_cloud_comp(l, i, k)             &
                    -ss_prop%k_ext_scat_cloud_comp(l, i, k))            &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF
      END IF

    END IF

  END DO



! Calculate the final optical properties.
! The scattering was included in the free total extinction earlier,
! but we have yet to divide the product of the phase function and
! the scattering by the mean scattering.

  DO i=1, n_cloud_top-1

    n_index=0
    DO l=1, n_profile
      IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
        n_index=n_index+1
        indx(n_index)=l
      END IF
    END DO

    IF (control%l_rescale) THEN
!CDIR NODEP
      DO k=1, n_index
        tmp_inv(k)=1.0_RealK/ss_prop%k_ext_scat_clr(indx(k), i)
        ss_prop%forward_scatter_clr(indx(k), i)                         &
          =ss_prop%forward_scatter_clr(indx(k), i)*tmp_inv(k)
        DO ls=1, n_order_phase
          ss_prop%phase_fnc_clr(indx(k), i, ls)                         &
            =ss_prop%phase_fnc_clr(indx(k), i, ls)*tmp_inv(k)
        END DO
      END DO
      IF (control%i_direct_tau == ip_direct_csr_scaling ) THEN
        CALL circumsolar_fraction(n_index                               &
           , indx, control%half_angle                                   &
           , ss_prop%phase_fnc_clr(:, i, 1)                             &
           , ss_prop%forward_scatter_clr_csr(:, i)                      &
           , nd_profile                                                 &
           )
      END IF
    ELSE
      DO ls=1, n_order_phase
!CDIR NODEP
        DO k=1, n_index
          ss_prop%phase_fnc_clr(indx(k), i, ls)                         &
            =ss_prop%phase_fnc_clr(indx(k), i, ls)                      &
            /ss_prop%k_ext_scat_clr(indx(k), i)
        END DO
      END DO
    END IF

    IF (l_solar_phf) THEN
      DO k=1, n_index
        ss_prop%forward_solar_clr(indx(k), i)                           &
          =ss_prop%forward_solar_clr(indx(k), i)                        &
          /ss_prop%k_ext_scat_clr(indx(k), i)
      END DO
      DO id=1, n_direction
        DO k=1, n_index
          ss_prop%phase_fnc_solar_clr(indx(k), i, id)                   &
            =ss_prop%phase_fnc_solar_clr(indx(k), i, id)                &
            /ss_prop%k_ext_scat_clr(indx(k), i)
        END DO
      END DO
    END IF

  END DO

  DO i=n_cloud_top, n_layer

    n_index=0
    DO l=1, n_profile
      IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
        n_index=n_index+1
        indx(n_index)=l
      END IF
    END DO

    IF (control%l_rescale) THEN
!CDIR NODEP
      DO k=1, n_index
        tmp_inv(k)=1.0_RealK/ss_prop%k_ext_scat(indx(k), i, 0)
        ss_prop%forward_scatter(indx(k), i, 0)                          &
          =ss_prop%forward_scatter(indx(k), i, 0)*tmp_inv(k)
        DO ls=1, n_order_phase
          ss_prop%phase_fnc(indx(k), i, ls, 0)                          &
            =ss_prop%phase_fnc(indx(k), i, ls, 0)*tmp_inv(k)
        END DO
      END DO
      IF (control%i_direct_tau == ip_direct_csr_scaling ) THEN
        CALL circumsolar_fraction(n_index                               &
           , indx, control%half_angle                                   &
           , ss_prop%phase_fnc(:, i, 1, 0)                              &
           , ss_prop%forward_scatter_csr(:, i, 0)                       &
           , nd_profile                                                 &
           )
      END IF
    ELSE
      DO ls=1, n_order_phase
!CDIR NODEP
        DO k=1, n_index
          ss_prop%phase_fnc(indx(k), i, ls, 0)                          &
            =ss_prop%phase_fnc(indx(k), i, ls, 0)                       &
            /ss_prop%k_ext_scat(indx(k), i, 0)
        END DO
      END DO
    END IF

    IF (l_solar_phf) THEN
      DO k=1, n_index
        ss_prop%forward_solar(indx(k), i, 0)                            &
          =ss_prop%forward_solar(indx(k), i, 0)                         &
          /ss_prop%k_ext_scat(indx(k), i, 0)
      END DO
      DO id=1, n_direction
        DO k=1, n_index
          ss_prop%phase_fnc_solar(indx(k), i, id, 0)                    &
            =ss_prop%phase_fnc_solar(indx(k), i, id, 0)                 &
            /ss_prop%k_ext_scat(indx(k), i, 0)
        END DO
      END DO
    END IF

  END DO

  IF (control%i_cloud /= ip_cloud_mcica .OR. control%l_avg_phase_fnc) THEN

! Repeat for clouds.
  DO k=1, n_cloud_type
    DO i=n_cloud_top, n_layer
      n_index=0
      DO l=1, n_profile
        IF (ss_prop%k_ext_scat(l, i, k) > tiny_k) THEN
          n_index=n_index+1
          indx(n_index)=l
        END IF
      END DO

      IF (control%l_rescale) THEN
!CDIR NODEP
        DO j=1, n_index
          tmp_inv(j)=1.0_RealK/ss_prop%k_ext_scat(indx(j), i, k)
          ss_prop%forward_scatter(indx(j), i, k)                        &
            =ss_prop%forward_scatter(indx(j), i, k)*tmp_inv(j)
          DO ls=1, n_order_phase
            ss_prop%phase_fnc(indx(j), i, ls, k)                        &
              =ss_prop%phase_fnc(indx(j), i, ls, k)*tmp_inv(j)
          END DO
        END DO
      ELSE
        DO ls=1, n_order_phase
!CDIR NODEP
          DO j=1, n_index
            ss_prop%phase_fnc(indx(j), i, ls, k)                        &
              =ss_prop%phase_fnc(indx(j), i, ls, k)                     &
              /ss_prop%k_ext_scat(indx(j), i, k)
          END DO
        END DO
      END IF

      IF (l_solar_phf) THEN
        DO j=1, n_index
          ss_prop%forward_solar(indx(j), i, k)                          &
            =ss_prop%forward_solar(indx(j), i, k)                       &
            /ss_prop%k_ext_scat(indx(j), i, k)
        END DO
        DO id=1, n_direction
          DO j=1, n_index
            ss_prop%phase_fnc_solar(indx(j), i, id, k)                  &
              =ss_prop%phase_fnc_solar(indx(j), i, id, k)               &
              /ss_prop%k_ext_scat(indx(j), i, k)
          END DO
        END DO
      END IF
      IF ( control%l_rescale .AND.                                      &
        control%i_direct_tau == ip_direct_csr_scaling ) THEN
        CALL circumsolar_fraction(n_index                               &
           , indx, control%half_angle                                   &
           , ss_prop%phase_fnc(:, i, 1, k)                              &
           , ss_prop%forward_scatter_csr(:, i, k)                       &
           , nd_profile                                                 &
           )
      END IF

    END DO
  END DO

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE grey_opt_prop
END MODULE grey_opt_prop_mod

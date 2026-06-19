! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the optical properties of aerosols.
!
! Method:
!   If the optical properties come from an observational
!   distribution a separate subroutine is called. Otherwise
!   appropriate mean quantities in the layer are calculated
!   as the parametrization requires and these values are
!   substituted into the parametrization to give the optical
!   properties. Aerosol properties may depend on the humidity.
!
!- ---------------------------------------------------------------------
MODULE opt_prop_aerosol_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OPT_PROP_AEROSOL_MOD'
CONTAINS
SUBROUTINE opt_prop_aerosol(ierr                                        &
    , n_profile, first_layer, last_layer                                &
    , n_order_phase, l_rescale, n_order_forward                         &
    , l_henyey_greenstein_pf                                            &
    , n_aerosol_mr, aerosol_mix_ratio                                   &
    , aerosol_mr_source, aerosol_mr_type_index                          &
    , i_aerosol_parametrization                                         &
    , i_humidity_pointer, humidities, delta_humidity                    &
    , mean_rel_humidity                                                 &
    , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc         &
    , l_solar_phf, n_order_phase_solar, n_direction, cos_sol_view       &
    , p, density                                                        &
    , n_opt_level_aerosol_prsc, aerosol_pressure_prsc                   &
    , aerosol_absorption_prsc, aerosol_scattering_prsc                  &
    , aerosol_phase_fnc_prsc                                            &
    , k_ext_tot, k_ext_scat, phase_fnc, forward_scatter                 &
    , forward_solar, phase_fnc_solar                                    &
    , nd_profile, nd_radiance_profile, nd_layer, id_lt, id_lb           &
    , nd_aerosol_species, nd_aerosol_mixratio, nd_humidities            &
    , nd_phase_term, nd_max_order, nd_direction                         &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    , nd_phf_term_prsc                                                  &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_aerosol_param_dry,                 &
                     ip_aerosol_param_moist, ip_aerosol_param_phf_dry,  &
                     ip_aerosol_param_phf_moist,                        &
                     ip_aerosol_unparametrized, ip_aersrc_arcl_ron,     &
                     ip_aersrc_classic_ron, ip_aersrc_cusack_ron
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE prsc_opt_prop_mod, ONLY: prsc_opt_prop

  IMPLICIT NONE


  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_radiance_profile                                               &
!       Size allocated for profiles of quantities only
!       required when radiances are wanted
    , nd_layer                                                          &
!       Size allocated for layers
    , id_lt                                                             &
!       Topmost declared layer for output optical properties
    , id_lb                                                             &
!       Bottom declared layer for output optical properties
    , nd_phase_term                                                     &
!       Size allocated for terms in phase function
    , nd_max_order                                                      &
!       Size allocated for orders of sperical harmonics
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_aerosol_species                                                &
!       Size allocated for aerosol species in spectral information
    , nd_aerosol_mixratio                                               &
!       Size allocated for aerosols in aerosol_mix_ratio
    , nd_humidities                                                     &
!       Size allocated for humidities
    , nd_profile_prsc                                                   &
!       Size allowed for profiles of prescribed properties
    , nd_opt_level_prsc                                                 &
!       Size allowed for levels of prescribed properties
    , nd_phf_term_prsc
!       Size allowed for terms in phase function prescribed properties

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , first_layer                                                       &
!       First layer where propeties are required
    , last_layer                                                        &
!       Last layer where propeties are required
    , n_order_phase                                                     &
!       Number of terms to retain in the phase function
    , n_order_phase_solar                                               &
!       Number of terms to retain in calculating the angular
!       scattering of solar radiation
    , n_direction                                                       &
!       Number of viewing directions
    , n_order_forward
!       Order used in forming the forward scattering parameter

  LOGICAL, INTENT(IN) ::                                                &
      l_rescale                                                         &
!       Flag for delta-rescaling
    , l_henyey_greenstein_pf                                            &
!       Flag to use a Henyey-Greenstein phase function
    , l_solar_phf
!       Flag to use calculate a separate solar phase function

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
      cos_sol_view(nd_radiance_profile, nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing direction

  REAL (RealK), INTENT(IN) ::                                           &
      aerosol_mix_ratio(nd_profile, nd_layer                            &
        , nd_aerosol_mixratio)                                          &
!       Number densty of aerosols
    , aerosol_absorption(nd_humidities, nd_aerosol_species)             &
!       Aerosol absorption in band/mix ratio
    , aerosol_scattering(nd_humidities, nd_aerosol_species)             &
!       Aerosol scattering in band/mix ratio
    , aerosol_phase_fnc(nd_humidities                                   &
        , nd_phase_term, nd_aerosol_species)                            &
!       Aerosol phase function in band
    , humidities(nd_humidities, nd_aerosol_species)                     &
!       Array of humidities
    , delta_humidity                                                    &
!       Increment in humidity
    , mean_rel_humidity(nd_profile, nd_layer)
!       Mixing ratio of water vapour

  INTEGER, INTENT(IN) ::                                                &
      n_opt_level_aerosol_prsc(nd_aerosol_species)
!       Number of aerosol data layers

  REAL (RealK), INTENT(IN) ::                                           &
      p(nd_profile, nd_layer)                                           &
!       Pressure
    , density(nd_profile, nd_layer)                                     &
!       Atmospheric density
    , aerosol_pressure_prsc(nd_profile_prsc, nd_opt_level_prsc          &
        , nd_aerosol_species)                                           &
!       Pressures at levels of prescribed aerosol properties
    , aerosol_absorption_prsc(nd_profile_prsc, nd_opt_level_prsc        &
        , nd_aerosol_species)                                           &
!       Prescribed aerosol absorption
    , aerosol_scattering_prsc(nd_profile_prsc, nd_opt_level_prsc        &
        , nd_aerosol_species)                                           &
!       Prescribed aerosol scattering
    , aerosol_phase_fnc_prsc(nd_profile_prsc, nd_opt_level_prsc         &
        , nd_phf_term_prsc, nd_aerosol_species)
!       Prescribed aerosol phase function

  REAL (RealK), INTENT(INOUT) ::                                        &
      k_ext_scat(nd_profile, id_lt: id_lb)                              &
!       Scattering extinction
    , k_ext_tot(nd_profile, id_lt: id_lb)                               &
!       Total extinction
    , phase_fnc(nd_profile, id_lt: id_lb, nd_max_order)                 &
!       Phase function
    , forward_scatter(nd_profile, id_lt: id_lb)                         &
!       Forward scattering
    , forward_solar(nd_profile, id_lt: id_lb)                           &
!       Forward scattering for the solar beam
    , phase_fnc_solar(nd_radiance_profile                               &
        , id_lt: id_lb, nd_direction)
!       Phase function relative to the solar beam

! Local variables.
  LOGICAL                                                               &
      l_use_hg_phf                                                      &
!       Flag to use Henyey-Greenstein phase functions
    , l_interpolate_hum
!       Flag to interpolate optical properties through a look-up
!       table of humidities
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , j_mr                                                              &
!       Loop variable over aerosol mixing rations
    , j                                                                 &
!       j_mr, converted to aerosol spectral information
    , i                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , ls                                                                &
!       Loop variable
    , i_pointer
!       Temporary pointer
  REAL (RealK) ::                                                       &
      k_scatter(nd_profile)                                             &
!       Scattering of current extinction of the current aerosol
    , asymmetry(nd_profile)                                             &
!       Asymmetry of the current aerosol
    , ks_phf(nd_profile)                                                &
!       Scattering coefficient multiplied by a coefficient in the
!       phase function
    , phf_coeff(nd_profile)
!       Coefficient in the phase function of the current aerosol

  REAL (RealK) ::                                                       &
      weight_upper(nd_profile)
!       Weighting towards the upper end of an interval
!       in a look-up table

! Optical properties interpolated from prescribed properties
! for each component
  REAL (RealK) ::                                                       &
      k_ext_scat_comp(nd_profile, id_lt: id_lb)                         &
!       Scattering extinction of component
    , k_ext_tot_comp(nd_profile, id_lt: id_lb)                          &
!       Total extinction of component
    , phase_fnc_comp(nd_profile, id_lt: id_lb, nd_max_order)            &
!       Phase function of component
    , phase_fnc_solar_comp(nd_radiance_profile                          &
        , id_lt: id_lb, nd_direction)                                   &
!       Solar phase function of component
    , forward_scatter_comp(nd_profile, id_lt: id_lb)                    &
!       Forward scattering of component
    , forward_solar_comp(nd_profile, id_lt: id_lb)
!       Forward scattering of the solar beam for the component

! Legendre polynomials:
  REAL (RealK) ::                                                       &
      cnst1                                                             &
!       Constant in recurrence for Legendre polynomials
    , p_legendre_ls(nd_radiance_profile)                                &
!       Legendre polynomial at the current order
    , p_legendre_ls_m1(nd_radiance_profile)                             &
!       Legendre polynomial at the previous order
    , p_legendre_tmp(nd_radiance_profile)                               &
!       Temporary Legendre polynomial
    , phase_fnc_solar_tmp(nd_radiance_profile)
!       Current contribution to the solar phase function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'OPT_PROP_AEROSOL'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO j_mr=1, n_aerosol_mr

!  Test to see if this aerosol is radiatively active:
   IF (aerosol_mr_source(j_mr) == ip_aersrc_cusack_ron .OR.            &
       aerosol_mr_source(j_mr) == ip_aersrc_classic_ron  .OR.          &
       aerosol_mr_source(j_mr) == ip_aersrc_arcl_ron ) THEN

!   Get the J relevant for the spectral information
    j=aerosol_mr_type_index(j_mr)

!   Optimise two stream cases most often used.
    IF ( (n_order_phase == 1) .AND.                                     &
      (i_aerosol_parametrization(j) == ip_aerosol_param_dry) .AND.      &
      l_rescale .AND. (n_order_forward == 2) ) THEN

      DO i=first_layer, last_layer
!CDIR NODEP
        DO l=1, n_profile
          k_ext_tot(l, i)=k_ext_tot(l, i)                               &
            +aerosol_mix_ratio(l, i, j_mr)                              &
            *aerosol_absorption(1, j)
          k_scatter(l)=aerosol_mix_ratio(l, i, j_mr)                    &
            *aerosol_scattering(1, j)
          k_ext_scat(l, i)=k_ext_scat(l, i)                             &
            +k_scatter(l)
          ks_phf(l)=k_scatter(l)*aerosol_phase_fnc(1, 1, j)
          phase_fnc(l, i, 1)=phase_fnc(l, i, 1)+ks_phf(l)
          forward_scatter(l, i)=forward_scatter(l, i)                   &
            +ks_phf(l)*aerosol_phase_fnc(1, 1, j)
        END DO
      END DO

    ELSE IF ( (n_order_phase == 1) .AND.                                &
      (i_aerosol_parametrization(j) == ip_aerosol_param_dry) .AND.      &
      .NOT. l_rescale ) THEN

      DO i=first_layer, last_layer
!CDIR NODEP
        DO l=1, n_profile
          k_ext_tot(l, i)=k_ext_tot(l, i)                               &
            +aerosol_mix_ratio(l, i, j_mr)                              &
            *aerosol_absorption(1, j)
          k_scatter(l)=aerosol_mix_ratio(l, i, j_mr)                    &
            *aerosol_scattering(1, j)
          k_ext_scat(l, i)=k_ext_scat(l, i)                             &
            +k_scatter(l)
          phase_fnc(l, i, 1)=phase_fnc(l, i, 1)                         &
            +k_scatter(l)*aerosol_phase_fnc(1, 1, j)
        END DO
      END DO

    ELSE IF ( (n_order_phase == 1) .AND.                                &
      (i_aerosol_parametrization(j) == ip_aerosol_param_moist) .AND.    &
      l_rescale .AND. (n_order_forward == 2) ) THEN

      DO i=first_layer, last_layer
!CDIR NODEP
        DO l=1, n_profile
          i_pointer=i_humidity_pointer(l, i)
          weight_upper(l)=(mean_rel_humidity(l, i)                      &
            -humidities(i_pointer, j))                                  &
            /delta_humidity
          k_ext_tot(l, i)=k_ext_tot(l, i)                               &
            +aerosol_mix_ratio(l, i, j_mr)                              &
            *(aerosol_absorption(i_pointer, j)                          &
            +weight_upper(l)                                            &
            *(aerosol_absorption(i_pointer+1, j)                        &
            -aerosol_absorption(i_pointer, j)))
          k_scatter(l)=aerosol_mix_ratio(l, i, j_mr)                    &
            *(aerosol_scattering(i_pointer, j)                          &
            +weight_upper(l)                                            &
            *(aerosol_scattering(i_pointer+1, j)                        &
            -aerosol_scattering(i_pointer, j)))
          k_ext_scat(l, i)=k_ext_scat(l, i)                             &
            +k_scatter(l)
          asymmetry(l)                                                  &
            =aerosol_phase_fnc(i_pointer, 1, j)                         &
            +weight_upper(l)                                            &
            *(aerosol_phase_fnc(i_pointer+1, 1, j)                      &
            -aerosol_phase_fnc(i_pointer, 1, j))
          ks_phf(l)=k_scatter(l)*asymmetry(l)
          phase_fnc(l, i, 1)=phase_fnc(l, i, 1)+ks_phf(l)
          forward_scatter(l, i)=forward_scatter(l, i)                   &
            +ks_phf(l)*asymmetry(l)
        END DO
      END DO

    ELSE IF ( (n_order_phase == 1) .AND.                                &
      (i_aerosol_parametrization(j) == ip_aerosol_param_moist) .AND.    &
      .NOT. l_rescale ) THEN

      DO i=first_layer, last_layer
!CDIR NODEP
        DO l=1, n_profile
          i_pointer=i_humidity_pointer(l, i)
          weight_upper(l)=(mean_rel_humidity(l, i)                      &
            -humidities(i_pointer, j))                                  &
            /delta_humidity
          k_ext_tot(l, i)=k_ext_tot(l, i)                               &
            +aerosol_mix_ratio(l, i, j_mr)                              &
            *(aerosol_absorption(i_pointer, j)                          &
            +weight_upper(l)                                            &
            *(aerosol_absorption(i_pointer+1, j)                        &
            -aerosol_absorption(i_pointer, j)))
          k_scatter(l)=aerosol_mix_ratio(l, i, j_mr)                    &
            *(aerosol_scattering(i_pointer, j)                          &
            +weight_upper(l)                                            &
            *(aerosol_scattering(i_pointer+1, j)                        &
            -aerosol_scattering(i_pointer, j)))
          k_ext_scat(l, i)=k_ext_scat(l, i)                             &
            +k_scatter(l)
          phase_fnc(l, i, 1)=phase_fnc(l, i, 1)                         &
            +k_scatter(l)*(aerosol_phase_fnc(i_pointer, 1, j)           &
            +weight_upper(l)                                            &
            *(aerosol_phase_fnc(i_pointer+1, 1, j)                      &
            -aerosol_phase_fnc(i_pointer, 1, j)))
        END DO
      END DO

!   Now treat all other cases in a general way.
    ELSE IF ( (i_aerosol_parametrization(j) ==                          &
          ip_aerosol_param_dry).OR.                                     &
         (i_aerosol_parametrization(j) ==                               &
          ip_aerosol_param_phf_dry).OR.                                 &
         (i_aerosol_parametrization(j) ==                               &
          ip_aerosol_param_moist).OR.                                   &
         (i_aerosol_parametrization(j) ==                               &
          ip_aerosol_param_phf_moist) ) THEN

!     Use the Henyey-Greenstein phase function if specifically
!     requested, or if using an old parametrization which gives
!     only an asymmetry.
      l_use_hg_phf=l_henyey_greenstein_pf.OR.                           &
        (i_aerosol_parametrization(j) == ip_aerosol_param_dry).OR.      &
        (i_aerosol_parametrization(j) == ip_aerosol_param_moist)

!     Interpolate from the look-up table if using moist properties.
      l_interpolate_hum=                                                &
        (i_aerosol_parametrization(j) == ip_aerosol_param_moist).OR.    &
        (i_aerosol_parametrization(j) == ip_aerosol_param_phf_moist)

      DO i=first_layer, last_layer
        IF (l_interpolate_hum) THEN

!         Calculate the required weights for interpolation
!         in this layer.
          DO l=1, n_profile
            i_pointer=i_humidity_pointer(l, i)
            weight_upper(l)=(mean_rel_humidity(l, i)                    &
              -humidities(i_pointer, j))                                &
              /delta_humidity
          END DO

!         Interpolate the absorption and scattering.
          DO l=1, n_profile
            i_pointer=i_humidity_pointer(l, i)
            k_ext_tot(l, i)=k_ext_tot(l, i)                             &
              +aerosol_mix_ratio(l, i, j_mr)                            &
              *(aerosol_absorption(i_pointer, j)                        &
              +weight_upper(l)                                          &
              *(aerosol_absorption(i_pointer+1, j)                      &
              -aerosol_absorption(i_pointer, j)))
            k_scatter(l)=aerosol_mix_ratio(l, i, j_mr)                  &
              *(aerosol_scattering(i_pointer, j)                        &
              +weight_upper(l)                                          &
              *(aerosol_scattering(i_pointer+1, j)                      &
              -aerosol_scattering(i_pointer, j)))
            k_ext_scat(l, i)=k_ext_scat(l, i)                           &
              +k_scatter(l)
          END DO

        ELSE

          DO l=1, n_profile

!           Calculate volume extinctions directly from the
!           mass extinction coefficients.
            k_ext_tot(l, i)=k_ext_tot(l, i)                             &
              +aerosol_mix_ratio(l, i, j_mr)                            &
              *aerosol_absorption(1, j)
            k_scatter(l)=aerosol_mix_ratio(l, i, j_mr)                  &
              *aerosol_scattering(1, j)
            k_ext_scat(l, i)=k_ext_scat(l, i)                           &
              +k_scatter(l)

          END DO

        END IF

!       The phase function:

        IF (l_use_hg_phf) THEN

!         Note that there is an ambiguity in the definition of a
!         Henyey-Greenstein phase function when humidity is included
!         since one could set up the lookup table externally with
!         all moments at the reference points set to powers of the
!         appropriate asymmetries, but then linear interpolation in
!         the humidity would not give a true Henyey-Greenstein
!         phase function at intermediate points. Here we adopt a
!         true Henyey-Greenstein approach, calculating just the
!         asymmetry.

!         Calculate the asymmetry:
          IF (l_interpolate_hum) THEN
            DO l=1, n_profile
              i_pointer=i_humidity_pointer(l, i)
              asymmetry(l)                                              &
                =aerosol_phase_fnc(i_pointer, 1, j)                     &
                +weight_upper(l)                                        &
                *(aerosol_phase_fnc(i_pointer+1, 1, j)                  &
                -aerosol_phase_fnc(i_pointer, 1, j))
            END DO
          ELSE
            DO l=1, n_profile
              asymmetry(l)=aerosol_phase_fnc(1, 1, j)
            END DO
          END IF

!         Set the lowest order in the phase function (required
!         for two-stream calculations and other quadratures).
          DO l=1, n_profile
            phase_fnc(l, i, 1)=phase_fnc(l, i, 1)                       &
              +k_scatter(l)*asymmetry(l)
          END DO

!         Initialize the product of the scattering and the
!         current moment of the phase function. This repeats
!         part of the preceeding loop, but separating it saves
!         an assignment in the case of two-stream calculations.
          IF (l_rescale.OR.(n_order_phase >= 2)) THEN
            DO l=1, n_profile
              ks_phf(l)=k_scatter(l)*asymmetry(l)
            END DO
          END IF

!         Calculate weighted higher moments recursively.
          DO ls=2, n_order_phase
            DO l=1, n_profile
              ks_phf(l)=ks_phf(l)*asymmetry(l)
              phase_fnc(l, i, ls)                                       &
                =phase_fnc(l, i, ls)+ks_phf(l)
            END DO
          END DO

!         Usually, we will retain terms as far as the order of
!         truncation, but in the case of two-stream methods the
!         order of truncation will exceed the order of retention
!         by 1.
          IF (l_rescale) THEN

            IF (n_order_forward == n_order_phase) THEN
              DO l=1, n_profile
                forward_scatter(l, i)                                   &
                  =forward_scatter(l, i)+ks_phf(l)
              END DO
            ELSE IF (n_order_forward == n_order_phase+1) THEN
              DO l=1, n_profile
                forward_scatter(l, i)                                   &
                  =forward_scatter(l, i)+ks_phf(l)*asymmetry(l)
              END DO
            ELSE
!             This case probably shouldn't arise so we use
!             inefficient explicit exponentiation.
              DO l=1, n_profile
                forward_scatter(l, i)                                   &
                  =forward_scatter(l, i)                                &
                  +k_scatter(l)*asymmetry(l)**n_order_forward
              END DO
            END IF

          END IF


        ELSE

!         Calculate the phase function generally. We don't
!         separate the first order here, because it is unlikely
!         that this block will be used in the case in a
!         two-stream calculation.
          DO ls=1, n_order_phase
            IF (l_interpolate_hum) THEN
              DO l=1, n_profile
                i_pointer=i_humidity_pointer(l, i)
                phf_coeff(l)=aerosol_phase_fnc(i_pointer, ls, j)        &
                  +weight_upper(l)                                      &
                  *(aerosol_phase_fnc(i_pointer+1, ls, j)               &
                  -aerosol_phase_fnc(i_pointer, ls, j))
              END DO
            ELSE
              DO l=1, n_profile
                phf_coeff(l)=aerosol_phase_fnc(1, ls, j)
              END DO
            END IF
            DO l=1, n_profile
              phase_fnc(l, i, ls)=phase_fnc(l, i, ls)                   &
                +k_scatter(l)*phf_coeff(l)
            END DO
          END DO

          IF (l_rescale) THEN
            IF (l_interpolate_hum) THEN
              DO l=1, n_profile
                i_pointer=i_humidity_pointer(l, i)
                phf_coeff(l)                                            &
                  =aerosol_phase_fnc(i_pointer, n_order_forward, j)     &
                  +weight_upper(l)                                      &
                  *(aerosol_phase_fnc(i_pointer+1                       &
                  , n_order_forward, j)                                 &
                  -aerosol_phase_fnc(i_pointer                          &
                  , n_order_forward, j))
              END DO
            ELSE
              DO l=1, n_profile
                phf_coeff(l)                                            &
                  =aerosol_phase_fnc(1, n_order_forward, j)
              END DO
            END IF
            DO l=1, n_profile
              forward_scatter(l, i)                                     &
                =forward_scatter(l, i)                                  &
                +k_scatter(l)*phf_coeff(l)
            END DO
          END IF
        END IF

        IF (l_solar_phf) THEN
!         Calculate the solar phase function to higher accuracy.
          DO id=1, n_direction
!           The Legendre polynomials are not stored so as to reduce
!           the requirement for memory at very high orders of solar
!           truncation.
            IF (l_interpolate_hum) THEN
              DO l=1, n_profile
                i_pointer=i_humidity_pointer(l, i)
                phf_coeff(l)=aerosol_phase_fnc(i_pointer, 1, j)         &
                  +weight_upper(l)                                      &
                  *(aerosol_phase_fnc(i_pointer+1, 1, j)                &
                  -aerosol_phase_fnc(i_pointer, 1, j))
              END DO
            ELSE
              DO l=1, n_profile
                phf_coeff(l)=aerosol_phase_fnc(1, 1, j)
              END DO
            END IF
            DO l=1, n_profile
!             Initialize the Legendre polynomials at the zeroth and
!             first orders.
              p_legendre_ls_m1(l)=1.0e+00_RealK
              p_legendre_ls(l)=cos_sol_view(l, id)
              phase_fnc_solar_tmp(l)=1.0e+00_RealK+phf_coeff(l)         &
                *p_legendre_ls(l)*REAL(2*1+1, RealK)
            END DO

            IF (l_use_hg_phf) THEN
!             Store the asymmetry in this case.
              DO l=1, n_profile
                asymmetry(l)=phf_coeff(l)
              END DO
            END IF

            DO ls=2, n_order_phase_solar
!             Calculate higher orders by recurrences.
              cnst1=1.0e+00_RealK-1.0e+00_RealK/REAL(ls, RealK)
              DO l=1, n_profile
                p_legendre_tmp(l)=p_legendre_ls(l)
                p_legendre_ls(l)                                        &
                  =(1.0e+00_RealK+cnst1)*p_legendre_ls(l)               &
                  *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
                p_legendre_ls_m1(l)=p_legendre_tmp(l)
              END DO

!             Calculate the next moment of the phase function.
              IF (l_use_hg_phf) THEN
                DO l=1, n_profile
                  phf_coeff(l)=phf_coeff(l)*asymmetry(l)
                END DO
              ELSE
                IF (l_interpolate_hum) THEN
                  DO l=1, n_profile
                    i_pointer=i_humidity_pointer(l, i)
                    phf_coeff(l)                                        &
                      =aerosol_phase_fnc(i_pointer, ls, j)              &
                      +weight_upper(l)                                  &
                      *(aerosol_phase_fnc(i_pointer+1, ls, j)           &
                      -aerosol_phase_fnc(i_pointer, ls, j))
                  END DO
                ELSE
                  DO l=1, n_profile
                    phf_coeff(l)=aerosol_phase_fnc(1, ls, j)
                  END DO
                END IF
              END IF
              DO l=1, n_profile
                phase_fnc_solar_tmp(l)=phase_fnc_solar_tmp(l)           &
                  +phf_coeff(l)                                         &
                  *REAL(2*ls+1, RealK)*p_legendre_ls(l)
              END DO
            END DO
!           Increment the stored phase function.
            DO l=1, n_profile
              phase_fnc_solar(l, i, id)                                 &
                =phase_fnc_solar(l, i, id)                              &
              +k_scatter(l)*phase_fnc_solar_tmp(l)
            END DO
          END DO

!         Continue to an extra order to find the rescaling
!         for the solar beam.
          IF (l_rescale) THEN
            IF (l_use_hg_phf) THEN
              DO l=1, n_profile
                phf_coeff(l)=phf_coeff(l)*asymmetry(l)
              END DO
            ELSE
              ls=n_order_phase_solar+1
              IF (l_interpolate_hum) THEN
                DO l=1, n_profile
                  i_pointer=i_humidity_pointer(l, i)
                  phf_coeff(l)                                          &
                    =aerosol_phase_fnc(i_pointer, ls, j)                &
                    +weight_upper(l)                                    &
                    *(aerosol_phase_fnc(i_pointer+1, ls, j)             &
                    -aerosol_phase_fnc(i_pointer, ls, j))
                END DO
              ELSE
                DO l=1, n_profile
                  phf_coeff(l)=aerosol_phase_fnc(1, ls, j)
                END DO
              END IF
            END IF
            DO l=1, n_profile
              forward_solar(l, i)=forward_solar(l, i)                   &
                +k_scatter(l)*phf_coeff(l)
            END DO
          END IF

        END IF

      END DO


    ELSE IF (i_aerosol_parametrization(j) ==                            &
             ip_aerosol_unparametrized) THEN
       CALL prsc_opt_prop(ierr                                          &
         , n_profile, first_layer, last_layer                           &
         , l_rescale, n_order_forward                                   &
         , l_henyey_greenstein_pf, n_order_phase                        &
         , p, density                                                   &
         , n_opt_level_aerosol_prsc(j)                                  &
         , aerosol_pressure_prsc(1, 1, j)                               &
         , aerosol_absorption_prsc(1, 1, j)                             &
         , aerosol_scattering_prsc(1, 1, j)                             &
         , aerosol_phase_fnc_prsc(1, 1, 1, j)                           &
         , k_ext_tot_comp, k_ext_scat_comp, phase_fnc_comp              &
         , forward_scatter_comp, forward_solar_comp                     &
         , l_solar_phf, n_order_phase_solar                             &
         , n_direction, cos_sol_view                                    &
         , phase_fnc_solar_comp                                         &
         , nd_profile, nd_radiance_profile, nd_layer                    &
         , id_lt, id_lb                                                 &
         , nd_direction                                                 &
         , nd_profile_prsc, nd_opt_level_prsc                           &
         , nd_phf_term_prsc, nd_max_order                               &
         )

      DO i=first_layer, last_layer
        DO l=1, n_profile
          k_ext_tot(l, i)=k_ext_tot(l, i)                               &
            +k_ext_tot_comp(l, i)
          k_ext_scat(l, i)=k_ext_scat(l, i)                             &
            +k_ext_scat_comp(l, i)
        END DO
      END DO
      DO ls=1, n_order_phase
        DO i=first_layer, last_layer
          DO l=1, n_profile
            phase_fnc(l, i, ls)=phase_fnc(l, i, ls)                     &
              +phase_fnc_comp(l, i, ls)
          END DO
        END DO
      END DO
      IF (l_rescale) THEN
        DO i=first_layer, last_layer
          DO l=1, n_profile
            forward_scatter(l, i)=forward_scatter(l, i)                 &
              +forward_scatter_comp(l, i)
          END DO
        END DO
      END IF
      IF (l_solar_phf) THEN
        IF (l_rescale) THEN
          DO i=first_layer, last_layer
            DO l=1, n_profile
              forward_solar(l, i)=forward_solar(l, i)                   &
                +forward_solar_comp(l, i)
            END DO
          END DO
        END IF
        DO id=1, n_direction
          DO i=first_layer, last_layer
            DO l=1, n_profile
              phase_fnc_solar(l, i, id)                                 &
                =phase_fnc_solar(l, i, id)                              &
                +phase_fnc_solar_comp(l, i, id)
            END DO
          END DO
        END DO
      END IF

    ELSE

      WRITE(cmessage, '(a, i3, a)')                                     &
        'i_aerosol_parametrization for species ', j,                    &
        ' has been set to an illegal value.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)

    END IF

   END IF ! Ends test for radiatively active aerosol

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE opt_prop_aerosol
END MODULE opt_prop_aerosol_mod

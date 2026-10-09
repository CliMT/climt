! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate optical properties of ice clouds.
!
! Method:
!   If the optical properties come from an observational
!   distribution a separate subroutine is called. Otherwise
!   appropriate mean quantities in the layer are calculated
!   as the parametrization requires and these values are
!   substituted into the parametrization to give the optical
!   properties.
!
!- ---------------------------------------------------------------------
MODULE opt_prop_ice_cloud_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OPT_PROP_ICE_CLOUD_MOD'
CONTAINS
SUBROUTINE opt_prop_ice_cloud(ierr                                      &
    , n_profile, n_layer, n_cloud_top                                   &
    , n_cloud_profile, i_cloud_profile                                  &
    , n_order_phase, l_rescale, n_order_forward                         &
    , l_henyey_greenstein_pf, l_solar_phf, l_lanczos                    &
    , n_order_phase_solar, n_direction, cos_sol_view                    &
    , i_parametrization_ice, ice_cloud_parameter                        &
    , ice_mass_frac, dim_char_ice                                       &
    , p, t, density                                                     &
    , n_opt_level_cloud_prsc, ice_pressure_prsc                         &
    , ice_absorption_prsc, ice_scattering_prsc                          &
    , ice_phase_fnc_prsc                                                &
    , k_ext_tot_cloud, k_ext_scat_cloud                                 &
    , phase_fnc_cloud, forward_scatter_cloud                            &
    , forward_solar_cloud, phase_fnc_solar_cloud                        &
    , nd_profile, nd_radiance_profile, nd_layer, id_ct                  &
    , nd_direction, nd_phase_term, nd_max_order                         &
    , nd_cloud_parameter                                                &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_ice_fu_phf, ip_ice_baran,          &
                     ip_aerosol_param_dry, ip_aerosol_param_moist,      &
                     ip_aerosol_param_phf_dry, ip_aersrc_classic_ron,   &
                     ip_aerosol_param_phf_moist, ip_aersrc_cusack_ron,  &
                     ip_aerosol_unparametrized, ip_aersrc_arcl_ron,     &
                     ip_sun_shine_vn2_ir, ip_slingo_schrecker_ice,      &
                     ip_slingo_schr_ice_phf, ip_ice_unparametrized,     &
                     ip_ice_t_iwc, ip_ice_iwc_only, ip_ice_fu_solar,    &
                     ip_ice_fu_ir, ip_ice_adt_10, ip_ice_adt,           &
                     ip_sun_shine_vn2_vis, ip_ice_pade_2_phf
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE opt_prop_pade_2_mod, ONLY: opt_prop_pade_2
  USE opt_prop_fu_phf_mod, ONLY: opt_prop_fu_phf
  USE opt_prop_baran_mod, ONLY: opt_prop_baran
  USE prsc_opt_prop_mod, ONLY: prsc_opt_prop

  IMPLICIT NONE


  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_radiance_profile                                               &
!       Size allocated for points where radiances are calculated
    , nd_layer                                                          &
!       Size allocated for layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_cloud_parameter                                                &
!       Size allocated for cloud parameters
    , nd_phase_term                                                     &
!       Size allocated for terms in the phase function
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_profile_prsc                                                   &
!       Size allowed for profiles of prescribed optical properties
    , nd_opt_level_prsc
!       Size allowed for levels of prescribed optical properties

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_order_phase                                                     &
!       Order of the phase function
    , n_order_phase_solar                                               &
!       Number of terms to retain in single scattered solar
!       phase function
    , n_order_forward                                                   &
!       Order used in forming the forward scattering parameter
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , i_parametrization_ice                                             &
!       Treatment of ice crystals
    , n_opt_level_cloud_prsc                                            &
!       Number of levels of prescribed data
    , n_cloud_profile(id_ct: nd_layer)                                  &
!       Number of cloudy profiles
    , i_cloud_profile(nd_profile, id_ct: nd_layer)
!       Profiles containing clouds
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale                                                         &
!       Delta-rescaling required
    , l_henyey_greenstein_pf                                            &
!       Flag to use a Henyey-Greenstein phase function
    , l_solar_phf                                                       &
!       Flag to use an extended solar phase function in
!       single scattering
    , l_lanczos
!       Flag to use Lanczos smoothing of solar phf

! Viewing directions:
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing dierctions
  REAL (RealK), INTENT(IN) ::                                           &
      cos_sol_view(nd_radiance_profile, nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing direction

  REAL (RealK), INTENT(IN) ::                                           &
      ice_cloud_parameter(nd_cloud_parameter)                           &
!       Ice cloud parameters
    , ice_mass_frac(nd_profile, id_ct: nd_layer)                        &
!       Ice mass fraction
    , dim_char_ice(nd_profile, id_ct: nd_layer)                         &
!       Characteristic dimension for crystals
    , p(nd_profile, nd_layer)                                           &
!       Pressure
    , t(nd_profile, nd_layer)                                           &
!       Temperature
    , density(nd_profile, nd_layer)                                     &
!       Atmospheric density
    , ice_pressure_prsc(nd_profile_prsc, nd_opt_level_prsc)             &
!       Pressure at which optical properties are prescribed
    , ice_absorption_prsc(nd_profile_prsc, nd_opt_level_prsc)           &
!       Prescribed absorption by ice crystals
    , ice_scattering_prsc(nd_profile_prsc, nd_opt_level_prsc)           &
!       Prescribed absorption by ice crystals
    , ice_phase_fnc_prsc(nd_profile_prsc, nd_opt_level_prsc             &
        , nd_phase_term)
!       Prescribed phase functions of ice crystals

  REAL (RealK), INTENT(OUT) ::                                          &
      k_ext_scat_cloud(nd_profile, id_ct: nd_layer)                     &
!       Scattering extinction
    , k_ext_tot_cloud(nd_profile, id_ct: nd_layer)                      &
!       Total extinction
    , phase_fnc_cloud(nd_profile, id_ct: nd_layer, nd_max_order)        &
!       Cloudy phase function
    , phase_fnc_solar_cloud(nd_radiance_profile, id_ct: nd_layer        &
        , nd_direction)                                                 &
!       Cloudy phase function for singly scattered solar radiation
    , forward_scatter_cloud(nd_profile, id_ct: nd_layer)                &
!       Cloudy forward scattering
    , forward_solar_cloud(nd_radiance_profile, id_ct: nd_layer)
!       Cloudy forward scattering for the solar beam

! Local variables.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , ll                                                                &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , ls
!       Loop variable
  REAL (RealK) ::                                                       &
      asymmetry_process(nd_profile, id_ct: nd_layer)                    &
!       Asymmetry factor for current process
    , omega                                                             &
!       Albedo of single scattering for the current process
    , x                                                                 &
!       Temporary algebraic variable
    , y                                                                 &
!       Temporary algebraic variable
    , t_celsius                                                         &
!       Temperature in celsius
    , temp_correction                                                   &
!       Temperature correction
    , phf_tmp                                                           &
!       Temporary phase function
    , smoothing
!       Lanczos (cosine) smoothing factor for truncated series

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
    , ks_phf(nd_radiance_profile)
!       Product of the scattering and the current moment of
!       the phase function

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'OPT_PROP_ICE_CLOUD'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


  IF ( (n_order_phase == 1) .AND.                                       &
    (i_parametrization_ice == ip_ice_pade_2_phf) .AND.                  &
    l_rescale .AND. (n_order_forward == 2) ) THEN

    CALL opt_prop_pade_2(id_ct, &
      n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
      ice_cloud_parameter, ice_mass_frac, dim_char_ice, &
      k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
      phase_fnc_cloud, forward_scatter_cloud)

  ELSE IF ( (n_order_phase == 1) .AND.                                  &
    (i_parametrization_ice == ip_ice_fu_phf) .AND.                      &
    l_rescale .AND. (n_order_forward == 2) ) THEN

    CALL opt_prop_fu_phf(id_ct, &
      n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
      ice_cloud_parameter, ice_mass_frac, dim_char_ice, &
      k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
      phase_fnc_cloud, forward_scatter_cloud)

  ELSE IF ( (n_order_phase == 1) .AND.                                  &
    (i_parametrization_ice == ip_ice_baran) .AND.                       &
    l_rescale .AND. (n_order_forward == 2) ) THEN

    CALL opt_prop_baran(id_ct, &
      n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
      ice_cloud_parameter, ice_mass_frac, dim_char_ice, t, &
      k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
      phase_fnc_cloud, forward_scatter_cloud)

  ELSE IF ( (i_parametrization_ice == ip_slingo_schrecker_ice).OR.      &
       (i_parametrization_ice == ip_ice_adt).OR.                        &
       (i_parametrization_ice == ip_ice_adt_10).OR.                     &
       (i_parametrization_ice == ip_ice_fu_solar).OR.                   &
       (i_parametrization_ice == ip_ice_fu_ir).OR.                      &
       (i_parametrization_ice == ip_sun_shine_vn2_ir).OR.               &
       (i_parametrization_ice == ip_sun_shine_vn2_ir).OR.               &
       (i_parametrization_ice == ip_ice_t_iwc).OR.                      &
       (i_parametrization_ice == ip_ice_iwc_only).OR.                   &
       (i_parametrization_ice == ip_ice_baran).OR.                      &
       ( l_henyey_greenstein_pf .AND.                                   &
         (i_parametrization_ice == ip_ice_pade_2_phf) ).OR.             &
       ( l_henyey_greenstein_pf .AND.                                   &
         (i_parametrization_ice == ip_slingo_schr_ice_phf) ).OR.        &
       ( l_henyey_greenstein_pf .AND.                                   &
         (i_parametrization_ice == ip_ice_fu_phf) ) ) THEN

!   Optical properties are calculated from parametrized data.
    SELECT CASE(i_parametrization_ice)

    CASE(ip_slingo_schrecker_ice, ip_slingo_schr_ice_phf)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)                                         &
            =ice_mass_frac(l, i)*(ice_cloud_parameter(1)                &
            +ice_cloud_parameter(2)/dim_char_ice(l, i))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0_RealK-ice_cloud_parameter(3)                          &
          -ice_cloud_parameter(4)*dim_char_ice(l, i))
          asymmetry_process(l, i)                                       &
            =ice_cloud_parameter(5)+ice_cloud_parameter(6)              &
            *dim_char_ice(l, i)
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE (ip_ice_adt)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          x=LOG(dim_char_ice(l, i)/ice_cloud_parameter(10))
          IF (x >  0.0_RealK) THEN
!           Large mode.
            k_ext_tot_cloud(l, i)=ice_mass_frac(l, i)                   &
              *EXP(ice_cloud_parameter(1)                               &
              +x*(ice_cloud_parameter(2)                                &
              +x*(ice_cloud_parameter(3)                                &
              +x*(ice_cloud_parameter(4)                                &
              +x*ice_cloud_parameter(5)))))
          ELSE IF (x <= 0.0_RealK) THEN
!           Small mode.
            k_ext_tot_cloud(l, i)=ice_mass_frac(l, i)                   &
              *EXP(ice_cloud_parameter(1)                               &
              +x*(ice_cloud_parameter(6)                                &
              +x*(ice_cloud_parameter(7)                                &
              +x*(ice_cloud_parameter(8)                                &
              +x*ice_cloud_parameter(9)))))
          END IF
          x=LOG(dim_char_ice(l, i)/ice_cloud_parameter(20))
          IF (x >  0.0_RealK) THEN
!           Large mode.
            k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                &
              *(1.0_RealK-(ice_cloud_parameter(11)                      &
              +x*(ice_cloud_parameter(12)                               &
              +x*(ice_cloud_parameter(13)                               &
              +x*(ice_cloud_parameter(14)                               &
              +x*ice_cloud_parameter(15))))))
          ELSE IF (x <= 0.0_RealK) THEN
!           Small mode.
            k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                &
              *(1.0_RealK-(ice_cloud_parameter(11)                      &
              +x*(ice_cloud_parameter(16)                               &
              +x*(ice_cloud_parameter(17)                               &
              +x*(ice_cloud_parameter(18)                               &
              +x*ice_cloud_parameter(19))))))
          END IF
          x=LOG(dim_char_ice(l, i)/ice_cloud_parameter(30))
          IF (x >  0.0_RealK) THEN
!           Large mode.
            asymmetry_process(l, i)=ice_cloud_parameter(21)             &
              +x*(ice_cloud_parameter(22)                               &
              +x*(ice_cloud_parameter(23)                               &
              +x*(ice_cloud_parameter(24)                               &
              +x*ice_cloud_parameter(25))))
          ELSE IF (x <= 0.0_RealK) THEN
!           Small mode.
            asymmetry_process(l, i)=ice_cloud_parameter(21)             &
              +x*(ice_cloud_parameter(26)                               &
              +x*(ice_cloud_parameter(27)                               &
              +x*(ice_cloud_parameter(28)                               &
              +x*ice_cloud_parameter(29))))
          END IF
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_ice_adt_10)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          x=dim_char_ice(l, i)/ice_cloud_parameter(12)
          y=ice_cloud_parameter(6)                                      &
            +x*(ice_cloud_parameter(7)                                  &
            +x*(ice_cloud_parameter(8)                                  &
            +x*(ice_cloud_parameter(9)                                  &
            +x*(ice_cloud_parameter(10)                                 &
            +x*ice_cloud_parameter(11)))))
          k_ext_tot_cloud(l, i)=ice_mass_frac(l, i)                     &
            *EXP(ice_cloud_parameter(1)                                 &
            +x*(ice_cloud_parameter(2)                                  &
            +x*(ice_cloud_parameter(3)                                  &
            +x*(ice_cloud_parameter(4)                                  &
            +x*(ice_cloud_parameter(5)                                  &
            +x*y)))))
          x=dim_char_ice(l, i)/ice_cloud_parameter(24)
          y=ice_cloud_parameter(18)                                     &
            +x*(ice_cloud_parameter(19)                                 &
            +x*(ice_cloud_parameter(20)                                 &
            +x*(ice_cloud_parameter(21)                                 &
            +x*(ice_cloud_parameter(22)                                 &
            +x*ice_cloud_parameter(23)))))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0_RealK-(ice_cloud_parameter(13)                        &
            +x*(ice_cloud_parameter(14)                                 &
            +x*(ice_cloud_parameter(15)                                 &
            +x*(ice_cloud_parameter(16)                                 &
            +x*(ice_cloud_parameter(17)                                 &
            +x*y))))))
          x=dim_char_ice(l, i)/ice_cloud_parameter(36)
          y=ice_cloud_parameter(30)                                     &
            +x*(ice_cloud_parameter(31)                                 &
            +x*(ice_cloud_parameter(32)                                 &
            +x*(ice_cloud_parameter(33)                                 &
            +x*(ice_cloud_parameter(34)                                 &
            +x*ice_cloud_parameter(35)))))
          asymmetry_process(l, i)=ice_cloud_parameter(25)               &
            +x*(ice_cloud_parameter(26)                                 &
            +x*(ice_cloud_parameter(27)                                 &
            +x*(ice_cloud_parameter(28)                                 &
            +x*(ice_cloud_parameter(29)                                 &
            +x*y))))
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE (ip_sun_shine_vn2_vis)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          t_celsius=t(l, i)-2.7316e+02_RealK
          temp_correction=1.047_RealK                                   &
            +t_celsius*(-9.13e-05_RealK+t_celsius                       &
            *(2.026e-04_RealK-1.056e-05_RealK*t_celsius))
          k_ext_tot_cloud(l, i)=temp_correction*ice_mass_frac(l, i)     &
            /(3.05548e-02_RealK                                         &
            +2.54802e+02_RealK*density(l, i)*ice_mass_frac(l, i))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0_RealK-ice_cloud_parameter(1)                          &
            *EXP(ice_cloud_parameter(2)                                 &
            *LOG(density(l, i)*ice_mass_frac(l, i)+1.0e-12_RealK)))     &
            *(1.0_RealK+ice_cloud_parameter(5)                          &
            *(temp_correction-1.0_RealK)/temp_correction)
          asymmetry_process(l, i)                                       &
            =ice_cloud_parameter(3)*EXP(ice_cloud_parameter(4)          &
            *LOG(density(l, i)*ice_mass_frac(l, i)+1.0e-12_RealK))      &
            *(1.0_RealK+ice_cloud_parameter(6)                          &
            *(temp_correction-1.0_RealK)/temp_correction)
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_sun_shine_vn2_ir)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          t_celsius=t(l, i)-2.7316e+02_RealK
          temp_correction=1.047_RealK+t_celsius                         &
            *(-9.13e-05_RealK+t_celsius                                 &
            *(2.026e-04_RealK-1.056e-05_RealK*t_celsius))
          k_ext_tot_cloud(l, i)=temp_correction*ice_mass_frac(l, i)     &
            /(6.30689e-02_RealK                                         &
            +2.65874e+02_RealK*density(l, i)*ice_mass_frac(l, i))
        END DO
      END DO

    CASE(ip_ice_fu_solar)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)                                         &
            =ice_mass_frac(l, i)*(ice_cloud_parameter(1)                &
            +ice_cloud_parameter(2)/dim_char_ice(l, i))
          omega=1.0_RealK-(ice_cloud_parameter(3)                       &
            +dim_char_ice(l, i)*(ice_cloud_parameter(4)                 &
            +dim_char_ice(l, i)*(ice_cloud_parameter(5)                 &
            +dim_char_ice(l, i)*ice_cloud_parameter(6))))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)*omega
          asymmetry_process(l, i)                                       &
            =ice_cloud_parameter(7)+dim_char_ice(l, i)                  &
            *(ice_cloud_parameter(8)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(9)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(10))))
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
!         The forward scattering will be limited later.
          forward_scatter_cloud(l, i)=k_ext_scat_cloud(l, i)            &
            *(1.0_RealK                                                 &
            /MAX(1.0_RealK, 2.0_RealK*omega)                            &
            +ice_cloud_parameter(11)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(12)+dim_char_ice(l, i)                &
            *(ice_cloud_parameter(13)+dim_char_ice(l, i)                &
            *(ice_cloud_parameter(14)))))
        END DO
      END DO

    CASE(ip_ice_fu_ir)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)=ice_mass_frac(l, i)                     &
            *((ice_cloud_parameter(3)/dim_char_ice(l, i)                &
            +ice_cloud_parameter(2))/dim_char_ice(l, i)                 &
            +ice_cloud_parameter(1))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            -(ice_mass_frac(l, i)/dim_char_ice(l, i))                   &
            *(ice_cloud_parameter(4)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(5)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(6)+dim_char_ice(l, i)                 &
            *ice_cloud_parameter(7))))
          asymmetry_process(l, i)                                       &
            =ice_cloud_parameter(8)+dim_char_ice(l, i)                  &
            *(ice_cloud_parameter(9)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(10)+dim_char_ice(l, i)                &
            *(ice_cloud_parameter(11))))
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_ice_fu_phf)
      CALL opt_prop_fu_phf(id_ct, &
        n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
        ice_cloud_parameter, ice_mass_frac, dim_char_ice, &
        k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
        phase_fnc_cloud)

    CASE(ip_ice_t_iwc)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          x=LOG(MIN(ice_mass_frac(l, i)+1.0e-12_RealK, 1.0e-01_RealK))
          k_ext_tot_cloud(l, i)=EXP(                                    &
             ice_cloud_parameter(1)                                     &
            +ice_cloud_parameter(2)*t(l, i)                             &
            +ice_cloud_parameter(3)*x )
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)*MIN(             &
             ice_cloud_parameter(4)                                     &
            +ice_cloud_parameter(5)*t(l, i)                             &
            +ice_cloud_parameter(6)*x , 1.0_RealK)
          asymmetry_process(l, i)=MIN(                                  &
             ice_cloud_parameter(7)                                     &
            +ice_cloud_parameter(8)*t(l, i)                             &
            +ice_cloud_parameter(9)*x , 1.0_RealK-EPSILON(1.0_RealK))
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_ice_iwc_only)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          x=MIN(ice_mass_frac(l, i)+1.0e-12_RealK, 1.0e-01_RealK)
          k_ext_tot_cloud(l, i)=                                        &
             ice_cloud_parameter(1)*x
          k_ext_scat_cloud(l, i)=                                       &
             ice_cloud_parameter(2)*x
          asymmetry_process(l, i)=MIN(                                  &
             ice_cloud_parameter(3)*EXP(ice_cloud_parameter(4)          &
             *LOG(x)), 1.0_RealK-EPSILON(1.0_RealK))
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_ice_baran)
      CALL opt_prop_baran(id_ct, &
        n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
        ice_cloud_parameter, ice_mass_frac, dim_char_ice, t, &
        k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
        phase_fnc_cloud)

    CASE(ip_ice_pade_2_phf)
      CALL opt_prop_pade_2(id_ct, &
        n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
        ice_cloud_parameter, ice_mass_frac, dim_char_ice, &
        k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
        phase_fnc_cloud)

    END SELECT


!   Parametrizations which do not include explicit
!   representation of the higher moments are extended using the
!   Henyey-Greenstein phase function.
    DO ls=2, n_order_phase
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          phase_fnc_cloud(l, i, ls)                                     &
            =phase_fnc_cloud(l, i, ls-1)*asymmetry_process(l, i)
        END DO
      END DO
    END DO

    IF (l_rescale) THEN
!     For most parameterizations the forward scattering
!     is determined from the asymmetry, but in the case of
!     Fu's parametrization it is defined specially, but must
!     be limited to avoid negative moments in the phase function.
      IF (i_parametrization_ice == ip_ice_fu_solar) THEN
        DO i=n_cloud_top, n_layer
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            forward_scatter_cloud(l, i)                                 &
              =MIN(forward_scatter_cloud(l, i)                          &
              , k_ext_scat_cloud(l, i)                                  &
              *asymmetry_process(l, i)**(n_order_forward-1))
          END DO
        END DO
      ELSE
        DO i=n_cloud_top, n_layer
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            forward_scatter_cloud(l, i)                                 &
              =k_ext_scat_cloud(l, i)                                   &
              *asymmetry_process(l, i)**n_order_forward
          END DO
        END DO
      END IF
    END IF

    IF (l_solar_phf) THEN
!     Calculate the solar phase function to higher accuracy.
      DO i=n_cloud_top, n_layer
        DO id=1, n_direction
!         The Legendre polynomials are not stored so as to reduce
!         the requirement for memory at very high orders of solar
!         truncation.
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
!           Initialize the Legendre polynomials at the zeroth and
!           first orders.
            p_legendre_ls_m1(l)=1.0_RealK
            p_legendre_ls(l)=cos_sol_view(l, id)
            ks_phf(l)=k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
            phase_fnc_solar_cloud(l, i, id)=k_ext_scat_cloud(l, i)      &
              +ks_phf(l)*p_legendre_ls(l)*REAL(2*1+1, RealK)
          END DO

          DO ls=2, n_order_phase_solar
!           Calculate higher orders by recurrences.
            cnst1=1.0_RealK-1.0_RealK/REAL(ls, RealK)
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              p_legendre_tmp(l)=p_legendre_ls(l)
              p_legendre_ls(l)                                          &
                =(1.0_RealK+cnst1)*p_legendre_ls(l)                     &
                *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
              p_legendre_ls_m1(l)=p_legendre_tmp(l)
              ks_phf(l)=ks_phf(l)*asymmetry_process(l, i)
              phase_fnc_solar_cloud(l, i, id)                           &
                =phase_fnc_solar_cloud(l, i, id)                        &
                +ks_phf(l)*p_legendre_ls(l)                             &
                *REAL(2*ls+1, RealK)
            END DO

          END DO
        END DO

!       Continue to an extra order to find the rescaling
!       for the solar beam.
        IF (l_rescale) THEN
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            forward_solar_cloud(l, i)                                   &
              =ks_phf(l)*asymmetry_process(l, i)
          END DO
        END IF
      END DO
    END IF

  ELSE IF ( .NOT.l_henyey_greenstein_pf .AND.                           &
       ( (i_parametrization_ice == ip_slingo_schr_ice_phf).OR.          &
         (i_parametrization_ice == ip_ice_fu_phf).OR.                   &
         (i_parametrization_ice == ip_ice_pade_2_phf) ) ) THEN

    DO i=n_cloud_top, n_layer


!     To avoid the repetition of blocks of code or excessive
!     use of memory it is easiest to have an outer loop over
!     layers.


      IF (i_parametrization_ice == ip_slingo_schr_ice_phf) THEN

        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)                                         &
            =ice_mass_frac(l, i)*(ice_cloud_parameter(1)                &
            +ice_cloud_parameter(2)/dim_char_ice(l, i))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0_RealK-ice_cloud_parameter(3)                          &
          -ice_cloud_parameter(4)*dim_char_ice(l, i))
        END DO
        DO ls=1, n_order_phase
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            phase_fnc_cloud(l, i, ls)                                   &
              =k_ext_scat_cloud(l, i)*(ice_cloud_parameter(2*ls+3)      &
              +ice_cloud_parameter(2*ls+4)*dim_char_ice(l, i))
          END DO
        END DO
        ls=n_order_forward
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          forward_scatter_cloud(l, i)                                   &
            =k_ext_scat_cloud(l, i)*(ice_cloud_parameter(2*ls+3)        &
            +ice_cloud_parameter(2*ls+4)*dim_char_ice(l, i))
        END DO

      ELSE IF (i_parametrization_ice == ip_ice_fu_phf) THEN

        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          x=dim_char_ice(l, i)/ice_cloud_parameter(4)
          k_ext_tot_cloud(l, i)=ice_mass_frac(l, i)                     &
            *((ice_cloud_parameter(3)/x                                 &
            +ice_cloud_parameter(2))/x                                  &
            +ice_cloud_parameter(1))
          x=dim_char_ice(l, i)/ice_cloud_parameter(9)
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0_RealK                                                 &
            -(ice_cloud_parameter(5)+x                                  &
            *(ice_cloud_parameter(6)+x                                  &
            *(ice_cloud_parameter(7)+x                                  &
            *ice_cloud_parameter(8)))))
        END DO
        DO ls=1, n_order_phase
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            x=dim_char_ice(l, i)/ice_cloud_parameter(5*ls+9)
            phase_fnc_cloud(l, i, ls)                                   &
              =k_ext_scat_cloud(l, i)*(ice_cloud_parameter(5*ls+5)      &
              +x*(ice_cloud_parameter(5*ls+6)                           &
              +x*(ice_cloud_parameter(5*ls+7)                           &
              +x*ice_cloud_parameter(5*ls+8))))
          END DO
        END DO
        ls=n_order_forward
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          x=dim_char_ice(l, i)/ice_cloud_parameter(5*ls+9)
          forward_scatter_cloud(l, i)                                   &
            =k_ext_scat_cloud(l, i)*(ice_cloud_parameter(5*ls+5)        &
            +x*(ice_cloud_parameter(5*ls+6)                             &
            +x*(ice_cloud_parameter(5*ls+7)                             &
            +x*ice_cloud_parameter(5*ls+8))))
        END DO

      ELSE IF (i_parametrization_ice == ip_ice_pade_2_phf) THEN

        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)=ice_mass_frac(l, i)                     &
            *(ice_cloud_parameter(1)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(2)+dim_char_ice(l, i)                 &
            *ice_cloud_parameter(3)))                                   &
            /(1.0e+00_RealK+dim_char_ice(l, i)                          &
            *(ice_cloud_parameter(4)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(5)+dim_char_ice(l, i)                 &
            *ice_cloud_parameter(6))))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0e+00_RealK                                             &
            -(ice_cloud_parameter(7)+dim_char_ice(l, i)                 &
            *(ice_cloud_parameter(8)+dim_char_ice(l, i)                 &
            *ice_cloud_parameter(9)))                                   &
            /(1.0e+00_RealK+dim_char_ice(l, i)                          &
            *(ice_cloud_parameter(10)+dim_char_ice(l, i)                &
            *ice_cloud_parameter(11))))
        END DO
        DO ls=1, n_order_phase
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            phase_fnc_cloud(l, i, ls)=k_ext_scat_cloud(l, i)            &
              *(ice_cloud_parameter(5*ls+7)+dim_char_ice(l, i)          &
              *(ice_cloud_parameter(5*ls+8)+dim_char_ice(l, i)          &
              *ice_cloud_parameter(5*ls+9)))                            &
              /(1.0e+00_RealK+dim_char_ice(l, i)                        &
              *(ice_cloud_parameter(5*ls+10)+dim_char_ice(l, i)         &
              *ice_cloud_parameter(5*ls+11)))
          END DO
        END DO
        ls=n_order_forward
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          forward_scatter_cloud(l, i)=k_ext_scat_cloud(l, i)            &
            *(ice_cloud_parameter(5*ls+7)+dim_char_ice(l, i)            &
            *(ice_cloud_parameter(5*ls+8)+dim_char_ice(l, i)            &
            *ice_cloud_parameter(5*ls+9)))                              &
            /(1.0e+00_RealK+dim_char_ice(l, i)                          &
            *(ice_cloud_parameter(5*ls+10)+dim_char_ice(l, i)           &
            *ice_cloud_parameter(5*ls+11)))
        END DO

      END IF


      IF (l_solar_phf) THEN
!       Calculate the solar phase function to higher accuracy.
        DO id=1, n_direction
!         The Legendre polynomials are not stored so as to reduce
!         the requirement for memory at very high orders of solar
!         truncation.
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
!           Initialize the Legendre polynomials at the zeroth and
!           first orders.
            p_legendre_ls_m1(l)=1.0_RealK
            p_legendre_ls(l)=cos_sol_view(l, id)
            phase_fnc_solar_cloud(l, i, id)=k_ext_scat_cloud(l, i)      &
              +phase_fnc_cloud(l, i, 1)                                 &
              *p_legendre_ls(l)*REAL(2*1+1, RealK)
          END DO

!         Calculate higher orders by recurrences. Moments of
!         the phase function cannot be taken from above as
!         we will typically require a much higher order here.
          DO ls=2, n_order_phase_solar

            IF (l_lanczos) THEN
!             Cosine filter
              smoothing = COS( REAL(ls, RealK) * pi /                   &
                         (2.0 * n_order_phase_solar) )
            ELSE
              smoothing = 1.0e+00_RealK
            END IF

            cnst1=1.0_RealK-1.0_RealK/REAL(ls, RealK)
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              p_legendre_tmp(l)=p_legendre_ls(l)
              p_legendre_ls(l)                                          &
                =(1.0_RealK+cnst1)*p_legendre_ls(l)                     &
                *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
              p_legendre_ls_m1(l)=p_legendre_tmp(l)

              SELECT CASE(i_parametrization_ice)

                CASE(ip_slingo_schr_ice_phf)
                  phf_tmp=ice_cloud_parameter(2*ls+3)                   &
                    +dim_char_ice(l, i)                                 &
                    *ice_cloud_parameter(2*ls+4)

                CASE(ip_ice_fu_phf)
                  x=dim_char_ice(l, i)/ice_cloud_parameter(5*ls+9)
                  phf_tmp                                               &
                    =(ice_cloud_parameter(5*ls+5)                       &
                    +x*(ice_cloud_parameter(5*ls+6)                     &
                    +x*(ice_cloud_parameter(5*ls+7)                     &
                    +x*ice_cloud_parameter(5*ls+8))))

                CASE(ip_ice_pade_2_phf)
                  phf_tmp                                               &
                    =(ice_cloud_parameter(5*ls+7)+dim_char_ice(l, i)    &
                    *(ice_cloud_parameter(5*ls+8)+dim_char_ice(l, i)    &
                    *ice_cloud_parameter(5*ls+9)))                      &
                    /(1.0e+00_RealK+dim_char_ice(l, i)                  &
                    *(ice_cloud_parameter(5*ls+10)+dim_char_ice(l, i)   &
                    *ice_cloud_parameter(5*ls+11)))

              END SELECT

              ks_phf(l)=k_ext_scat_cloud(l, i)*phf_tmp
              phase_fnc_solar_cloud(l, i, id)                           &
                =phase_fnc_solar_cloud(l, i, id)                        &
                +ks_phf(l)*p_legendre_ls(l)                             &
                *REAL(2*ls+1, RealK)                                    &
                *smoothing
            END DO


          END DO
        END DO

!       Continue to an extra order to find the rescaling
!       for the solar beam.
        IF (l_rescale) THEN
          ls=n_order_phase_solar+1

          IF (l_lanczos) THEN
!           Cosine filter
            smoothing = COS( REAL(n_order_phase_solar, RealK)           &
                        * pi / (2.0 * n_order_phase_solar) )
          ELSE
            smoothing = 1.0e+00_RealK
          END IF

          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)

            SELECT CASE(i_parametrization_ice)

              CASE(ip_slingo_schr_ice_phf)
                phf_tmp=ice_cloud_parameter(2*ls+3)                     &
                  +dim_char_ice(l, i)                                   &
                  *ice_cloud_parameter(2*ls+4)

              CASE(ip_ice_fu_phf)
                x=dim_char_ice(l, i)/ice_cloud_parameter(5*ls+9)
                phf_tmp                                                 &
                  =(ice_cloud_parameter(5*ls+5)                         &
                  +x*(ice_cloud_parameter(5*ls+6)                       &
                  +x*(ice_cloud_parameter(5*ls+7)                       &
                  +x*ice_cloud_parameter(5*ls+8))))

              CASE(ip_ice_pade_2_phf)
                phf_tmp                                                 &
                  =(ice_cloud_parameter(5*ls+7)+dim_char_ice(l, i)      &
                  *(ice_cloud_parameter(5*ls+8)+dim_char_ice(l, i)      &
                  *ice_cloud_parameter(5*ls+9)))                        &
                  /(1.0e+00_RealK+dim_char_ice(l, i)                    &
                  *(ice_cloud_parameter(5*ls+10)+dim_char_ice(l, i)     &
                  *ice_cloud_parameter(5*ls+11)))

            END SELECT

            forward_solar_cloud(l, i)                                   &
              =k_ext_scat_cloud(l, i)*phf_tmp                           &
              *smoothing
          END DO
        END IF

      END IF

    END DO

  ELSE IF (i_parametrization_ice == ip_ice_unparametrized) THEN
    CALL prsc_opt_prop(ierr                                             &
      , n_profile, n_cloud_top, n_layer                                 &
      , l_rescale, n_order_forward                                      &
      , l_henyey_greenstein_pf, n_order_phase                           &
      , p, density                                                      &
      , n_opt_level_cloud_prsc                                          &
      , ice_pressure_prsc, ice_absorption_prsc, ice_scattering_prsc     &
      , ice_phase_fnc_prsc                                              &
      , k_ext_tot_cloud, k_ext_scat_cloud, phase_fnc_cloud              &
      , forward_scatter_cloud, forward_solar_cloud                      &
      , l_solar_phf, n_order_phase_solar, n_direction, cos_sol_view     &
      , phase_fnc_solar_cloud                                           &
      , nd_profile, nd_radiance_profile, nd_layer, id_ct, nd_layer      &
      , nd_direction                                                    &
      , nd_profile_prsc, nd_opt_level_prsc                              &
      , nd_phase_term, nd_max_order                                     &
      )

!   The absorption is returned from prsc_opt_prop in k_ext_tot_cloud.
!   The scattering is added here to give the correct extinction.
    DO i=n_cloud_top, n_layer
      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i) = k_ext_tot_cloud(l, i)                 &
            + k_ext_scat_cloud(l, i)
      END DO
    END DO

  ELSE

    cmessage = '*** Error: An invalid parametrization '                 &
      //'of ice crystals has been used.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE opt_prop_ice_cloud
END MODULE opt_prop_ice_cloud_mod

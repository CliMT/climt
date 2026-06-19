! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate optical properties of water clouds.
!
! Method:
!   If the optical properties come from an observational
!   distribution a separate subroutine is called. Otherwise
!   appropriate mean quantities in the layer are calculated
!   as the parametrization requires and these values are
!   substituted into the parametrization to give the optical
!   properties.
!
!   Note that this routine produces optical propeties for a
!   single condensed component of the cloud.
!
!- ---------------------------------------------------------------------
MODULE opt_prop_water_cloud_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OPT_PROP_WATER_CLOUD_MOD'
CONTAINS
SUBROUTINE opt_prop_water_cloud(ierr                                    &
    , n_profile, n_layer, n_cloud_top                                   &
    , n_cloud_profile, i_cloud_profile                                  &
    , n_order_phase, l_rescale, n_order_forward                         &
    , l_henyey_greenstein_pf, l_solar_phf, l_lanczos                    &
    , n_order_phase_solar, n_direction, cos_sol_view                    &
    , i_parametrization_drop, cloud_parameter                           &
    , liq_water_mass_frac, radius_effect                                &
    , p, density                                                        &
    , n_opt_level_cloud_prsc                                            &
    , drop_pressure_prsc, drop_absorption_prsc                          &
    , drop_scattering_prsc, drop_phase_fnc_prsc                         &
    , k_ext_tot_cloud, k_ext_scat_cloud                                 &
    , phase_fnc_cloud, forward_scatter_cloud                            &
    , forward_solar_cloud, phase_fnc_solar_cloud                        &
    , nd_profile, nd_radiance_profile, nd_layer, id_ct                  &
    , nd_direction                                                      &
    , nd_phase_term, nd_max_order, nd_cloud_parameter                   &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_slingo_schrecker,                  &
                     ip_ackerman_stephens, ip_drop_pade_2,              &
                     ip_drop_unparametrized, ip_ps_size_phf,            &
                     ip_slingo_schr_phf
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE opt_prop_pade_2_mod, ONLY: opt_prop_pade_2
  USE prsc_opt_prop_mod, ONLY: prsc_opt_prop

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_radiance_profile                                               &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_phase_term                                                     &
!       Size allocated for terms in phase function
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_cloud_parameter                                                &
!       Size allocated for cloud parameters
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
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_order_phase                                                     &
!       Number of terms to retain in the phase function
    , n_order_phase_solar                                               &
!       Number of terms to retain in single scattered solar
!       phase function
    , n_order_forward                                                   &
!       Order used in forming the forward scattering parameter
    , i_parametrization_drop                                            &
!       Treatment of droplets
    , n_opt_level_cloud_prsc                                            &
!       Number of levels of prescribed optical properties
    , n_cloud_profile(id_ct: nd_layer)                                  &
!       Number of cloudy profiles
    , i_cloud_profile(nd_profile, id_ct: nd_layer)
!       Profiles containing clouds
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale                                                         &
!       Flag for delta-rescaling
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
      cloud_parameter(nd_cloud_parameter)                               &
!       Cloud parameters
    , liq_water_mass_frac(nd_profile, id_ct: nd_layer)                  &
!       Liquid water content
    , radius_effect(nd_profile, id_ct: nd_layer)                        &
!       Effective radius
    , p(nd_profile, nd_layer)                                           &
!       Pressure
    , density(nd_profile, nd_layer)                                     &
!       Atmospheric density
    , drop_pressure_prsc(nd_profile_prsc, nd_opt_level_prsc)            &
!       Pressure levels where optical data are prescribed
    , drop_absorption_prsc(nd_profile_prsc, nd_opt_level_prsc)          &
!       Prescribed absorption by droplets
    , drop_scattering_prsc(nd_profile_prsc, nd_opt_level_prsc)          &
!       Prescribed scattering by droplets
    , drop_phase_fnc_prsc(nd_profile_prsc, nd_opt_level_prsc            &
        , nd_phase_term)
!       Prescribed phase function of droplets

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
!       Asymmetry of current process.
    , phf_tmp                                                           &
!       Temporary Phase Function
    , sz(nd_radiance_profile)                                           &
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
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'OPT_PROP_WATER_CLOUD'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF ( (n_order_phase == 1) .AND.                                       &
    (i_parametrization_drop == ip_drop_pade_2) .AND.                    &
    l_rescale .AND. (n_order_forward == 2) ) THEN

    CALL opt_prop_pade_2(id_ct, &
      n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
      cloud_parameter, liq_water_mass_frac, radius_effect, &
      k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
      phase_fnc_cloud, forward_scatter_cloud)

  ELSE IF ( (i_parametrization_drop == ip_slingo_schrecker).OR.         &
       (i_parametrization_drop == ip_ackerman_stephens).OR.             &
       (i_parametrization_drop == ip_drop_pade_2) .OR.                  &
       ( l_henyey_greenstein_pf .AND.                                   &
         (i_parametrization_drop == ip_slingo_schr_phf) ) ) THEN

!   Optical properties are calculated from parametrized data.
    SELECT CASE(i_parametrization_drop)

    CASE(ip_slingo_schrecker, ip_slingo_schr_phf)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)                                         &
            =liq_water_mass_frac(l, i)*(cloud_parameter(1)              &
            +cloud_parameter(2)/radius_effect(l, i))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0e+00_RealK-cloud_parameter(3)                          &
            -cloud_parameter(4)*radius_effect(l, i))
          asymmetry_process(l, i)=                                      &
            cloud_parameter(5)+cloud_parameter(6)                       &
            *radius_effect(l, i)
          phase_fnc_cloud(l, i, 1)=                                     &
            k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_ackerman_stephens)
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          k_ext_tot_cloud(l, i)=liq_water_mass_frac(l, i)               &
            *(cloud_parameter(1)+cloud_parameter(2)                     &
            *EXP(cloud_parameter(3)*LOG(radius_effect(l, i))))
          k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                  &
            *(1.0e+00_RealK-cloud_parameter(4)                          &
            -cloud_parameter(5)*EXP(cloud_parameter(6)                  &
            *LOG(radius_effect(l, i))))
          asymmetry_process(l, i)                                       &
            =cloud_parameter(7)+cloud_parameter(8)                      &
            *EXP(cloud_parameter(9)*LOG(radius_effect(l, i)))
          phase_fnc_cloud(l, i, 1)                                      &
            =k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
        END DO
      END DO

    CASE(ip_drop_pade_2)
      CALL opt_prop_pade_2(id_ct, &
        n_cloud_top, n_layer, n_cloud_profile, i_cloud_profile, &
        cloud_parameter, liq_water_mass_frac, radius_effect, &
        k_ext_tot_cloud, k_ext_scat_cloud, asymmetry_process, &
        phase_fnc_cloud)

    END SELECT


!   Since these parametrizations include only the asymmetry,
!   it seems reasonable to extend them to higher
!   truncations using the Henyey-Greenstein phase function.
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
      DO i=n_cloud_top, n_layer
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
            forward_scatter_cloud(l, i)                                 &
              =k_ext_scat_cloud(l, i)                                   &
              *asymmetry_process(l, i)**n_order_forward
        END DO
      END DO
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
            p_legendre_ls_m1(l)=1.0e+00_RealK
            p_legendre_ls(l)=cos_sol_view(l, id)
            ks_phf(l)=k_ext_scat_cloud(l, i)*asymmetry_process(l, i)
            phase_fnc_solar_cloud(l, i, id)=k_ext_scat_cloud(l, i)      &
              +ks_phf(l)*p_legendre_ls(l)*REAL(2*1+1, RealK)
          END DO

          DO ls=2, n_order_phase_solar
!           Calculate higher orders by recurrences.
            cnst1=1.0e+00_RealK-1.0e+00_RealK/REAL(ls, RealK)
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              p_legendre_tmp(l)=p_legendre_ls(l)
              p_legendre_ls(l)                                          &
                =(1.0e+00_RealK+cnst1)*p_legendre_ls(l)                 &
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


  ELSE IF (.NOT. l_henyey_greenstein_pf .AND.                           &
          (i_parametrization_drop == ip_slingo_schr_phf) ) THEN

    DO i=n_cloud_top, n_layer

!     To avoid the repetition of blocks of code or excessive
!     use of memory it is easiest to have an outer loop over
!     layers

      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        k_ext_tot_cloud(l, i)                                           &
          =liq_water_mass_frac(l, i)*(cloud_parameter(1)                &
          +cloud_parameter(2)/radius_effect(l, i))
        k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                    &
          *(1.0e+00_RealK-cloud_parameter(3)                            &
          -cloud_parameter(4)*radius_effect(l, i))
      END DO

      DO ls=1, n_order_phase
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          phase_fnc_cloud(l, i, ls)                                     &
            =k_ext_scat_cloud(l,i)*(cloud_parameter(2*ls+3)             &
            +cloud_parameter(2*ls+4)*radius_effect(l,i))
        END DO
      END DO

      ls=n_order_forward

      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        forward_scatter_cloud(l,i)                                      &
          =k_ext_scat_cloud(l,i)*(cloud_parameter(2*ls+3)               &
          +cloud_parameter(2*ls+4)*radius_effect(l,i))
      END DO

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
            p_legendre_ls_m1(l)=1.0e+00_RealK
            p_legendre_ls(l)=cos_sol_view(l, id)
            phase_fnc_solar_cloud(l, i, id)=k_ext_scat_cloud(l, i)      &
               + phase_fnc_cloud(l, i, 1)                               &
               * p_legendre_ls(l)*REAL(2*1+1, RealK)
          END DO

!         Calculate higher orders by recurrences.
          DO ls=2, n_order_phase_solar

            IF (l_lanczos) THEN
!             Cosine filter
              smoothing = COS( REAL(ls, RealK) * pi /                   &
                         (2.0 * n_order_phase_solar) )
            ELSE
              smoothing = 1.0e+00_RealK
            END IF

            cnst1=1.0e+00_RealK-1.0e+00_RealK/REAL(ls, RealK)
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              p_legendre_tmp(l)=p_legendre_ls(l)
              p_legendre_ls(l)                                          &
                =(1.0e+00_RealK+cnst1)*p_legendre_ls(l)                 &
                *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
              p_legendre_ls_m1(l)=p_legendre_tmp(l)
              phf_tmp=cloud_parameter(2*ls+3)                           &
                     + radius_effect(l,i)                               &
                     * cloud_parameter(2*ls+4)
              IF (ls == n_order_phase_solar) phf_tmp=0.5*phf_tmp
              ks_phf(l)=k_ext_scat_cloud(l,i)*phf_tmp

              phase_fnc_solar_cloud(l, i, id)                           &
                = phase_fnc_solar_cloud(l, i, id)                       &
                + ks_phf(l)*p_legendre_ls(l)                            &
                * REAL(2*ls+1, RealK)                                   &
                * smoothing
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
            phf_tmp=cloud_parameter(2*ls+3)                             &
                  +radius_effect(l, i)                                  &
                  *cloud_parameter(2*ls+4)

!           Cosine weighting.
            forward_solar_cloud(l, i)                                   &
              = k_ext_scat_cloud(l, i)*0.5*phf_tmp                      &
              * smoothing
          END DO
        END IF

      END IF

    END DO

  ELSE IF (.NOT. l_henyey_greenstein_pf .AND.                           &
          (i_parametrization_drop == ip_ps_size_phf) ) THEN

    DO i=n_cloud_top, n_layer

!     To avoid the repetition of blocks of code or excessive
!     use of memory it is easiest to have an outer loop over
!     layers

      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        sz(l)=radius_effect(l, i)/cloud_parameter(1)
        sz(l)=MAX(0.2_RealK,MIN(sz(l),1.0_RealK))
        k_ext_tot_cloud(l, i)                                           &
          =liq_water_mass_frac(l, i)*(                                  &
          cloud_parameter(3)/sz(l)**3 +                                 &
          cloud_parameter(4)/sz(l)**2 +                                 &
          cloud_parameter(5)/sz(l)**1 +                                 &
          cloud_parameter(6) )
        k_ext_scat_cloud(l, i)=k_ext_tot_cloud(l, i)                    &
          *(1.0e+00_RealK-cloud_parameter(7)                            &
          -cloud_parameter(8)*sz(l)                                     &
          -cloud_parameter(9)*sz(l)**2                                  &
          -cloud_parameter(10)*sz(l)**3                                 &
          -cloud_parameter(11)*sz(l)**4)
      END DO

      DO ls=1, n_order_phase
        DO ll=1, n_cloud_profile(i)
          l=i_cloud_profile(ll, i)
          phase_fnc_cloud(l, i, ls)                                     &
            =k_ext_scat_cloud(l,i)*(                                    &
            cloud_parameter(5*ls+7)/sz(l)**4 +                          &
            cloud_parameter(5*ls+8)/sz(l)**3 +                          &
            cloud_parameter(5*ls+9)/sz(l)**2 +                          &
            cloud_parameter(5*ls+10)/sz(l) +                            &
            cloud_parameter(5*ls+11) ) *                                &
            EXP(-cloud_parameter(2)*REAL(ls)/sz(l))
        END DO
      END DO

      ls=n_order_forward

      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        forward_scatter_cloud(l,i)                                      &
          =k_ext_scat_cloud(l,i)*(                                      &
          cloud_parameter(5*ls+7)/sz(l)**4 +                            &
          cloud_parameter(5*ls+8)/sz(l)**3 +                            &
          cloud_parameter(5*ls+9)/sz(l)**2 +                            &
          cloud_parameter(5*ls+10)/sz(l) +                              &
          cloud_parameter(5*ls+11) ) *                                  &
          EXP(-cloud_parameter(2)*REAL(ls)/sz(l))
      END DO

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
            p_legendre_ls_m1(l)=1.0e+00_RealK
            p_legendre_ls(l)=cos_sol_view(l, id)
            phase_fnc_solar_cloud(l, i, id)=k_ext_scat_cloud(l, i)      &
               + phase_fnc_cloud(l, i, 1)                               &
               * p_legendre_ls(l)*REAL(2*1+1, RealK)
          END DO

!         Calculate higher orders by recurrences.
          DO ls=2, n_order_phase_solar

            IF (l_lanczos) THEN
!             Cosine filter
              smoothing = COS( REAL(ls, RealK) * pi /                   &
                         (2.0 * n_order_phase_solar) )
            ELSE
              smoothing = 1.0e+00_RealK
            END IF

            cnst1=1.0e+00_RealK-1.0e+00_RealK/REAL(ls, RealK)
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              p_legendre_tmp(l)=p_legendre_ls(l)
              p_legendre_ls(l)                                          &
                =(1.0e+00_RealK+cnst1)*p_legendre_ls(l)                 &
                *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
              p_legendre_ls_m1(l)=p_legendre_tmp(l)

              phf_tmp=(cloud_parameter(5*ls+7)/sz(l)**4 +               &
                cloud_parameter(5*ls+8)/sz(l)**3 +                      &
                cloud_parameter(5*ls+9)/sz(l)**2 +                      &
                cloud_parameter(5*ls+10)/sz(l) +                        &
                cloud_parameter(5*ls+11) ) *                            &
                EXP(-cloud_parameter(2)*REAL(ls)/sz(l))

              IF (ls.EQ.n_order_phase_solar) phf_tmp=0.5*phf_tmp

              ks_phf(l)=k_ext_scat_cloud(l,i)*phf_tmp
              phase_fnc_solar_cloud(l, i, id)                           &
                = phase_fnc_solar_cloud(l, i, id)                       &
                + ks_phf(l)*p_legendre_ls(l)                            &
                * REAL(2*ls+1, RealK)                                   &
                * smoothing
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

            phf_tmp=(cloud_parameter(5*ls+7)/sz(l)**4 +                 &
              cloud_parameter(5*ls+8)/sz(l)**3 +                        &
              cloud_parameter(5*ls+9)/sz(l)**2 +                        &
              cloud_parameter(5*ls+10)/sz(l) +                          &
              cloud_parameter(5*ls+11) ) *                              &
              EXP(-cloud_parameter(2)*REAL(ls)/sz(l))
!
            forward_solar_cloud(l, i)                                   &
              = k_ext_scat_cloud(l, i)*0.5*phf_tmp                      &
              * smoothing
          END DO
        END IF
      END IF
    END DO

  ELSE IF (i_parametrization_drop == ip_drop_unparametrized) THEN
    CALL prsc_opt_prop(ierr                                             &
      , n_profile, n_cloud_top, n_layer                                 &
      , l_rescale, n_order_forward                                      &
      , l_henyey_greenstein_pf, n_order_phase                           &
      , p, density                                                      &
      , n_opt_level_cloud_prsc                                          &
      , drop_pressure_prsc, drop_absorption_prsc                        &
      , drop_scattering_prsc, drop_phase_fnc_prsc                       &
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
      //'of cloud droplets has been selected.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE opt_prop_water_cloud
END MODULE opt_prop_water_cloud_mod

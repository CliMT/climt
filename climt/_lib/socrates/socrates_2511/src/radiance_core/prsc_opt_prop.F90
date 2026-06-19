! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set observational optical properties of aerosols.
!
! Method:
!   Splined fits to the given data at the corresponding pressure
!   levels are carried out. Optical properties at the required
!   pressure levels are calculated from the splined fits.
!
!- ---------------------------------------------------------------------
MODULE prsc_opt_prop_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE prsc_opt_prop(ierr                                           &
    , n_profile, i_first_layer, i_last_layer                            &
    , l_rescale, n_order_forward                                        &
    , l_henyey_greenstein_pf, n_order_phase                             &
    , p, density                                                        &
    , n_opt_level_prsc, prsc_pressure                                   &
    , prsc_absorption, prsc_scattering, prsc_phase_fnc                  &
    , k_ext_tot, k_ext_scat, phase_fnc                                  &
    , forward_scatter, forward_solar                                    &
    , l_solar_phf, n_order_phase_solar, n_direction, mu_v               &
    , phase_fnc_solar                                                   &
    , nd_profile, nd_radiance_profile, nd_layer, id_1, id_2             &
    , nd_direction                                                      &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    , nd_phase_term, nd_max_order                                       &
    )

  USE realtype_rd, ONLY: RealK
  USE prsc_gather_spline_mod, ONLY: prsc_gather_spline

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_radiance_profile                                               &
!       Size allocated for points where radiances are calculated
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_1                                                              &
!       Declared topmost layer for optical properties
    , id_2                                                              &
!       Declared bottom layer for optical properties
    , nd_profile_prsc                                                   &
!       Size allowed for profiles of prescribed properties
    , nd_opt_level_prsc                                                 &
!       Size allowed for levels of prescribed properties
    , nd_phase_term                                                     &
!       Size allowed for terms in the phase function
    , nd_max_order                                                      &
!       Size allowed orders of spherical harmonics
    , nd_direction
!       Size allowed for viewing directions

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_first_layer                                                     &
!       First layer where properties are required
    , i_last_layer                                                      &
!       Last layer where properties are required
    , n_order_forward                                                   &
!       Order used for forward scattering
    , n_order_phase                                                     &
!       Number of terms required in the phase function
    , n_order_phase_solar                                               &
!       Number of terms retained in the calculation of the
!       single scattering of solar radiation
    , n_direction
!       Number of viewing directions
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale                                                         &
!       Flag for rescaling
    , l_henyey_greenstein_pf                                            &
!       Flag to use Henyey-Greenstein phase functions
    , l_solar_phf
!       Flag to calculate the cwsingel scattering of solar
!       radiation directly
  REAL (RealK), INTENT(IN) ::                                           &
      p(nd_profile, nd_layer)                                           &
!       Pressure field
    , density(nd_profile, nd_layer)                                     &
!       Density field
    , mu_v(nd_profile, nd_direction)
!       Viewing directions

  INTEGER, INTENT(IN) ::                                                &
      n_opt_level_prsc
!       Number of levels of prescribed optical data
  REAL (RealK), INTENT(IN) ::                                           &
      prsc_pressure(nd_opt_level_prsc)                                  &
!       Pressure at prescribed levels
    , prsc_absorption(nd_profile_prsc, nd_opt_level_prsc)               &
!       Prescribed absorption
    , prsc_scattering(nd_profile_prsc, nd_opt_level_prsc)               &
!       Prescribed scattering
    , prsc_phase_fnc(nd_profile_prsc, nd_opt_level_prsc                 &
        , nd_phase_term)
!       Prescribed phase function

! Optical properties:
  REAL (RealK), INTENT(OUT) ::                                          &
      k_ext_tot(nd_profile, id_1: id_2)                                 &
!       Extinction
    , k_ext_scat(nd_profile, id_1: id_2)                                &
!       Scattering
    , phase_fnc(nd_profile, id_1: id_2, nd_max_order)                   &
!       Phase function: on exit this will be weighted by
!       the scattering.
    , forward_scatter(nd_profile, id_1: id_2)                           &
!       Forward scattering: on exit this will be weighted by
!       the scattering.
    , forward_solar(nd_profile, id_1: id_2)                             &
!       Forward scattering for the solar beam: on exit this
!       will be weighted by the scattering.
    , phase_fnc_solar(nd_radiance_profile, id_1: id_2                   &
        , nd_direction)
!       Current contribution to the solar phase function:
!       on exit this will be weighted by the scattering.


! Local variables.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , ls                                                                &
!       Loop variable
    , n_order_required
!       Order of terms which are required in the phase function

! Legendre polynomials:
  REAL (RealK) ::                                                       &
      cnst1                                                             &
!       Constant in recurrence for Legendre polynomials
    , p_legendre_ls(nd_profile)                                         &
!       Legendre polynomial at the current order
    , p_legendre_ls_m1(nd_profile)                                      &
!       Legendre polynomial at the previous order
    , p_legendre_tmp(nd_profile)
!       Temporary Legendre polynomial

! Temporary variables related to the phase function
  REAL (RealK) ::                                                       &
      asymmetry(nd_profile, id_1: id_2)                                 &
!       Asymmetry
    , phf_coeff(nd_profile, id_1: id_2)
!       Current coefficient in the phase function



! Absorption:
! Put this into the total extinction for now.
  CALL prsc_gather_spline(ierr                                          &
    , n_profile, i_first_layer, i_last_layer, p                         &
    , n_opt_level_prsc, prsc_pressure, prsc_absorption                  &
    , k_ext_tot                                                         &
    , nd_profile, nd_layer, id_1, id_2                                  &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    )
! Scattering:
  CALL prsc_gather_spline(ierr                                          &
    , n_profile, i_first_layer, i_last_layer, p                         &
    , n_opt_level_prsc, prsc_pressure, prsc_scattering                  &
    , k_ext_scat                                                        &
    , nd_profile, nd_layer, id_1, id_2                                  &
    , nd_profile_prsc, nd_opt_level_prsc                                &
    )
! Prescribed optical properties are not given as volume extinction
! coefficients, so they must be scaled by the density.
! Limit to positive values in case the spline fit returns negatives.
  DO i=i_first_layer, i_last_layer
    DO l=1, n_profile
      k_ext_scat(l, i)=MAX(0.0_RealK, k_ext_scat(l, i)/density(l, i))
      k_ext_tot(l, i)=MAX(0.0_RealK, k_ext_tot(l, i)/density(l, i))
    END DO
  END DO


! Phase function:

  IF (l_henyey_greenstein_pf) THEN

!   Interpolate only the first moment of the phase function.
    CALL prsc_gather_spline(ierr                                        &
      , n_profile, i_first_layer, i_last_layer, p                       &
      , n_opt_level_prsc, prsc_pressure, prsc_phase_fnc(1, 1, 1)        &
      , asymmetry                                                       &
      , nd_profile, nd_layer, id_1, id_2                                &
      , nd_profile_prsc, nd_opt_level_prsc                              &
      )

!   Initialize at the first order including the weighting.
    DO i=i_first_layer, i_last_layer
      DO l=1, n_profile
        phase_fnc(l, i, 1)=k_ext_scat(l, i)*asymmetry(l, i)
      END DO
    END DO

!   Expand all other moments.
    DO ls=2, n_order_phase
      DO i=i_first_layer, i_last_layer
        DO l=1, n_profile
          phase_fnc(l, i, ls)                                           &
            =phase_fnc(l, i, ls-1)*asymmetry(l, i)
        END DO
      END DO
    END DO

!   Calculate the forward scattering using special code for
!   the common cases.
    IF (l_rescale) THEN
      IF (n_order_forward == n_order_phase) THEN
        DO i=i_first_layer, i_last_layer
          DO l=1, n_profile
            forward_scatter(l, i)=phase_fnc(l, i, n_order_phase)
          END DO
        END DO
      ELSE IF (n_order_forward == n_order_phase+1) THEN
        DO i=i_first_layer, i_last_layer
          DO l=1, n_profile
            forward_scatter(l, i)=phase_fnc(l, i, n_order_phase)        &
              *asymmetry(l, i)
          END DO
        END DO
      ELSE
!       This case is unlikely so inefficient code is used.
        DO i=i_first_layer, i_last_layer
          DO l=1, n_profile
            forward_scatter(l, i)=k_ext_scat(l, i)                      &
              *asymmetry(l, i)**n_order_forward
          END DO
        END DO
      END IF
    END IF

  ELSE

    IF (l_rescale) THEN
      n_order_required=MAX(n_order_phase, n_order_forward)
    ELSE
      n_order_required=n_order_phase
    END IF

    DO ls=1, n_order_required

      CALL prsc_gather_spline(ierr                                      &
        , n_profile, i_first_layer, i_last_layer, p                     &
        , n_opt_level_prsc, prsc_pressure, prsc_phase_fnc(1, 1, ls)     &
        , phase_fnc(1, 1, ls)                                           &
        , nd_profile, nd_layer, id_1, id_2                              &
        , nd_profile_prsc, nd_opt_level_prsc                            &
        )

!     The phase function must be weighted by the scattering to
!     calculate the correct overall phase function later.
      DO i=i_first_layer, i_last_layer
        DO l=1, n_profile
          phase_fnc(l, i, ls)=k_ext_scat(l, i)*phase_fnc(l, i, ls)
        END DO
      END DO

    END DO

!   The forward scattering must also be weighted by the
!   scattering extinction, but this was done within the foregoing
!   loop so here we may simply copy the terms.
    IF (l_rescale) THEN
      DO i=i_first_layer, i_last_layer
        DO l=1, n_profile
          forward_scatter(l, i)=phase_fnc(l, i, n_order_forward)
        END DO
      END DO
    END IF

  END IF


! Higher orders of solar truncation:
  IF (l_solar_phf) THEN

!   It is somewhat inefficient to recalculate the lower orders
!   of the phase function, but the coding is simpler and the
!   penalty in practice will not be too great. To avoid using
!   excessive amounts of memory the moments of the phase function
!   are not stored for re-use with different directions.
!   Note that PHF_COEFF is the real coefficient in the phase
!   function, without any weighting by the scattering.
    DO id=1, n_direction

!     Calculate the asymmetry at all points and levels
      CALL prsc_gather_spline(ierr                                      &
        , n_profile, i_first_layer, i_last_layer, p                     &
        , n_opt_level_prsc, prsc_pressure, prsc_phase_fnc(1, 1, 1)      &
        , phf_coeff                                                     &
        , nd_profile, nd_layer, id_1, id_2                              &
        , nd_profile_prsc, nd_opt_level_prsc                            &
        )
      DO l=1, n_profile
!       Initialize the Legendre polynomials at the zeroth and
!       first orders.
        p_legendre_ls_m1(l)=1.0e+00_RealK
        p_legendre_ls(l)=mu_v(l, id)
      END DO
      DO i=i_first_layer, i_last_layer
        DO l=1, n_profile
          phase_fnc_solar(l, i, id)=1.0e+00_RealK+phf_coeff(l, i)       &
            *p_legendre_ls(l)*REAL(2*1+1, RealK)
        END DO
      END DO

      DO ls=2, n_order_phase_solar

!       Calculate the current moment of the phase function.
        IF (l_henyey_greenstein_pf) THEN
!         Calculate higher moments using the asymmetry which
!         is available from earlier computations.
          DO i=i_first_layer, i_last_layer
            DO l=1, n_profile
              phf_coeff(l, i)=phf_coeff(l, i)*asymmetry(l, i)
            END DO
          END DO
        ELSE
          CALL prsc_gather_spline(ierr                                  &
            , n_profile, i_first_layer, i_last_layer, p                 &
            , n_opt_level_prsc, prsc_pressure                           &
            , prsc_phase_fnc(1, 1, ls), phf_coeff                       &
            , nd_profile, nd_layer, id_1, id_2                          &
            , nd_profile_prsc, nd_opt_level_prsc                        &
            )
        END IF

!       Calculate higher Legendre polynomials by recurrences.
        cnst1=1.0e+00_RealK-1.0e+00_RealK/REAL(ls, RealK)
        DO l=1, n_profile
          p_legendre_tmp(l)=p_legendre_ls(l)
          p_legendre_ls(l)                                              &
            =(1.0e+00_RealK+cnst1)*p_legendre_ls(l)*mu_v(l, id)         &
            -cnst1*p_legendre_ls_m1(l)
          p_legendre_ls_m1(l)=p_legendre_tmp(l)
        END DO

        DO i=i_first_layer, i_last_layer
          DO l=1, n_profile
            phase_fnc_solar(l, i, id)=phase_fnc_solar(l, i, id)         &
              +phf_coeff(l, i)*p_legendre_ls(l)                         &
              *REAL(2*ls+1, RealK)
          END DO
        END DO

      END DO

!     Weight the phase function with the scattering extinction
!     to perform correct averaging later.
      DO i=i_first_layer, i_last_layer
        DO l=1, n_profile
          phase_fnc_solar(l, i, id)=phase_fnc_solar(l, i, id)           &
            *k_ext_scat(l, i)
        END DO
      END DO

      IF (l_rescale) THEN
!       Calculate one extra moment of the phase function to find
!       the forward scattering for theh solar beam.
        IF (l_henyey_greenstein_pf) THEN
!         Calculate higher moments using the asymmetry which
!         is available from earlier computations.
          DO i=i_first_layer, i_last_layer
            DO l=1, n_profile
              phf_coeff(l, i)=phf_coeff(l, i)*asymmetry(l, i)
            END DO
          END DO
        ELSE
          ls=n_order_phase_solar+1
          CALL prsc_gather_spline(ierr                                  &
            , n_profile, i_first_layer, i_last_layer, p                 &
            , n_opt_level_prsc, prsc_pressure                           &
            , prsc_phase_fnc(1, 1, ls), phf_coeff                       &
            , nd_profile, nd_layer, id_1, id_2                          &
            , nd_profile_prsc, nd_opt_level_prsc                        &
            )
        END IF

!       Apply the weighting by the scattering to the forward
!       scattering fraction.
        DO i=i_first_layer, i_last_layer
          DO l=1, n_profile
            forward_solar(l, i)                                         &
              =k_ext_scat(l, i)*phf_coeff(l, i)
          END DO
        END DO

      END IF

    END DO

  END IF

END SUBROUTINE prsc_opt_prop
END MODULE prsc_opt_prop_mod

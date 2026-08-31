! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set weights for the C.F. along a direction.
!
! Purpose:
!   The complementary function for the radiation involves unknown
!   coefficients: we set the weights for these coefficients in the
!   current layer here.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
MODULE set_dirn_weights_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_DIRN_WEIGHTS_MOD'
CONTAINS
SUBROUTINE set_dirn_weights(n_profile                                   &
    , ms, ls_trunc, up_lm                                               &
    , n_direction, mu_v, azim_factor                                    &
    , n_viewing_level, i_rad_layer, frac_rad_layer, i                   &
    , n_red_eigensystem, mu, eig_vec                                    &
    , isolir, z_sol, mu_0                                               &
    , l_ir_source_quad, diff_planck                                     &
    , upm_c, k_sol                                                      &
    , tau, omega, phase_fnc                                             &
    , weight_u, radiance                                                &
    , nd_profile, nd_layer, nd_direction, nd_viewing_level              &
    , nd_red_eigensystem, nd_max_order                                  &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_infra_red, ip_solar
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Dummy arguments:
! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_red_eigensystem                                                &
!       Size allocated for the reduced eigensystem
    , nd_max_order
!       Size allocated for orders of spherical harmonics

  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric profiles
    , n_direction                                                       &
!       Number of directions
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , i_rad_layer(nd_viewing_level)                                     &
!       Indices of layers containing viewing levels
    , n_red_eigensystem
!       Size of the reduced eigensystem
  INTEGER, INTENT(IN) ::                                                &
      i                                                                 &
!       Current layer
    , ms                                                                &
!       Current azimuthal order
    , ls_trunc                                                          &
!       Order of polar truncation
    , isolir
!       Index of spectral region

  REAL (RealK), INTENT(IN) ::                                           &
      mu_v(nd_profile, nd_direction)                                    &
!       Cosines of polar viewing angles
    , azim_factor(nd_profile, nd_direction)                             &
!       Azimuthal factors
    , frac_rad_layer(nd_viewing_level)                                  &
!       Fraction optical depth into its layer of the
!       viewing level
    , mu_0(nd_profile)                                                  &
!       Cosines of solar zenith angle
    , tau(nd_profile, nd_layer)                                         &
!       Optical depths
    , omega(nd_profile, nd_layer)                                       &
!       Albedos of single scattering
    , phase_fnc(nd_profile, nd_layer, nd_max_order)                     &
!       Phase function
    , mu(nd_profile, nd_red_eigensystem)                                &
!       Eigenvalues of the reduced eigensystem
    , eig_vec(nd_profile, 2*nd_red_eigensystem                          &
        , nd_red_eigensystem)                                           &
!       Eigenvalues of the full eigensystem scaled by
!       the s-parameters
    , z_sol(nd_profile, ls_trunc+1-ms)
!       Coefficient of the solar source function
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Flag for quadratic IR-sources
  REAL (RealK), INTENT(IN) ::                                           &
      diff_planck(nd_profile, nd_layer)
!       Differences in the hemispheric Planckian FLUX (bottom-top)
!       across the layer
  INTEGER, INTENT(IN) ::                                                &
      k_sol(nd_profile)
!       Indices of eigenvalues used to restore solar conditioning
  REAL (RealK), INTENT(IN) ::                                           &
      upm_c(nd_profile, 2*nd_red_eigensystem)
!       Coefficients of homogeneous solution used to restore
!       conditioning

  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_profile, nd_viewing_level, nd_direction)
!       Radiances (to be incremented by the contribution of
!       the particular integral)
  REAL (RealK), INTENT(OUT) ::                                          &
      weight_u(nd_profile, nd_viewing_level                             &
        , nd_direction, 2*nd_red_eigensystem)
!       Weights for the coefficients in the complementary
!       function


! Local variables
  INTEGER                                                               &
      l                                                                 &
!       Loop variable (points)
    , ir                                                                &
!       Loop variable (radiative levels)
    , id                                                                &
!       Loop variable (directions)
    , ls                                                                &
!       Loop variable (polar orders)
    , ll                                                                &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , ii
!       Loop variable
  INTEGER                                                               &
      n_list_up                                                         &
!       Numbers of points where the current viewing direction
!       points up
    , list_up(nd_profile)                                               &
!       List up points with where the current viewing direction
!       points up
    , n_list_down                                                       &
!       Numbers of points where the current viewing direction
!       points down
    , list_down(nd_profile)
!       List up points with where the current viewing direction
!       points up
  REAL (RealK) ::                                                       &
      geom_solar(nd_profile)                                            &
!       Geometrical factor for solar radiation
    , geom_integ_m(nd_profile)                                          &
!       Geometrical factor for negative eigenvalues
    , geom_integ_p(nd_profile)                                          &
!       Geometrical factor for positive eigenvalues
    , m_slant_depth_near(nd_profile)                                    &
!       Minus slantwise optical distance between the radiance
!       level and the nearer boundary of the current layer
    , m_slant_depth_inc(nd_profile)                                     &
!       Minus the increment in the slantwise optical distance
!       between the boundaries of the current layer or partial
!       layer when the viewing level lies within it
    , up_lm(nd_profile, nd_max_order+1, nd_direction)                   &
!       Spherical harmonics at a fixed azimuthal order
    , ls_sum_s(nd_profile)                                              &
!       Sum of terms over polar orders in the solar integral
    , ls_sum_p(nd_profile, nd_red_eigensystem)                          &
!       Sum of terms over polar orders in the integral over
!       eigenvalues
    , ls_sum_m(nd_profile, nd_red_eigensystem)                          &
!       Sum of terms over polar orders in the integral over
!       eigenvalues
    , tau_i(nd_profile)                                                 &
!       Optical depth of the relevant part of the current layer
    , frac_tau_i(nd_profile)                                            &
!       Fractional of the optical depth of the current layer in
!       the relevant part
    , trans_top(nd_profile)                                             &
!       Solar transmission from the top of the current layer to the
!       viewing level within the current layer
    , d_mu                                                              &
!       Difference in cosines of directions
    , x                                                                 &
!       Temporary variable
    , m1lsms
!       -1^(l+m)

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r                                                             &
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision
    , sq_eps_r                                                          &
!       The square root of the above
    , eta                                                               &
!       The conditioning weight
    , eta_nm
!       The conditioning multiplier applied in the numerator

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  INTEGER                       :: ierr
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SET_DIRN_WEIGHTS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(mu_0(1))
  sq_eps_r=SQRT(eps_r)

! Consider each level where radiances are required in turn
! and calculate appropriate weightings. The expressions for
! these weightings will look slightly different if the radiance
! level lies within the current layer.
  DO id=1, n_direction

!   Assemble the list of points where the direction is upward
!   or downard. Viewing along a horizontal direction is not
!   considered valid (and should be filtered out before this).
    n_list_up=0
    n_list_down=0
    DO l=1, n_profile
      IF (mu_v(l, id) >  0.0e+00_RealK) THEN
        n_list_up=n_list_up+1
        list_up(n_list_up)=l
      ELSE IF (mu_v(l, id) <  0.0e+00_RealK) THEN
        n_list_down=n_list_down+1
        list_down(n_list_down)=l
      END IF
    END DO

!   Sum the terms which vary with the polar order, but
!   are independent of the viewing level.
!   First the contributions to the particular integral:
    IF (isolir == ip_solar) THEN

      IF (ms == 0) THEN
!       The zeroth moment of the phase function is not stored
!       because it is 1, but this means that we need some
!       special code to treat the exception.
        DO l=1, n_profile
          ls_sum_s(l)=up_lm(l, 1, id)*z_sol(l, 1)
        END DO
      ELSE
        DO l=1, n_profile
          ls_sum_s(l)=phase_fnc(l, i, ms)*up_lm(l, 1, id)               &
            *z_sol(l, 1)
        END DO
      END IF
      DO ls=ms+1, ls_trunc
        DO l=1, n_profile
          ls_sum_s(l)=ls_sum_s(l)+z_sol(l, ls+1-ms)                     &
            *phase_fnc(l, i, ls)*up_lm(l, ls+1-ms, id)
        END DO
      END DO

    END IF


    DO k=1, n_red_eigensystem
      IF (ms == 0) THEN
        DO l=1, n_profile
          ls_sum_p(l, k)=up_lm(l, 1, id)*eig_vec(l, 1, k)
          ls_sum_m(l, k)=ls_sum_p(l, k)
        END DO
      ELSE
        DO l=1, n_profile
          ls_sum_p(l, k)=phase_fnc(l, i, ms)*up_lm(l, 1, id)            &
            *eig_vec(l, 1, k)
          ls_sum_m(l, k)=ls_sum_p(l, k)
        END DO
      END IF
      DO ls=ms+1, ls_trunc
        m1lsms=REAL(1-2*MOD((ls+ms), 2), RealK)
        DO l=1, n_profile
          x=phase_fnc(l, i, ls)*up_lm(l, ls+1-ms, id)                   &
            *eig_vec(l, ls+1-ms, k)
          ls_sum_p(l, k)=ls_sum_p(l, k)+x
          ls_sum_m(l, k)=ls_sum_m(l, k)+x*m1lsms
        END DO
      END DO
    END DO


    DO ir=1, n_viewing_level

!     Initialize the weights:
      DO k=1, 2*n_red_eigensystem
        DO l=1, n_profile
          weight_u(l, ir, id, k)=0.0e+00_RealK
        END DO
      END DO

!     Upward radiation:
!     No layers above the viewing level contribute here.

      IF (i >= i_rad_layer(ir)) THEN

!       Calculate minus the slantwise optical depths to the
!       boundaries of the layer. If the level where the radiance
!       is required lies in the current layer we perform the
!       calculation for a temporary layer reaching from the
!       viewing level to the bottom of the current layer.
        IF (i >  i_rad_layer(ir)) THEN
!         Full layers are required.
          DO ll=1, n_list_up
            l=list_up(ll)
            m_slant_depth_near(l)                                       &
              =(1.0e+00_RealK-frac_rad_layer(ir))                       &
              *tau(l, i_rad_layer(ir))
          END DO
          DO ii=i_rad_layer(ir)+1, i-1
            DO ll=1, n_list_up
              l=list_up(ll)
              m_slant_depth_near(l)                                     &
                =m_slant_depth_near(l)+tau(l, ii)
            END DO
          END DO
          DO ll=1, n_list_up
            l=list_up(ll)
            m_slant_depth_near(l)                                       &
              =-m_slant_depth_near(l)/mu_v(l, id)
            m_slant_depth_inc(l)=-tau(l, i)/mu_v(l, id)
            tau_i(l)=tau(l, i)
            frac_tau_i(l)=1.0e+00_RealK
          END DO
          IF (isolir == ip_solar) THEN
            DO ll=1, n_list_up
              l=list_up(ll)
              trans_top(l)=1.0e+00_RealK
            END DO
          END IF
        ELSE IF (i == i_rad_layer(ir)) THEN
!         The viewing level lies in the current layer.
          DO ll=1, n_list_up
            l=list_up(ll)
            m_slant_depth_near(l)=0.0e+00_RealK
            frac_tau_i(l)=1.0e+00_RealK-frac_rad_layer(ir)
            tau_i(l)=frac_tau_i(l)*tau(l, i)
            m_slant_depth_inc(l)=-tau_i(l)/mu_v(l, id)
          END DO
          IF (isolir == ip_solar) THEN
            DO ll=1, n_list_up
              l=list_up(ll)
              trans_top(l)                                              &
                =EXP(-frac_rad_layer(ir)*tau(l, i)/mu_0(l))
            END DO
          END IF
        END IF


        IF (isolir == ip_solar) THEN
!         Set the geometrical terms for the solar integral.
          DO ll=1, n_list_up
            l=list_up(ll)
            geom_solar(l)=(mu_0(l)/(mu_0(l)+mu_v(l, id)))               &
              *EXP(m_slant_depth_near(l))*(1.0e+00_RealK                &
              -EXP(m_slant_depth_inc(l)-tau_i(l)/mu_0(l)))
          END DO
!         Add the contribution of the particular integral to the
!         radiance. TRANS_TOP is required to adjust the solar
!         particular integral from its value at the top of the
!         actual layer to its value at the top of the notional
!         layer when the viewing level lies within the current
!         layer.
          DO ll=1, n_list_up
            l=list_up(ll)
            radiance(l, ir, id)                                         &
              =radiance(l, ir, id)+azim_factor(l, id)                   &
              *ls_sum_s(l)*omega(l, i)*geom_solar(l)*trans_top(l)
          END DO
        ELSE IF (isolir == ip_infra_red) THEN
          IF (ms == 0) THEN
            IF (l_ir_source_quad) THEN
              cmessage = 'Quadratic variation of the source function '  &
                //'not implemented for radiances.'
              ierr=i_err_fatal
              CALL ereport(RoutineName, ierr, cmessage)
            ELSE
              DO ll=1, n_list_up
                l=list_up(ll)
!               The azimuthal factor is omitted since it will be 1.
!               Numerical ill-conditioning can arise when the
!               optical depth is small, necessitating special
!               treatment.
                IF (m_slant_depth_inc(l) <  -sq_eps_r) THEN
                  x=EXP(m_slant_depth_near(l))                          &
                  *(1.0e+00_RealK-EXP(m_slant_depth_inc(l)))            &
                  /m_slant_depth_inc(l)
                ELSE
!                 Keep the first couple of terms from a power
!                 series.
                  x=-EXP(m_slant_depth_near(l))*(1.0e+00_RealK          &
                    +0.5e+00_RealK*m_slant_depth_inc(l))
                END IF
                radiance(l, ir, id)                                     &
                  =radiance(l, ir, id)                                  &
                  -(diff_planck(l, i)/pi)*x*frac_tau_i(l)               &
                  /(1.0_RealK-omega(l, i)*phase_fnc(l, i, 1))
              END DO
            END IF
          END IF
        END IF

!       Determine the contribution from each eigenvalue.
        DO k=1, n_red_eigensystem
          DO ll=1, n_list_up
            l=list_up(ll)
!           The term for u^+:
!           This may exhibit ill-conditioning, so it is perturbed
!           using L'Hopital's rule (actually we keep two terms in
!           the expansion).
            d_mu=mu(l, k)-mu_v(l, id)
            eta=eps_r/(d_mu+SIGN(sq_eps_r, d_mu))
            x=tau_i(l)/(mu(l, k)*mu_v(l, id))
            eta_nm=1.0_RealK-eta*x*(1.0_RealK+0.5_RealK*x*d_mu)
            geom_integ_p(l)                                             &
              =(mu(l, k)/(d_mu+eta))*EXP(m_slant_depth_near(l))         &
              *(EXP(-tau_i(l)/mu(l, k))*eta_nm                          &
              -EXP(m_slant_depth_inc(l)))
            geom_integ_m(l)                                             &
              =(mu(l, k)/(mu(l, k)+mu_v(l, id)))                        &
              *EXP(m_slant_depth_near(l))                               &
              *(EXP(-(tau(l, i)-tau_i(l))/mu(l, k))                     &
              -EXP(m_slant_depth_inc(l)-tau(l, i)/mu(l, k)))
          END DO

!         Combine to form the weights for each element of the
!         solution. Only a poprtion of WEIGHT_U is passed to this
!         routine, so the offsetting over layers takes
!         care of itself.
          DO ll=1, n_list_up
            l=list_up(ll)
            weight_u(l, ir, id, k)                                      &
              =omega(l, i)*ls_sum_m(l, k)*geom_integ_m(l)
            weight_u(l, ir, id, k+n_red_eigensystem)                    &
              =omega(l, i)*ls_sum_p(l, k)*geom_integ_p(l)
          END DO
        END DO

      END IF


!     Downward Radiation:
!     No layers below the current viewing level contribute here.
      IF (i <= i_rad_layer(ir)) THEN

!       Calculate the slantwise optical depths to the
!       boundaries of the layer. If the observing level lies
!       within the current layer we perform the calculation for
!       a layer reaching from the top of the current layer to
!       the observing level.
        IF (i <  i_rad_layer(ir)) THEN
          DO ll=1, n_list_down
            l=list_down(ll)
            m_slant_depth_near(l)                                       &
              =frac_rad_layer(ir)*tau(l, i_rad_layer(ir))
          END DO
          DO ii=i_rad_layer(ir)-1, i+1, -1
            DO ll=1, n_list_down
              l=list_down(ll)
              m_slant_depth_near(l)                                     &
                =m_slant_depth_near(l)+tau(l, ii)
            END DO
          END DO
          DO ll=1, n_list_down
            l=list_down(ll)
            m_slant_depth_near(l)                                       &
              =m_slant_depth_near(l)/mu_v(l, id)
            m_slant_depth_inc(l)=tau(l, i)/mu_v(l, id)
            tau_i(l)=tau(l, i)
            frac_tau_i(l)=1.0e+00_RealK
          END DO
        ELSE
!         The viewing level lies in the current layer.
          DO ll=1, n_list_down
            l=list_down(ll)
            tau_i(l)=frac_rad_layer(ir)*tau(l, i)
            m_slant_depth_near(l)=0.0e+00_RealK
            m_slant_depth_inc(l)=tau_i(l)/mu_v(l, id)
            frac_tau_i(l)=frac_rad_layer(ir)
          END DO
        END IF


        IF (isolir == ip_solar) THEN
!         Set the geometrical terms for the solar integral.
          DO ll=1, n_list_down
            l=list_down(ll)
!           This may exhibit ill-conditioning, so it is perturbed
!           using L'Hopital's rule (actually we keep two terms in
!           the expansion).
            d_mu=mu_0(l)+mu_v(l, id)
            eta=eps_r/(d_mu+SIGN(sq_eps_r, d_mu))
            x=tau_i(l)/(mu_0(l)*mu_v(l, id))
            eta_nm=1.0_RealK-eta*x*(1.0_RealK+0.5_RealK*x*d_mu)
            geom_solar(l)=(mu_0(l)/(d_mu+eta))                          &
              *EXP(m_slant_depth_near(l))                               &
              *(EXP(-tau_i(l)/mu_0(l))*eta_nm                           &
              -EXP(m_slant_depth_inc(l)))
          END DO
!         Add the contribution of the particular integral to the
!         radiance. In this case there is no factor representing
!         transmission from the top of the layer, since that is
!         intrinsically 1.
          DO ll=1, n_list_down
            l=list_down(ll)
            radiance(l, ir, id)                                         &
              =radiance(l, ir, id)+azim_factor(l, id)                   &
              *ls_sum_s(l)*omega(l, i)*geom_solar(l)
          END DO
        ELSE IF (isolir == ip_infra_red) THEN
          IF (ms == 0) THEN
            IF (l_ir_source_quad) THEN
              cmessage = 'Quadratic variation of the source function '  &
                //'not implemented for radiances.'
              ierr=i_err_fatal
              CALL ereport(RoutineName, ierr, cmessage)
            ELSE
              DO ll=1, n_list_down
                l=list_down(ll)
!               The azimuthal factor is omitted since it will be 1.
                IF (m_slant_depth_inc(l) <  -sq_eps_r) THEN
                  x=EXP(m_slant_depth_near(l))                          &
                  *(1.0e+00_RealK-EXP(m_slant_depth_inc(l)))            &
                  /m_slant_depth_inc(l)
                ELSE
!                 Keep the first couple of terms from a power
!                 series.
                  x=-EXP(m_slant_depth_near(l))*(1.0e+00_RealK          &
                    +0.5e+00_RealK*m_slant_depth_inc(l))
                END IF
                radiance(l, ir, id)                                     &
                  =radiance(l, ir, id)                                  &
                  +(diff_planck(l, i)/pi)*x*frac_tau_i(l)               &
                  /(1.0_RealK-omega(l, i)*phase_fnc(l, i, 1))
              END DO
            END IF
          END IF
        END IF

!       Determine the contribution from each eigenvalue.
        DO k=1, n_red_eigensystem
!         The term for u^+:
          DO ll=1, n_list_down
            l=list_down(ll)
            geom_integ_p(l)                                             &
              =(mu(l, k)/(mu(l, k)-mu_v(l, id)))                        &
              *EXP(m_slant_depth_near(l))                               &
              *(EXP(-(tau(l, i)-tau_i(l))/mu(l, k))                     &
              -EXP(m_slant_depth_inc(l)-tau(l, i)/mu(l, k)))
!           The term for u^- may exhibit ill-conditioning,
!           so it is perturbed using L'Hopital's rule
!           (actually we keep two terms in the expansion).
            d_mu=mu(l, k)+mu_v(l, id)
            eta=eps_r/(d_mu+SIGN(sq_eps_r, d_mu))
            x=tau_i(l)/(mu(l, k)*mu_v(l, id))
            eta_nm=1.0_RealK-eta*x*(1.0_RealK+0.5_RealK*x*d_mu)
            geom_integ_m(l)                                             &
              =(mu(l, k)/(d_mu+eta))                                    &
              *EXP(m_slant_depth_near(l))                               &
              *(EXP(-tau_i(l)/mu(l, k))*eta_nm                          &
              -EXP(m_slant_depth_inc(l)))
          END DO

!         Combine to form the weights for each element of the
!         solution.
          DO ll=1, n_list_down
            l=list_down(ll)
            weight_u(l, ir, id, k)                                      &
              =omega(l, i)*ls_sum_m(l, k)*geom_integ_m(l)
            weight_u(l, ir, id, k+n_red_eigensystem)                    &
              =omega(l, i)*ls_sum_p(l, k)*geom_integ_p(l)
          END DO

        END DO

      END IF

!     Add on the contribution from the conditioning homogeneous
!     solution. This includes some redundant calculation when
!     the weights will be zero for layers which cannot contribute.
!     Eventually, it may be better to tidy this up.
      IF (isolir == ip_solar) THEN

        DO l=1, n_profile
          radiance(l, ir, id)=radiance(l, ir, id)                       &
            +azim_factor(l, id)                                         &
            *weight_u(l, ir, id, k_sol(l))*upm_c(l, k_sol(l))
        END DO

      ELSE IF (isolir == ip_infra_red) THEN

        DO k=1, n_red_eigensystem
          DO l=1, n_profile
            radiance(l, ir, id)=radiance(l, ir, id)                     &
              +azim_factor(l, id)                                       &
              *(weight_u(l, ir, id, k)*upm_c(l, k)                      &
              +weight_u(l, ir, id, k+n_red_eigensystem)                 &
              *upm_c(l, k+n_red_eigensystem))
          END DO
        END DO

      END IF

    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_dirn_weights
END MODULE set_dirn_weights_mod

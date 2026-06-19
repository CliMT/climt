! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the particular integral for one layer
!
! Purpose:
!   This routine calculates the particular integral in the
!   requested spectral region for the current layer.
!
! Method:
!   The solar particular integral is calculated using a recurrence,
!   while the particular integral in the infra-red is calculated
!   from the Planckian terms. A complementary function is added to
!   the naive form of the particular integral to maintain
!   numerical conditioning.
!
!- ---------------------------------------------------------------------
MODULE layer_part_integ_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LAYER_PART_INTEG_MOD'
CONTAINS
SUBROUTINE layer_part_integ(                                            &
!                 Basic sizes
      n_profile, ls_trunc, ms, n_red_eigensystem                        &
!                 Numerical arrays of spherical terms
    , cg_coeff, mu, eig_vec, theta                                      &
!                 Solar variables
    , isolir, i_direct_top, mu_0, uplm_sol                              &
!                 Infra-red variables
    , diff_planck, l_ir_source_quad, diff_planck_2                      &
!                 Optical properies
    , tau, sqs2                                                         &
!                 Output variables
    , source_top, source_bottom, upm_c, k_sol, z_sol, q_0, q_1          &
!                 Dimensions
    , nd_profile, nd_max_order, nd_red_eigensystem                      &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_infra_red, ip_solar
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_red_eigensystem
!       Size allocated for the reduced eigensystem


! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric layers
    , n_red_eigensystem
!       Size of the reduced eigensystem
  INTEGER, INTENT(IN) ::                                                &
      ms                                                                &
!       Azimuthal order
    , ls_trunc
!       The truncating order of the system of equations
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Flag for spectral region
  REAL (RealK), INTENT(IN) ::                                           &
      cg_coeff(ls_trunc+1-ms)                                           &
!       Clebsch-Gordan coefficients
    , mu(nd_profile, nd_red_eigensystem)                                &
!       (Positive) Eigenvalues
    , eig_vec(nd_profile, 2*nd_red_eigensystem, nd_red_eigensystem)     &
!       Eigenvectors of the full systems for positive eigenvalues
!       (these are scaled by the s-coefficients in the routine
!       EIG_SYS)
    , theta(nd_profile, nd_red_eigensystem)
!       Exponentials of optical depths along slant paths defined
!       by the eigenvalues
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile)                                                   &
!       Optical depths of the layers
    , sqs2(nd_profile, 0: nd_max_order)
!       S-coefficients
  REAL (RealK), INTENT(IN) ::                                           &
      mu_0(nd_profile)                                                  &
!       Cosine of solar zenith angle
    , i_direct_top(nd_profile)                                          &
!       The direct solar radiance at the top of the current layer
    , uplm_sol(nd_profile, ls_trunc+2-ms)
!       Spherical harmonics of the solar direction
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Flag for quadratic source function in the IR
  REAL (RealK), INTENT(IN) ::                                           &
      diff_planck(nd_profile)                                           &
!       Differences in the hemispheric Planckian FLUX (bottom-top)
!       across the layer
    , diff_planck_2(nd_profile)
!       Twice the second differences in the hemispheric Planckian
!       FLUX

  INTEGER, INTENT(OUT) ::                                               &
      k_sol(nd_profile)
!       Index of eigenvalue used for solar conditioning
  REAL (RealK), INTENT(OUT) ::                                          &
      source_top(nd_profile, ls_trunc+1-ms)                             &
!       Source function at the top of the layer
    , source_bottom(nd_profile, ls_trunc+1-ms)                          &
!       Source function at the bottom of the layer
    , z_sol(nd_profile, ls_trunc+1-ms)                                  &
!       Coefficient of the solar particular integral
    , q_0(nd_profile)                                                   &
!       Term for thermal particular integral
    , q_1(nd_profile)                                                   &
!       Term for thermal particular integral
    , upm_c(nd_profile, 2*nd_red_eigensystem)
!       Arrays for coefficients of the complementary function
!       used to condition the particular integral


! Local variables
  INTEGER                                                               &
      lsr                                                               &
!       Reduced polar order
    , k                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      v_dif(nd_profile, ls_trunc+2-ms)                                  &
!       Difference between particular integral and eigenvector
    , gamma(nd_profile)                                                 &
!       Constant used in the solar particular integral
    , x(nd_profile)                                                     &
!       Temporary variable
    , m1ls                                                              &
!       -1 raised to the power l+m
    , eig_sep(nd_profile)                                               &
!       Separation of the eigenvalue from the cosine of the
!       solar zenith angle
    , eig_sep_tmp(nd_profile)                                           &
!       Temporary version of the above used in searching
    , eig_diff                                                          &
!       Difference between eigenvalue and cosine of zenith angle
    , const
!       A 'working' constant

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='LAYER_PART_INTEG'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  upm_c=0.0

  IF (isolir == ip_solar) THEN

!   If an eigenvalue approaches the cosine of the zenith angle
!   the solution will become ill-conditioned. We reduce this
!   ill conditioning by subtracting a multiple of the
!   eigensolution with the eigenvalue closest to the cosine
!   of the zenith angle.

!   Fine the closest eigenvalue.
    DO l=1, n_profile
      k_sol(l)=1
      eig_sep(l)=ABS(mu(l, 1)-mu_0(l))
    END DO
    DO k=1, n_red_eigensystem
      DO l=1, n_profile
        eig_sep_tmp(l)=ABS(mu(l, k)-mu_0(l))
        IF (eig_sep_tmp(l) <  eig_sep(l)) THEN
          k_sol(l)=k
          eig_sep(l)=eig_sep_tmp(l)
        END IF
      END DO
    END DO

!   Determine the particular integral for this layer.
!   Upward recurrence is stable here.
    DO l=1, n_profile

      v_dif(l, 1)=0.0e+00_RealK

      m1ls=REAL(1-2*MOD(1, 2), RealK)

      v_dif(l, 2)=( -mu_0(l)*sqs2(l, ms)*v_dif(l, 1)                    &
                 + sqs2(l, ms)*m1ls*eig_vec(l,1,k_sol(l)) )             &
        /cg_coeff(1)

    END DO
!   V_SOL is required one order beyond the truncation to
!   complete the solution.
    DO lsr=3, ls_trunc+2-ms
      DO l=1, n_profile

        m1ls=REAL(1-2*MOD(lsr-1, 2), RealK)

        v_dif(l, lsr)                                                   &
          =(-mu_0(l)*sqs2(l, lsr+ms-2)*v_dif(l, lsr-1)                  &
          -cg_coeff(lsr-2)*v_dif(l, lsr-2)                              &
          +sqs2(l, lsr+ms-2)*m1ls*eig_vec(l,lsr-1,k_sol(l)))            &
          /cg_coeff(lsr-1)

      END DO
    END DO
    DO l=1, n_profile
      gamma(l)=uplm_sol(l, ls_trunc+2-ms)/v_dif(l, ls_trunc+2-ms)
    END DO


!   Calculate the source function at the top and bottom
!   of this layer.
    DO lsr=1, ls_trunc+1-ms
      DO l=1, n_profile

        z_sol(l, lsr)=i_direct_top(l)                                   &
          *(gamma(l)*v_dif(l, lsr)-uplm_sol(l, lsr))

        source_top(l, lsr)=i_direct_top(l)                              &
          *(gamma(l)*v_dif(l, lsr)-uplm_sol(l, lsr))


        IF( eig_sep(l) <  1.0e-06) THEN


          eig_diff=tau(l)/(mu_0(l)*mu(l,k_sol(l)))

          const = - gamma(l)*EXP(-tau(l)/mu_0(l))                       &
                  *( eig_diff+0.5*eig_diff**2                           &
                  *(mu(l,k_sol(l))-mu_0(l)))

        ELSE

          eig_diff=tau(l)*(1.0/mu(l,k_sol(l))-1.0/mu_0(l))

          IF (eig_diff <  0.0 ) THEN

             const=gamma(l)*EXP(-tau(l)/mu(l,k_sol(l)))                 &
                  *(EXP(eig_diff)-1.0)                                  &
                  /(mu(l,k_sol(l))-mu_0(l))

          ELSE

             const=gamma(l)*EXP(-tau(l)/mu_0(l))                        &
                  *(1.0-EXP(-eig_diff))                                 &
                  /(mu(l,k_sol(l))-mu_0(l))

          END IF
        END IF

        m1ls=REAL(1-2*MOD(lsr-1, 2), RealK)

        source_bottom(l, lsr)                                           &
          = i_direct_top(l) * EXP(-tau(l)/mu_0(l))                      &
          *( gamma(l)*v_dif(l, lsr)-uplm_sol(l, lsr))                   &
          + i_direct_top(l) *eig_vec(l, lsr, k_sol(l))                  &
          *m1ls*const

      END DO
    END DO


  ELSE IF (isolir == ip_infra_red) THEN

!   The variation of the Planckian across the layer can be either
!   linear or quadratic in the optical depth. The particular
!   integrals tend to infinity as the optical depth tends to 0, so
!   a particular form of the complementary function must be added
!   to cancel off the singularity; otherwise ill-conditioning will
!   arise. Since the Planckian is azimuthally symmetric only terms
!   with m=0 are affected. Linear variations in the Planckian
!   produce a term in the particular integral with polar order 1.
!   More complicated variations produce terms at higher orders.
!   Note that ill-conditioning has been removed only in the case of
!   linear variations so far as the quadratic case is more
!   complicated. To deal with the case when TAU is 0, we add a
!   tolerance: it is therefore essential that Q_0 should be used
!   consistently to avoid errors in this limit.

    IF (ms == 0) THEN

      DO l=1, n_profile
        q_0(l)=SQRT(4.0e+00_RealK/(3.0e+00_RealK*pi))                   &
          *diff_planck(l)/(sqs2(l, 1)*tau(l)+EPSILON(tau))
      END DO

      IF (l_ir_source_quad) THEN

        DO l=1, n_profile
          q_1(l)=2.0e+00_RealK                                          &
            *SQRT(4.0e+00_RealK/(3.0e+00_RealK*pi))                     &
            *diff_planck_2(l)/(sqs2(l, 1)*tau(l)**2+EPSILON(tau))
          source_top(l, 1)                                              &
            =cg_coeff(1)*q_1(l)/sqs2(l, 0)
          source_bottom(l, 1)=source_top(l, 1)
          source_top(l, 2)=q_0(l)-0.5e+00_RealK*q_1(l)
          source_bottom(l, 2)=q_0(l)+0.5e+00_RealK*q_1(l)
        END DO
        IF (ls_trunc >  1) THEN
          DO l=1, n_profile
            source_top(l, 3)=cg_coeff(2)*q_1(l)/sqs2(l, 2)
            source_bottom(l, 3)=source_top(l, 3)
          END DO
        END IF

      ELSE

        DO l=1, n_profile
          source_top(l, 1)=0.0e+00_RealK
          source_bottom(l, 1)=0.0e+00_RealK
          source_top(l, 2)=q_0(l)
          source_bottom(l, 2)=source_top(l, 2)
        END DO
        IF (ls_trunc >  1) THEN
          DO l=1, n_profile
            source_top(l, 3)=0.0e+00_RealK
            source_bottom(l, 3)=0.0e+00_RealK
          END DO
        END IF

      END IF

!     Higher orders are unaffected.
      DO lsr=4, ls_trunc+1-ms
        DO l=1, n_profile
          source_top(l, lsr)=0.0e+00_RealK
          source_bottom(l, lsr)=0.0e+00_RealK
        END DO
      END DO

!     Now define the part of the complementary function to
!     restore conditioning.
      DO k=1, n_red_eigensystem
        DO l=1, n_profile
          upm_c(l, k+n_red_eigensystem)                                 &
            =-q_0(l)*sqs2(l, 1)*eig_vec(l, 2, k)
          upm_c(l, k)=-upm_c(l, k+n_red_eigensystem)
        END DO
      END DO

!     We take advantage of the relationship between the formats
!     of the positive and negative exponentials to reduce the
!     number of operations.
      DO lsr=1, ls_trunc+1-ms
        m1ls=REAL(1-2*MOD(lsr-1, 2), RealK)
        DO l=1, n_profile
          x(l)=upm_c(l, 1+n_red_eigensystem)*(theta(l, 1)-m1ls)         &
            *eig_vec(l, lsr, 1)
        END DO
        DO k=2, n_red_eigensystem
          DO l=1, n_profile
            x(l)=x(l)                                                   &
              +upm_c(l, k+n_red_eigensystem)*(theta(l, k)-m1ls)         &
              *eig_vec(l, lsr, k)
          END DO
        END DO
        DO l=1, n_profile
          source_top(l, lsr)=source_top(l, lsr)+x(l)
          source_bottom(l, lsr)=source_bottom(l, lsr)-m1ls*x(l)
        END DO
      END DO

    ELSE

!     This code should never be executed as non-zero azimuthal
!     orders are not relevant in the IR, but it is here for
!     safety.
      DO lsr=1, ls_trunc+1-ms
        DO l=1, n_profile
          source_top(l, lsr)=0.0e+00_RealK
          source_bottom(l, lsr)=0.0e+00_RealK
        END DO
      END DO
    END IF

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE layer_part_integ
END MODULE layer_part_integ_mod

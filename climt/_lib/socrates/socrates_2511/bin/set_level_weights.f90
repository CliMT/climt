! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set level weights for calculating radiances.
!
! Purpose:
!   This routine yields the weights to be applied to the
!   solution of the equation for the complementary function.
!
! Method:
!   The particular integral is evaluated at each viewing level
!   within the current layer and weights are calculated for each
!   unknown.
!
!- ---------------------------------------------------------------------
MODULE set_level_weights_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_LEVEL_WEIGHTS_MOD'
CONTAINS
SUBROUTINE set_level_weights(i                                          &
!                 Basic sizes
    , n_profile, ls_trunc, ms, n_red_eigensystem                        &
!                 Numerical arrays of spherical terms
    , cg_coeff, mu, eig_vec                                             &
!                 Solar variables
    , isolir, z_sol, mu_0                                               &
!                 Infra-red variables
    , q_0, l_ir_source_quad, q_1                                        &
!                 Conditioning terms
    , upm_c, k_sol                                                      &
!                 Optical properies
    , tau, sqs2                                                         &
!                 Levels where radiances are calculated
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
    , l_assign, i_assign_level                                          &
!                 Output variables
    , c_ylm, weight_u                                                   &
!                 Dimensions
    , nd_profile, nd_viewing_level                                      &
    , nd_max_order                                                      &
    , nd_red_eigensystem, nd_sph_cf_weight                              &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_infra_red, ip_solar
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_viewing_level                                                  &
!       Allocated size for levels where radiances are calculated
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_red_eigensystem                                                &
!       Size allocated for the reduced eigensystem
    , nd_sph_cf_weight
!       Size allocated for entities to be weighted by the C. F.


! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric layers
    , n_red_eigensystem                                                 &
!       Size of the reduced eigensystem
    , i
!       Current layer
  INTEGER, INTENT(IN) ::                                                &
      ms                                                                &
!       Azimuthal order
    , ls_trunc
!       The truncating order of the system of equations
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Flag for spectral region
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile)
!       Optical depths of the layers
  REAL (RealK), INTENT(IN) ::                                           &
      mu_0(nd_profile)                                                  &
!       Cosine of solar zenith angle
    , z_sol(nd_profile, ls_trunc+1-ms)
!       The direct solar radiance
  REAL (RealK), INTENT(IN) ::                                           &
      q_0(nd_profile)                                                   &
!       Term for thermal particular integral
    , q_1(nd_profile)
!       Term for thermal particular integral

  INTEGER, INTENT(IN) ::                                                &
      k_sol(nd_profile)
!       Index of eigenvalue closest to the cosine of the solar
!       zenith angle
  REAL (RealK), INTENT(IN) ::                                           &
      upm_c(nd_profile, 2*nd_red_eigensystem)
!       Weights for exponentials in conditioning term

  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Flag for quadratic source function in the IR
  LOGICAL, INTENT(INOUT) ::                                             &
      l_assign
!       Controlling logical for assigning levels
  INTEGER, INTENT(INOUT) ::                                             &
      i_assign_level
!       Current level where radiances are to be assigned
  REAL (RealK), INTENT(IN) ::                                           &
      sqs2(nd_profile, 0: nd_max_order)                                 &
!       S-coefficients
    , cg_coeff(ls_trunc+1-ms)                                           &
!       Clebsch-Gordan coefficients
    , mu(nd_profile, nd_red_eigensystem)                                &
!       Eigenvaluse of the reduced system
    , eig_vec(nd_profile, 2*nd_red_eigensystem                          &
        , nd_red_eigensystem)
!       Eigenvectors of the full systems for positive eigenvalues
!       (these are scaled by the s-coefficients in the routine
!       EIG_SYS)

  INTEGER, INTENT(IN) ::                                                &
      n_viewing_level                                                   &
!       Number of levels where radiances are calculated
    , i_rad_layer(nd_viewing_level)
!       Layers in which radiances are calculated
  REAL (RealK), INTENT(IN) ::                                           &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers


  REAL (RealK), INTENT(OUT) ::                                          &
      c_ylm(nd_profile, nd_viewing_level, ls_trunc+1-ms)                &
!       Coefficients for radiances
    , weight_u(nd_profile, nd_viewing_level, nd_sph_cf_weight           &
        , 2*nd_red_eigensystem)
!       Weights to be applied to the vector U containing the
!       complementary functions


! Local variables
  INTEGER                                                               &
      ivm                                                               &
!       Index for u^-
    , ivp
!       Index for u^+
  INTEGER                                                               &
      lsr                                                               &
!       Reduced polar order
    , k                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) :: m1ls
!       -1^(l+m)

  REAL (RealK) ::                                                       &
      exp_minus(nd_profile, nd_red_eigensystem)                         &
!       Exponentials on viewing levels for negative terms
    , exp_plus(nd_profile, nd_red_eigensystem)                          &
!       Exponentials on viewing levels for positive terms
    , x_m(nd_profile)                                                   &
!       Work array connected with negative exponentials
    , x_p(nd_profile)
!       Work array connected with positive exponentials

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_LEVEL_WEIGHTS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO WHILE (l_assign)

!   Calculate exponential terms for the whole routine for speed.
    DO k=1, n_red_eigensystem
      DO l=1, n_profile
        exp_minus(l, k)                                                 &
          =EXP(-frac_rad_layer(i_assign_level)*tau(l)                   &
          /mu(l, k))
        exp_plus(l, k)                                                  &
          =EXP((frac_rad_layer(i_assign_level)-1.0e+00_RealK)           &
          *tau(l)/mu(l, k))
      END DO
    END DO

!   Add on the particular integral.
    IF (isolir == ip_solar) THEN
      DO l=1, n_profile
        x_m(l)                                                          &
          =EXP(-frac_rad_layer(i_assign_level)*tau(l)/mu_0(l))
      END DO
      DO lsr=1, ls_trunc-ms+1
        DO l=1, n_profile
          c_ylm(l, i_assign_level, lsr)                                 &
            =c_ylm(l, i_assign_level, lsr)+x_m(l)*z_sol(l, lsr)
        END DO
      END DO

!     Add on the homogeneous conditioning solution.
      DO l=1, n_profile
        x_m(l)=upm_c(l, k_sol(l))*exp_minus(l, k_sol(l))
      END DO
      DO lsr=1, ls_trunc-ms+1
        m1ls=REAL(1-2*MOD((lsr-1),2), RealK)
        DO l=1, n_profile
          c_ylm(l, i_assign_level, lsr)                                 &
            =c_ylm(l, i_assign_level, lsr)+x_m(l)*m1ls                  &
            *eig_vec(l, lsr, k_sol(l))
        END DO
      END DO

    ELSE IF (isolir == ip_infra_red) THEN

      IF (ms == 0) THEN

        IF (l_ir_source_quad) THEN

          DO l=1, n_profile
            c_ylm(l, i_assign_level, 1)                                 &
              =c_ylm(l, i_assign_level, 1)                              &
              +cg_coeff(1)*q_1(l)/sqs2(l, 0)
            c_ylm(l, i_assign_level, 2)                                 &
              =c_ylm(l, i_assign_level, 2)                              &
              +q_0(l)+q_1(l)                                            &
              *(frac_rad_layer(i_assign_level)-0.5e+00_RealK)
          END DO
          IF (ls_trunc >  1) THEN
            DO l=1, n_profile
              c_ylm(l, i_assign_level, 3)                               &
                =c_ylm(l, i_assign_level, 3)                            &
                *cg_coeff(2)*q_1(l)/sqs2(l, 2)
            END DO
          END IF

        ELSE

          DO l=1, n_profile
            c_ylm(l, i_assign_level, 2)                                 &
              =c_ylm(l, i_assign_level, 2)+q_0(l)
          END DO

!         Now add on the homogeneous conditioning solution.
          DO k=1, n_red_eigensystem
            DO l=1, n_profile
              x_m(l)=upm_c(l, k)*exp_minus(l, k)
              x_p(l)=upm_c(l, k+n_red_eigensystem)*exp_plus(l, k)
            END DO
            DO lsr=1, ls_trunc+1-ms
              m1ls=REAL(1-2*MOD(lsr-1, 2), RealK)
!             Increment subsequent terms.
              DO l=1, n_profile
                c_ylm(l, i_assign_level, lsr)                           &
                  =c_ylm(l, i_assign_level, lsr)                        &
                  +(x_m(l)*m1ls+x_p(l))*eig_vec(l, lsr, k)
              END DO
            END DO
          END DO
        END IF
      END IF

    END IF

!   Calculate the appropriate weights.
    DO k=1, n_red_eigensystem
!     Variable numbers:
      ivm=k
      ivp=ivm+n_red_eigensystem
      DO lsr=1, ls_trunc+1-ms
        m1ls=REAL(1-2*MOD((lsr-1),2), RealK)
        DO l=1, n_profile
          weight_u(l, i_assign_level, lsr, ivm)                         &
            =eig_vec(l, lsr, k)*exp_minus(l, k)*m1ls
          weight_u(l, i_assign_level, lsr, ivp)                         &
            =eig_vec(l, lsr, k)*exp_plus(l, k)
        END DO
      END DO
    END DO

!   Increment the level for assignments:
    i_assign_level=i_assign_level+1
    IF (i_assign_level <= n_viewing_level) THEN
      IF (i_rad_layer(i_assign_level) >  i) l_assign=.FALSE.
    ELSE
      l_assign=.FALSE.
    END IF

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_level_weights
END MODULE set_level_weights_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define residuals for fitting ice parametrizations.
!
MODULE cloud_fitting
!
! Description:
!   This module contains the functions to calculate the residuals
! and their gradients for cloud fitting. This module is written
! in terms of the scaled variables used for properly conditioned
! fitting and not the raw variables used in the main code, so it
! is not directly interchangeable within anything there.
!
! N.B. 15/08/03: I'm not happy with the results from the Ackerman-
! Stephens scheme in the IR. Perhaps it's only suitable for the SW?
!
!
!
!
CONTAINS
!
!
!
  FUNCTION cloud_residual( &
    species, d, actual, &
    i_fit, property, parm) &
!
    RESULT(residual)
!
!
!   Modules used:
    USE realtype_rd
    USE rad_pcf
    USE def_std_io_icf
!
!
    IMPLICIT NONE
!
!
!
    INTEGER :: i_fit
!     Fitting scheme selected
    CHARACTER  (LEN=5), Intent(IN) :: species
!     Phase of condensate to be fitted
    CHARACTER  (LEN=10), Intent(IN) :: property
!     Name of the optical propety to be fitted: this varies from
!     scheme to scheme
    REAL  (RealK), Intent(IN), Dimension(:) :: actual
!     Actual values of quantities to be fitted
    REAL  (RealK), Intent(IN), Dimension(:) :: d
!     Sizes of particles
    REAL  (RealK), Intent(IN), Dimension(:) :: parm
!     Parameters in fitting function
!
    REAL  (RealK) :: residual
!     Value of residual
!
!   Local variables
    REAL  (RealK), Allocatable, Dimension(:) :: fit
!     Calculated fitted quantities
!
!
!
!
    ALLOCATE(fit(SIZE(d)))
!
    CALL cloud_fit_eval(species, d, i_fit, property, parm, fit) 
!
    residual = SUM( (actual - fit)**2 ) / &
      REAL(SIZE(actual), RealK)
!
    DEALLOCATE(fit)
!
!
!
  END FUNCTION cloud_residual
!
!
!
  SUBROUTINE cloud_fit_eval( &
    species, d, &
    i_fit, property, parm, fit)

!   Modules used:
    USE realtype_rd
    USE rad_pcf
    USE def_std_io_icf
    USE polynomialroots, ONLY: CubicRoots, QuadraticRoots
!
    IMPLICIT NONE
!
!
!
    INTEGER :: i_fit
!     Fitting scheme selected
    CHARACTER  (LEN=5), Intent(IN) :: species
!     Phase of condensate to be fitted
    CHARACTER  (LEN=10), Intent(IN) :: property
!     Name of the optical propety to be fitted: this varies from
!     scheme to scheme
    REAL  (RealK), Intent(IN), Dimension(:) :: d
!     Sizes of particles
    REAL  (RealK), Intent(IN), Dimension(:) :: parm
!     Parameters in fitting function
!
    REAL  (RealK), Intent(OUT), Dimension(:) :: fit
!     Calculated fitted quantities
!
!   Local variables
    INTEGER :: i
!     Loop variable
    INTEGER :: j
!     Loop variable
    INTEGER, Allocatable, Dimension(:) :: offset
!     Offset

    INTEGER, PARAMETER           :: dp=KIND(1.0D0)
    INTEGER                      :: i_root
    REAL (dp),    Dimension(0:3) :: poly
    COMPLEX (dp), Dimension(3)   :: roots
!     Variables used in the calculation of cubic and quadratic roots


    IceWater: SELECT CASE(species)
!
      CASE("Water")
!
        Water_Fit: SELECT CASE(i_fit)
!
!
          CASE (IP_Slingo_Schrecker)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
              CASE ("Coalbedo  ")
                fit = parm(1) + parm(2) * d
              CASE ("Asymmetry ")
                fit = parm(1) + parm(2) * d
            END SELECT
!
          CASE (IP_Ackerman_Stephens)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) * EXP( parm(3) * LOG(d) )
              CASE ("Coalbedo  ")
                fit = parm(1) + parm(2) * EXP( parm(3) * LOG(d) )
              CASE ("Asymmetry ")
                fit = parm(1) + parm(2) * EXP( parm(3) * LOG(d) )
            END SELECT
!
          CASE (IP_Drop_Pade_2, IP_Drop_Pade_2_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = (parm(1) + d * (parm(2) + d * parm(3)) ) / &
                  (1.0_RealK + d * (parm(4) + d * (parm(5) + &
                  d * parm(6))))
                ! Check to see if the fit would give NaNs for any intermediate
                ! values of d between the fitted values. This will happen if
                ! the cubic equation in the denominator has roots in the fitted
                ! range of d.
                poly(0)   = REAL(1.0, dp)
                poly(1:3) = REAL(parm(4:6), dp)
                CALL CubicRoots(poly, roots)
                DO i_root = 1, 3
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d)) THEN
                      ! If this fit would produce NaNs then artificially
                      ! increase the value of the fit for the point closest
                      ! to the root. This will then give a large residual and
                      ! the fit should not converge to these values.
                      fit(MINLOC(ABS(d-roots(i_root)))) = &
                        fit(MINLOC(ABS(d-roots(i_root)))) / EPSILON(1.0_RealK)
                    END IF
                  END IF 
                END DO
              CASE ("Coalbedo  ")
                fit = (parm(1) + d * (parm(2) + d * parm(3)) ) / &
                  (1.0_RealK + d * (parm(4) + d * parm(5)))
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d)) THEN
                      fit(MINLOC(ABS(d-roots(i_root)))) = &
                        fit(MINLOC(ABS(d-roots(i_root)))) / EPSILON(1.0_RealK)
                    END IF
                  END IF 
                END DO
              CASE ("Asymmetry ", "Moment    ")
                fit = (parm(1) + d * (parm(2) + d * parm(3)) ) / &
                  (1.0_RealK + d * (parm(4) + d * parm(5)))
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d)) THEN
                      fit(MINLOC(ABS(d-roots(i_root)))) = &
                        fit(MINLOC(ABS(d-roots(i_root)))) / EPSILON(1.0_RealK)
                    END IF
                  END IF 
                END DO
            END SELECT
!
          CASE (IP_Slingo_Schr_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
              CASE ("Coalbedo  ")
                fit = parm(1) + parm(2) * d
              CASE ("Moment    ")
                fit = parm(1) + parm(2) * d
            END SELECT
!
          CASE DEFAULT
!
            WRITE(iu_err, '(/A)') &
              '*** Error in calculating water parametrizations.'
!
        END SELECT Water_Fit
!
!
      CASE("Ice")
!
        Ice_Fit: SELECT CASE(i_fit)
!
!
          CASE (IP_Slingo_Schrecker_ice)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
              CASE ("Coalbedo  ")
                fit = parm(1) + parm(2) * d
              CASE ("Asymmetry ")
                fit = parm(1) + parm(2) * d
            END SELECT
!
!
          CASE (IP_ice_adt)
!
            ALLOCATE(offset(SIZE(d)))
            WHERE (d > 0)
              offset = 4
            ELSEWHERE
              offset = 0
            ENDWHERE
            SELECT CASE(property)
              CASE ("Extinction", &
                    "Coalbedo  ", &
                    "Asymmetry ")
                DO i =1, SIZE(d)
                  fit(i) = parm(5+offset(i))
                ENDDO
                DO j=4, 2, -1
                  DO i =1, SIZE(d)
                    fit(i) = parm(j+offset(i)) + d(i) * fit(i)
                  ENDDO
                ENDDO
                fit = parm(1) + d * fit
            END SELECT
            DEALLOCATE(offset)
!
          CASE (IP_ice_adt_10)
!
            SELECT CASE(property)
              CASE ("Extinction", &
                    "Coalbedo  ", &
                    "Asymmetry ")
                fit = parm(11)
                DO j=10, 1, -1
                  fit = parm(j) + d * fit
                ENDDO
            END SELECT
!
          CASE (IP_Slingo_Schr_ice_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
              CASE ("Coalbedo  " , &
                    "Moment   ")
                fit = parm(1) + parm(2) * d
            END SELECT
!
          CASE (IP_ice_Fu_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(3)
                DO j=2, 1, -1
                  fit = parm(j) + fit / d
                ENDDO
              CASE ("Coalbedo  " , &
                    "Moment   ")
                fit = parm(4)
                DO j=3, 1, -1
                  fit = parm(j) + d * fit
                ENDDO
            END SELECT
!
          CASE (IP_Ice_Pade_2_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = (parm(1) + d * (parm(2) + d * parm(3)) ) / &
                  (1.0_RealK + d * (parm(4) + d * (parm(5) + &
                  d * parm(6))))
                ! Check to see if the fit would give NaNs for any intermediate
                ! values of d between the fitted values. This will happen if
                ! the cubic equation in the denominator has roots in the fitted
                ! range of d.
                poly(0)   = REAL(1.0, dp)
                poly(1:3) = REAL(parm(4:6), dp)
                CALL CubicRoots(poly, roots)
                DO i_root = 1, 3
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d)) THEN
                      ! If this fit would produce NaNs then artificially
                      ! increase the value of the fit for the point closest
                      ! to the root. This will then give a large residual and
                      ! the fit should not converge to these values.
                      fit(MINLOC(ABS(d-roots(i_root)))) = &
                        fit(MINLOC(ABS(d-roots(i_root)))) / EPSILON(1.0_RealK)
                    END IF
                  END IF 
                END DO
              CASE ("Coalbedo  ")
                fit = (parm(1) + d * (parm(2) + d * parm(3)) ) / &
                  (1.0_RealK + d * (parm(4) + d * parm(5)))
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d)) THEN
                      fit(MINLOC(ABS(d-roots(i_root)))) = &
                        fit(MINLOC(ABS(d-roots(i_root)))) / EPSILON(1.0_RealK)
                    END IF
                  END IF 
                END DO
              CASE ("Asymmetry ", "Moment    ")
                fit = (parm(1) + d * (parm(2) + d * parm(3)) ) / &
                  (1.0_RealK + d * (parm(4) + d * parm(5)))
                poly(0)   = REAL(1.0, dp)
                poly(1:2) = REAL(parm(4:5), dp)
                CALL QuadraticRoots(poly, roots)
                DO i_root = 1, 2
                  IF (ABS(AIMAG(roots(i_root))) < TINY(1.0d0)) THEN
                    IF (DBLE(roots(i_root)) > MINVAL(d) .AND. &
                        DBLE(roots(i_root)) < MAXVAL(d)) THEN
                      fit(MINLOC(ABS(d-roots(i_root)))) = &
                        fit(MINLOC(ABS(d-roots(i_root)))) / EPSILON(1.0_RealK)
                    END IF
                  END IF 
                END DO
            END SELECT
!
          CASE DEFAULT
!
            fit = 0.0_RealK
            WRITE(iu_err, '(/A)') &
              '*** Error in calculating ice parametrizations.'
!
        END SELECT Ice_Fit
!
    END SELECT IceWater
!
!
!
!
  END SUBROUTINE cloud_fit_eval
!
!
!
  SUBROUTINE cloud_residual_gradient( &
    species, d, actual, i_fit, property, parm, gradient)

!   Modules used:
    USE realtype_rd
    USE rad_pcf
    USE def_std_io_icf
!
!
    IMPLICIT NONE
!
!
!
    INTEGER :: i_fit
!     Fitting scheme selected
    CHARACTER  (LEN=5), Intent(IN) :: species
!     Phase of condensate to be fitted
    CHARACTER  (LEN=10), Intent(IN) :: property
!     Name of the optical propety to be fitted: this varies from
!     scheme to scheme
    REAL  (RealK), Intent(IN), Dimension(:) :: actual
!     Actual values of quantities to be fitted
    REAL  (RealK), Intent(IN), Dimension(:) :: d
!     Sizes of particles
    REAL  (RealK), Intent(IN), Dimension(:) :: parm
!     Parameters in fitting function
!
    REAL  (RealK), Intent(OUT), Dimension(:) :: gradient
!     Gradient of the residual with respect to the free parameters
!
!
!
!   Local variables
    INTEGER :: n_parm_free
!     Number of free parameters 
!     (the last often holds a size for scaling)
    INTEGER :: i
!     Loop variable
    INTEGER :: j
!     Loop variable
    REAL  (RealK), Allocatable, Dimension(:) :: fit
!     Calculated fitted quantities
    INTEGER, Allocatable, Dimension(:) :: offset
!     Offset
    REAL  (RealK), Allocatable, Dimension(:) :: p
!     Temporary array
    REAL  (RealK), Allocatable, Dimension(:) :: q
!     Temporary array
    REAL  (RealK), Allocatable, Dimension(:, :) :: gradient_local
!     Gradient of fitted value with respect to parameters at individual
!     points
!
!
!
    gradient=0.0
    ALLOCATE(fit(SIZE(d)))
!
    IceWater: SELECT CASE(species)
!
      CASE("Water")
!
        Water_Fit: SELECT CASE(i_fit)
!
!
          CASE (IP_Slingo_Schrecker)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = 1.0_RealK /d
              CASE ("Coalbedo  ")
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                fit = parm(1) + parm(2) * d
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
              CASE ("Asymmetry ")
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                fit = parm(1) + parm(2) * d
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
            END SELECT
!
          CASE (IP_Ackerman_Stephens)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) * EXP( parm(3) * LOG(d) )
                n_parm_free = 3
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = EXP( parm(3) * LOG(d) )
                gradient_local(:, 3) = parm(2) * LOG(d) * &
                  EXP( parm(3) * LOG(d) )
              CASE ("Coalbedo  ")
                fit = parm(1) + parm(2) * EXP( parm(3) * LOG(d) )
                n_parm_free = 3
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = EXP( parm(3) * LOG(d) )
                gradient_local(:, 3) = parm(2) * LOG(d) * &
                  EXP( parm(3) * LOG(d) )
              CASE ("Asymmetry ")
                fit = parm(1) + parm(2) * EXP( parm(3) * LOG(d) )
                n_parm_free = 3
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = EXP( parm(3) * LOG(d) )
                gradient_local(:, 3) = parm(2) * LOG(d) * &
                  EXP( parm(3) * LOG(d) )
            END SELECT
!
          CASE (IP_Drop_Pade_2, IP_Drop_Pade_2_PHF)
!
            ALLOCATE(p(SIZE(d)))
            ALLOCATE(q(SIZE(d)))
!
            SELECT CASE(property)
              CASE ("Extinction")
                p = parm(1) + d * (parm(2) + d * parm(3))
                q = 1.0_RealK + &
                  d * (parm(4) + d * (parm(5) + parm(6) * d))
                fit = p / q
                n_parm_free = 6
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK / q
                gradient_local(:, 2) = d / q
                gradient_local(:, 3) = d**2 / q
                gradient_local(:, 4) = - fit * d / q
                gradient_local(:, 5) = - fit * d**2 / q
                gradient_local(:, 6) = - fit * d**3 / q
              CASE ("Coalbedo  ")
                p = parm(1) + d * (parm(2) + d * parm(3))
                q = 1.0_RealK + d * (parm(4) + d * parm(5))
                fit = p / q
                n_parm_free = 5
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK / q
                gradient_local(:, 2) = d / q
                gradient_local(:, 3) = d**2 / q
                gradient_local(:, 4) = - fit * d / q
                gradient_local(:, 5) = - fit * d**2 / q
              CASE ("Asymmetry ", "Moment    ")
                p = parm(1) + d * (parm(2) + d * parm(3))
                q = 1.0_RealK + d * (parm(4) + d * parm(5))
                fit = p / q
                n_parm_free = 5
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK / q
                gradient_local(:, 2) = d / q
                gradient_local(:, 3) = d**2 / q
                gradient_local(:, 4) = - fit * d / q
                gradient_local(:, 5) = - fit * d**2 / q
            END SELECT
!
            DEALLOCATE(p)
            DEALLOCATE(q)
!
          CASE (IP_Slingo_Schr_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = 1.0_RealK /d
              CASE ("Coalbedo  ")
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                fit = parm(1) + parm(2) * d
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
              CASE ("Moment    ")
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                fit = parm(1) + parm(2) * d
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
            END SELECT
!
          CASE DEFAULT
!
            WRITE(iu_err, '(/A)') &
              '*** Error in calculating water parametrizations.'
!
        END SELECT Water_Fit
!
!
      CASE("Ice  ")
!
        Ice_Fit: SELECT CASE(i_fit)
!
!
          CASE (IP_Slingo_Schrecker_ice)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = 1.0_RealK /d
              CASE ("Coalbedo  ")
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                fit = parm(1) + parm(2) * d
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
              CASE ("Asymmetry ")
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                fit = parm(1) + parm(2) * d
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
            END SELECT
!
          CASE (IP_ice_adt)
!
            ALLOCATE(offset(SIZE(d)))
            WHERE (d > 0)
              offset = 4
            ELSEWHERE
              offset = 0
            ENDWHERE
            SELECT CASE(property)
              CASE ("Extinction", &
                    "Coalbedo  ", &
                    "Asymmetry ")
                DO i=1, SIZE(d)
                  fit(i) = parm(5+offset(i))
                ENDDO
                DO j=4, 2, -1
                  DO i=1, SIZE(d)
                    fit(i) = parm(j+offset(i)) + d(i) * fit(i)
                  ENDDO
                ENDDO
                fit = parm(1) + d * fit
!
                n_parm_free = 9
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = 0.0_RealK
                gradient_local(:, 6) = 0.0_RealK
                DO i=1, SIZE(d)
                  gradient_local(i, 2+offset(i)) = d(i)
                ENDDO
                DO j=1, 5
                  gradient_local(:, j)   = gradient_local(:, j-1) * d
                  gradient_local(:, j+4) = gradient_local(:, j+3) * d
                ENDDO
            END SELECT
            DEALLOCATE(offset)
!
          CASE (IP_ice_adt_10)
!
            SELECT CASE(property)
              CASE ("Extinction", &
                    "Coalbedo  ", &
                    "Asymmetry ")
                fit = parm(11)
                DO j=10, 1, -1
                  fit = parm(j) + d * fit
                ENDDO
                n_parm_free = 11
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                DO j=2, 11
                  gradient_local(:, j) = gradient_local(:, j-1) * d
                ENDDO
            END SELECT
!
          CASE (IP_Slingo_Schr_ice_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(1) + parm(2) / d
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = 1.0_RealK /d
              CASE ("Coalbedo  " , &
                    "Moment   ")
                fit = parm(1) + parm(2) * d
                n_parm_free = 2
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
            END SELECT
!
          CASE (IP_ice_Fu_PHF)
!
            SELECT CASE(property)
              CASE ("Extinction")
                fit = parm(3)
                DO j=2, 1, -1
                  fit = parm(j) + fit / d
                ENDDO
                n_parm_free = 3
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = 1.0_RealK / d
                gradient_local(:, 3) = 1.0_RealK / (d*d)
              CASE ("Coalbedo  " , &
                    "Moment    ")
                fit = parm(4)
                DO j=3, 1, -1
                  fit = parm(j) + d * fit
                ENDDO
                n_parm_free = 4
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK
                gradient_local(:, 2) = d
                gradient_local(:, 3) = d*d
                gradient_local(:, 4) = d*d*d
            END SELECT
!
          CASE (IP_Ice_Pade_2_PHF)
!
            ALLOCATE(p(SIZE(d)))
            ALLOCATE(q(SIZE(d)))
!
            SELECT CASE(property)
              CASE ("Extinction")
                p = parm(1) + d * (parm(2) + d * parm(3))
                q = 1.0_RealK + &
                  d * (parm(4) + d * (parm(5) + parm(6) * d))
                fit = p / q
                n_parm_free = 6
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK / q
                gradient_local(:, 2) = d / q
                gradient_local(:, 3) = d**2 / q
                gradient_local(:, 4) = - fit * d / q
                gradient_local(:, 5) = - fit * d**2 / q
                gradient_local(:, 6) = - fit * d**3 / q
              CASE ("Coalbedo  ")
                p = parm(1) + d * (parm(2) + d * parm(3))
                q = 1.0_RealK + d * (parm(4) + d * parm(5))
                fit = p / q
                n_parm_free = 5
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK / q
                gradient_local(:, 2) = d / q
                gradient_local(:, 3) = d**2 / q
                gradient_local(:, 4) = - fit * d / q
                gradient_local(:, 5) = - fit * d**2 / q
              CASE ("Asymmetry ", "Moment    ")
                p = parm(1) + d * (parm(2) + d * parm(3))
                q = 1.0_RealK + d * (parm(4) + d * parm(5))
                fit = p / q
                n_parm_free = 5
                ALLOCATE(gradient_local(SIZE(d), n_parm_free))
                gradient_local(:, 1) = 1.0_RealK / q
                gradient_local(:, 2) = d / q
                gradient_local(:, 3) = d**2 / q
                gradient_local(:, 4) = - fit * d / q
                gradient_local(:, 5) = - fit * d**2 / q
            END SELECT
!
            DEALLOCATE(p)
            DEALLOCATE(q)
!
          CASE DEFAULT
!
            WRITE(iu_err, '(/A)') &
              '*** Error in calculating ice parametrizations.'
!
        END SELECT Ice_Fit
!
    END SELECT IceWater
!
    DO i=1, n_parm_free
      gradient(i) = &
        ( 2.0_RealK / REAL(SIZE(actual), RealK) ) * &
        DOT_PRODUCT( (fit - actual), gradient_local(:, i) )  
    ENDDO
!
    DEALLOCATE(fit)
    DEALLOCATE(gradient_local)
!
!
!
  END SUBROUTINE cloud_residual_gradient
!
!
!
END MODULE cloud_fitting

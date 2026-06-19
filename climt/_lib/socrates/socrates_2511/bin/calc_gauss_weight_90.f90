! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate points for Gaussian integration.
!
! Description:
!   This subroutine calculates the zeros of a Legendre polynomial
!   of the appropriate order and determines corresponding Gaussian
!   weights.
!
! Method:
!   The algorithm is not very efficient, but that is unlikely to
!   matter in the cases in which this subroutine will probably be
!   used. Zeros are determined for all orders up to that required
!   by an iterative process. The zeros of the Legendre polynomial 
!   of order n+1 lie in the intervals defined by the end points
!   -1 and 1 and the zeros of the polynomial of order n. The root
!   can then be isolated by false position.
!
!------------------------------------------------------------------------------
MODULE calc_gauss_weight_90_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE calc_gauss_weight_90(ierr, n, root, weight)

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal
  USE errormessagelength_mod, ONLY: errormessagelength
  USE ereport_mod, ONLY: ereport
  USE legendre_weight_mod, ONLY: legendre_weight

  IMPLICIT NONE


  INTEGER, INTENT(INOUT) :: ierr
!   Error flag

  INTEGER, INTENT(IN) :: n
!   Number of points required
  REAL  (RealK), INTENT(OUT) :: root(n)
!   Roots of the Legendre polynomial (quadrature points)
  REAL  (RealK), INTENT(OUT) :: weight(n)
!   Weights at quadrature points


! Local variables
  INTEGER :: order
!   Order of polynomial
  INTEGER :: n_iteration
!   Number of iterations in fiding root
  INTEGER, PARAMETER :: NP_Gauss_maxit = 100
!   Maximum number of iterations to find Gaussian weights
  INTEGER :: i
!   Loop variable
  REAL  (RealK) :: bracket(n+1)
!   Values bracketting roots
  REAL  (RealK) :: xl
!   Lower limit of interval containing the root
  REAL  (RealK) :: xh
!   Upper limit of interval containing the root
  REAL  (RealK) :: yl
!   Value of function at lower limit
  REAL  (RealK) :: yh
!   Value of function at upper limit
  REAL  (RealK) :: x
!   New trial position
  REAL  (RealK) :: y
!   Value of function at trial position
  REAL  (RealK) :: temp
!   Temporary variable used in swapping
  REAL  (RealK) :: shift
!   Absolute shift in the position of the root in this iteration
  REAL  (RealK) :: wt
!   Temporary weight
  REAL  (RealK) :: sum_weight
!   Sum of weights

! Variables related to the treatment of ill-conditioning
  REAL  (RealK) :: tol_zero
!   The tolerance for detecting zeros

  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CALC_GAUSS_WEIGHT_90'



! Set the tolerances used in avoiding ill-conditioning.
  tol_zero = 1.6e+01_RealK*EPSILON(tol_zero)

! First order:
  order = 1
  root(1) = 0.0_RealK
  weight(1) = 2.0_RealK

  DO order = 2, n

!   Set intervals bracketting the roots.
    bracket(1) = -1.0_RealK + tol_zero
    bracket(2:order) = root(1:order-1)
    bracket(order+1) = 1.0_RealK-tol_zero

    DO i = 1, order

      xl = bracket(i)
      CALL legendre_weight(xl, order, yl, wt)
      xh = bracket(i+1)
      CALL legendre_weight(xh, order, yh, wt)
      IF (yl*yh > 0.0_RealK) THEN
        cmessage = 'Root not bracketted.'
        ierr=i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      ENDIF

!     Swap so that the function is negative at XL
      IF (yl > 0.0_RealK) THEN
        temp = xh
        xh = xl
        xl = temp
        temp = yh
        yh = yl
        yl = temp
      ENDIF

!     Find the root by false position.
      n_iteration = 0
      shift = abs(xh-xl)
      y = yl

      DO WHILE ( (shift > tol_zero) .AND. (ABS(y) > EPSILON(shift)) )
        n_iteration = n_iteration + 1
        IF (n_iteration < NP_gauss_maxit) THEN
          x = xl-yl*(xh-xl)/(yh-yl)
          CALL legendre_weight(x, order, y, wt)
          IF (y < 0.0_RealK) THEN
            shift = abs(x-xl)
            xl = x
            yl = y
          ELSE
            shift = abs(xh-x)
            xh = x
            yh = y
          ENDIF
        ELSE
          cmessage = 'Too many interations in ' // &
            'calculating Legendre functions'
          ierr = i_err_fatal
          CALL ereport(RoutineName, ierr, cmessage)
        ENDIF
      ENDDO

      root(i) = x
      weight(i) = wt

    ENDDO

  ENDDO

  sum_weight = SUM(weight(1:n))

  IF (ABS(sum_weight-2.0_RealK) > 1.0e+04_RealK*EPSILON(sum_weight)) THEN
    cmessage = 'Gaussian weights do not sum correctly.'
    ierr = i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)
  ENDIF

END SUBROUTINE calc_gauss_weight_90
END MODULE calc_gauss_weight_90_mod

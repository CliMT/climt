! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate a Legendre polynomial and a weight.
!
! Purpose:
!   This subroutine uses a recurrence relation to calculate the
!   value of a Legendre polynomial and determines the weight for
!   Gaussian integration on the assumption that the polynomial
!   is actually zero.
!
!------------------------------------------------------------------------------
MODULE legendre_weight_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE legendre_weight(x, n, lpy, wt)

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE


  INTEGER, Intent(IN) :: n
!   Order of polynomial
  REAL  (RealK), Intent(IN) :: x
!   Point of evaluation
  REAL  (RealK), Intent(OUT) :: lpy
!   Value of polynomial
  REAL  (RealK), Intent(OUT) :: wt
!   Corresponding weight

! Local variables
  INTEGER ::  j
!   Loop variable 
  REAL  (RealK) :: y
!   Temporary value of polynomial
  REAL  (RealK) :: lpy_m2
!   Polynomial at order J-2
  REAL  (RealK) :: lpy_m1
!   Polynomial at order J-1
  REAL  (RealK) :: d_lpy
!   Derivative of polynomial at order N


! Initialize at the zeroth and first orders to start the recurrence.
  lpy_m2=1.0_RealK
  lpy_m1=x
  DO j=2, n
    y=(x*REAL(2*j-1, RealK)*lpy_m1-REAL(j-1, RealK)*lpy_m2) &
      /REAL(j, RealK)
    lpy_m2=lpy_m1
    lpy_m1=y
  ENDDO
  lpy=y
  d_lpy=-DBLE(n)*lpy_m2/(x*x-1.0_RealK)
  wt=2.0_RealK/MAX((REAL(n, RealK)*d_lpy*lpy_m2),TINY(d_lpy))

END SUBROUTINE legendre_weight
END MODULE legendre_weight_mod

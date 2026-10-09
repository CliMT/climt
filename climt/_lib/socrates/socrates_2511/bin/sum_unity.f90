! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to normalize the weights in an ESFT fit.
!
!- ---------------------------------------------------------------------
SUBROUTINE sum_unity(n_term, a, k)

! Modules to set types of variables:
  USE realtype_rd
  USE dimensions_spec_ucf

  IMPLICIT NONE

! Dummy arguments.
  INTEGER     , INTENT(IN)    :: n_term        ! Number of terms in fit
  REAL (RealK), INTENT(INOUT) :: a(npd_k_term) ! Weights of fit
  REAL (RealK), INTENT(INOUT) :: k(npd_k_term) ! Exponents of fit

! Local variables.
  INTEGER      :: i   ! Loop variable
  REAL (RealK) :: sum ! Sum of weights


! Find the sum of the present weights, then scale them to sum to 1.
! Adjust the exponents correspondingly so as to preserve the limit
! of small amounts of absorber.
  sum=0.0E+00_RealK
  DO i=1, n_term
    sum=sum+a(i)
  END DO
  DO i=1, n_term
    a(i)=a(i)/sum
    k(i)=k(i)*sum
  END DO

END SUBROUTINE sum_unity

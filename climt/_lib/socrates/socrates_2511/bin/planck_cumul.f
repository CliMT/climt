! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the cumulative Planck function.
!
! Method:
!  The integral of x^3/(e^x-1) from 0 to x is evaluated using the
!  formula obtained from a binomial expansion of the denoimator.
!
!- ---------------------------------------------------------------------
      FUNCTION planck_cumul(x, tol)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments.
      REAL  (RealK) ::
     &    planck_cumul
!           Name of function
      REAL  (RealK), Intent(IN) ::
     &    x
!           Argument of the function
     &  , tol
!           Relative accuracy required
!
!     Local variables.
      INTEGER
     &    n_term
!           Number of the current term
      REAL  (RealK) ::
     &    nx
!           Number of term multiplied by the argument
     &  , term
!           Current term
!
!
!
      term=0.0_RealK
      planck_cumul=0.0_RealK
      n_term=0
      DO WHILE (term >= tol*planck_cumul)
!
        n_term=n_term+1
        nx=real(n_term, RealK)*x
        term=(1.0_RealK/(real(n_term, RealK))**4)*(6.0_RealK
     &    -exp(-nx)
     &    *(nx*(nx*(nx+3.0_RealK)+6.0_RealK)+6.0_RealK))
        planck_cumul=planck_cumul+term
!
      ENDDO
!
!
!
      RETURN
      END

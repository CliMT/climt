! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve a set of banded matrix equations.
!
! Method:
!       A set of bands matrix equations is solved using the
!       standard method of Gaussian elimination. Diagonals are
!       numbered downward (i.e. upper diagonals first).
!
!- ---------------------------------------------------------------------
MODULE band_solver_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BAND_SOLVER_MOD'
CONTAINS
SUBROUTINE band_solver(n_matrix, n_equation                             &
    , iu, il                                                            &
    , a, b                                                              &
    , x                                                                 &
    , rho                                                               &
    , nd_matrix, nd_diagonal, nd_equation                               &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_matrix                                                         &
!       Size alloacted for matrices
    , nd_diagonal                                                       &
!       Size allocated for diagonals
    , nd_equation
!       Size allocated for equations

! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      n_matrix                                                          &
!       Number of matrices
    , n_equation                                                        &
!       Number of equations
    , iu                                                                &
!       Number of superdiagonals
    , il
!       Number of subdiagonals
  REAL (RealK), INTENT(INOUT) ::                                        &
      a(nd_matrix, nd_diagonal, nd_equation)                            &
!       Matrices of coefficients
    , b(nd_matrix, nd_equation)
!       Righthand sides
  REAL (RealK), INTENT(OUT) ::                                          &
       x(nd_matrix, nd_equation)
!       Solution vector
  REAL (RealK) ::                                                       &
                     !, INTENT(WORK)
       rho(nd_matrix)
!       Temporary array

! Local variables
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , iu1
!       Local scalar

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='BAND_SOLVER'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  iu1=iu+1
! Eliminative phase.
  DO i=n_equation, 2, -1
    DO j=1, MIN(iu, i-1)
      DO l=1, n_matrix
        rho(l)=a(l, iu1-j, i-j)/a(l, iu1, i)
        b(l, i-j)=b(l, i-j)-rho(l)*b(l, i)
      END DO
      DO k=1, MIN(il, i-1)
        DO l=1, n_matrix
          a(l, iu1+k-j, i-j)=a(l, iu1+k-j, i-j)                         &
            -rho(l)*a(l, iu1+k, i)
        END DO
      END DO
    END DO
  END DO

! Solution and back-substitution:

  IF ( (iu == 2).AND.(il == 2) ) THEN
!   A special version is used for the pentadiagonal case to allow
!   us to chain operations together for efficiency on the CRAY
!   vector machines, as this particular case arises quite often.

!   First equation:
    DO l=1, n_matrix
      x(l, 1)=b(l, 1)/a(l, 3, 1)
    END DO
!   Second equation:
    DO l=1, n_matrix
      x(l, 2)=(b(l, 2)-a(l, 4, 2)*x(l, 1))/a(l, 3, 2)
    END DO
!   Remaining equations:
    DO i=3, n_equation
      DO l=1, n_matrix
        x(l, i)=(b(l, i)-a(l, 4, i)*x(l, i-1)                           &
          -a(l, 5, i)*x(l, i-2))/a(l, 3, i)
      END DO
    END DO
  ELSE

!   General case:
    DO i=1, n_equation
      DO l=1, n_matrix
           x(l, i)=b(l, i)
      END DO
      DO k=1, MIN(il, i-1)
        DO l=1, n_matrix
          x(l, i)=x(l, i)-a(l, iu1+k, i)*x(l, i-k)
        END DO
      END DO
      DO l=1, n_matrix
        x(l, i)=x(l, i)/a(l, iu1, i)
      END DO
    END DO

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE band_solver
END MODULE band_solver_mod

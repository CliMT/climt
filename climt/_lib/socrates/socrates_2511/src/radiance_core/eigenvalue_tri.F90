! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to find the eigenvalues of a symmetric tridiagonal matrix.
!
! Purpose:
!   To caulate the eigenvalues of a symmetric tridiagonal matrix.
!
! Method:
!   The standard QR-algorithm with shifts is used, though this routine
!   is not a fully general implementation. The algorithm is based on the
!   pseudo-code and description given in "Numerical Analysis" by
!   R. L. Burden and D. J. Faires (PWS-Kent 1989).
!
!- ---------------------------------------------------------------------
MODULE eigenvalue_tri_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'EIGENVALUE_TRI_MOD'
CONTAINS
SUBROUTINE eigenvalue_tri(n_matrix, n_in, d, e                          &
     , tol, n_max_iteration                                             &
     , nd_matrix)


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_warning
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
     nd_matrix
!      Size allocated for matrices treated together.


! Dummy arguments

  INTEGER, INTENT(IN) ::                                                &
     n_matrix                                                           &
!      Number of matrices treated together
   , n_in                                                               &
!      Order of input matrix
   , n_max_iteration
!      Maximum number of iterations

  REAL (RealK), INTENT(IN) ::                                           &
     tol
!      Tolerance for setting the subdiagonal elements to 0.

  REAL (RealK), INTENT(INOUT) ::                                        &
     d(nd_matrix, n_in)                                                 &
!      Main diagonal of the matrix: this will hold the eigenvalues
!      on output.
   , e(nd_matrix, n_in)
!      Subdiagonal of the matrix: E(1) is not used. E is reduced
!      to below the tolerance by the routine.


! Local Variables:

  INTEGER ::                                                            &
     n                                                                  &
!      Current working size of the problem
   , l                                                                  &
!      Loop variable
   , j                                                                  &
!      Loop variable
   , iteration
!      Current iteration
  REAL (RealK) ::                                                       &
     shift(n_matrix)                                                    &
!      Accumulated `shift'
   , d_shift(n_matrix)                                                  &
!      Increment in `shift'
   , b                                                                  &
!      Temporary variable used in solving quadratic
   , c                                                                  &
!      Temporary variable used in solving quadratic
   , discr                                                              &
!      Discriminant used in solving quadratic
   , kappa_1                                                            &
!      First root of quadratic
   , kappa_2
!      Second root of quadratic
  REAL (RealK) ::                                                       &
     abs_e                                                              &
!      Maximum absolute value of diagonal elements of the
!      current rows of the matrices
   , sinr(n_matrix)                                                     &
!      Sine of current rotation
   , cosr(n_matrix)                                                     &
!      Cosine of current rotations
   , cosr_temp                                                          &
!      Temporary cosine
   , sq                                                                 &
!      Temporary square root
   , sup_diag(n_matrix)                                                 &
!      Element of first superdiagonal of matrix on the J-1st row
   , sup_diag_old(n_matrix)
!      Previous value of SUP_DIAG

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  INTEGER                       :: ierr
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'EIGENVALUE_TRI'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The algorithm proceeds iteratively. The matrix supplied, A, is
! decomposed as A=QR where Q is orthogonal and R is upper
! triangular. A'=RQ is then formed and the process is repeated with
! A'. This leads to a sequence of matrices which converge to one
! with the eigenvalues down the diagonal.

! Initialization:
! Reduce the working size of the matrix if the off-diagonal
! elements are small enough.
  n=n_in
  abs_e=0.0e+00_RealK
  DO l=1, n_matrix
    abs_e=MAX(abs_e, ABS(e(l, n)))
  END DO
  DO WHILE ( (n >  1).AND.(abs_e <  tol) )
    n=n-1
    DO l=1, n_matrix
      abs_e=MAX(abs_e, ABS(e(l, n)))
    END DO
  END DO

  iteration=0
  DO l=1, n_matrix
    shift(l)=0.0e+00_RealK
  END DO


  DO WHILE ( (n >  1).AND.(iteration <  n_max_iteration) )


    iteration=iteration+1

!   Form an estimate of the first eigenvalue to be found by
!   finding the eigenvalues of the 2x2 matrix at the bottom
!   right-hand corner.
    DO l=1, n_matrix
      b=d(l, n-1)+d(l, n)
      c=d(l, n-1)*d(l, n)-e(l, n)*e(l, n)
      discr=SQRT(b*b-4.0e+00_RealK*c)
!     For reasons of conditioning we calculate the root of largest
!     magnitude and determine the other from the product of the
!     roots.
      kappa_1=0.5e+00_RealK*(b+SIGN(discr, b))
      kappa_2=c/kappa_1

!     Calculate the `shift' so as to accelerate convergence to the
!     last eigenvalue. A simple two-branch IF-test should be
!     amenable to vectorization if the vector CPU has a vector
!     mask register.
      IF ( ABS(kappa_1-d(l, n)) <                                       &
           ABS(kappa_2-d(l, n)) ) THEN
        d_shift(l)=kappa_1
      ELSE
        d_shift(l)=kappa_2
      END IF
      shift(l)=shift(l)+d_shift(l)
    END DO

!   Shift the diagonal elements.
    DO j=1, n
      DO l=1, n_matrix
        d(l, j)=d(l, j)-d_shift(l)
      END DO
    END DO


!   Form the QR-decompostion of the matrix by constructing
!   rotations to eliminate the sub-diagonal elements. COSR(J)
!   and SINR(J) are the cosine and sine of the rotations to
!   eliminate the element (J, J-1) of the input matrix: these
!   values specify the transpose of Q as we really construct
!   R=Qt.A by this procedure. The upper triangular matrix, R,
!   has two superdiagonals, but in practice only the first
!   is required. As the resulting matrix, RQ, will be a
!   symmetric tridaigonal matrix only its diagonal, D, and
!   the sub-diagonal, E, need be formed.

!   Inintialize:
    DO l=1, n_matrix
      sup_diag(l)=e(l, 2)
      cosr(l)=1.0e+00_RealK
      sinr(l)=0.0e+00_RealK
    END DO

    DO j=2, n

      DO l=1, n_matrix

!       This block of code is a little opaque as the variables
!       SINR and COSR are re-used to avoid the need to declare
!       them explicitly as vectors. We form the rotation to
!       elminate E(J) and also calculate E(J-1) of the new matrix
!       RQ using SINR(J-1) for the last time. The new cosine of
!       the rotation must be stored because we still need
!       the old one.
        sq=SQRT(d(l, j-1)*d(l, j-1)+e(l, j)*e(l, j))
        e(l, j-1)=sinr(l)*sq
        sinr(l)=e(l, j)/sq
        cosr_temp=d(l, j-1)/sq

!       Adjust the superdiagonal of the previous row of the matrix
!       as required by the elimination. The calculation of D(J-1)
!       actually belongs to the formation of RQ, but is done here
!       before we overwrite COSR.
        sup_diag_old(l)=sup_diag(l)
        sup_diag(l)=cosr_temp*sup_diag(l)+sinr(l)*d(l, j)
        d(l, j-1)=cosr(l)*d(l, j-1)+sinr(l)*sup_diag(l)
        cosr(l)=cosr_temp

!       Adjustments to the current row:
        d(l, j)=-sinr(l)*sup_diag_old(l)+cosr(l)*d(l, j)
        IF (j <  n) sup_diag(l)=cosr(l)*e(l, j+1)

      END DO

    END DO

    DO l=1, n_matrix
      e(l, n)=sinr(l)*d(l, n)
      d(l, n)=cosr(l)*d(l, n)
    END DO


!   Test for convergence and `shift' the converged eigenvalues.
!   back to their true values.
    abs_e=0.0e+00_RealK
    DO l=1, n_matrix
      abs_e=MAX(abs_e, ABS(e(l, n)))
    END DO
    DO WHILE ( (n >  1).AND.(abs_e <  tol) )
      DO l=1, n_matrix
        d(l, n)=d(l, n)+shift(l)
      END DO
      n=n-1
      DO l=1, n_matrix
        abs_e=MAX(abs_e, ABS(e(l, n)))
      END DO
    END DO


  END DO


! Check that convergence has occurred.
  IF (n >  1) THEN
    cmessage = '*** Warning: Convergence has not occurred while '       &
      //'calculating eigenvalues. The calculation continues.'
    ierr=i_warning
    CALL ereport(RoutineName, ierr, cmessage)
  ELSE
!   Shift the first eigenvalue back to its true value.
    DO l=1, n_matrix
      d(l, 1)=d(l, 1)+shift(l)
    END DO
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eigenvalue_tri
END MODULE eigenvalue_tri_mod

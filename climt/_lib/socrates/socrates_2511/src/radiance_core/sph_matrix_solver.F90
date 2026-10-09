! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve a set of stepped block equations.
!
! Method:
!   A set of linear equations of the stepped block form
!   is solved using Gaussian elimination with pivoting by rows
!   over a restricted range. In this application it should
!   not be necessary to consider all potential equations
!   when choosing pivots and this helps to reduce the
!   band-width of the system.
!
!- ---------------------------------------------------------------------
MODULE sph_matrix_solver_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SPH_MATRIX_SOLVER_MOD'
CONTAINS
SUBROUTINE sph_matrix_solver(n_matrix, n_step, n_block                  &
    , a, b                                                              &
    , x                                                                 &
    , nd_matrix, nd_equation, nd_diagonal                               &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_matrix                                                         &
!       Size alloacted for matrices
    , nd_equation                                                       &
!       Size allocated for equations
    , nd_diagonal
!       Size allocated for diagonals

! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      n_matrix                                                          &
!       Number of matrices
    , n_step                                                            &
!       Number of steps in the matrix
    , n_block
!       Size of each block
  REAL (RealK), INTENT(INOUT) ::                                        &
      a(nd_matrix, nd_equation, nd_diagonal)                            &
!       Matrices of coefficients
    , b(nd_matrix, nd_equation)
!       Righthand sides
  REAL (RealK), INTENT(OUT) ::                                          &
       x(nd_matrix, nd_equation)
!       Solution vector


! Local variables
  INTEGER                                                               &
      i_step                                                            &
!       Counter for steps in the matrix
    , i_phase                                                           &
!       Counter for phase of elimination or back-substitution
    , ie                                                                &
!       Number of equation
    , ic                                                                &
!       Column of the compressed matrix
    , right                                                             &
!       Rightmost column of the compressed matrix
    , row_first                                                         &
!       First row of the matrix considered in the second phase
!       of elimination
    , row_last
!       Last row of the matrix considered in the second phase
!       of elimination
  INTEGER                                                               &
      i_pivot(nd_matrix)                                                &
!       Index of pivot
    , offset_pivot(nd_matrix)
!       Offset of the pivoting element relative to the current row
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      pivot(nd_matrix)                                                  &
!       Absolute value of the pivoting element
    , rho(nd_matrix)                                                    &
!       Scaling applied to the pivotal row in elimination
    , aabs                                                              &
!       Absolute value of the current potential pivot
    , tmp
!       Temporary variable used in swapping rows

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SPH_MATRIX_SOLVER'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Eliminative phase:
  i_step=1
  DO WHILE (i_step <= n_step)


!   The elminination falls naturally into two phases given
!   the structure of the matrix.
    DO i_phase=1, 2

      DO i=1, n_block

!       Set the number of the equation for elimination
!       and its column in the reduced matrix.
        ie=n_block*(2*i_step+i_phase-3)+i
        ic=n_block*(3-i_phase)+i
        IF (i_step <  n_step) THEN
          right=n_block*(8-2*i_phase)
        ELSE
          right=n_block*(6-2*i_phase)
        END IF

!       Choose the row for pivoting.
        DO l=1, n_matrix
          pivot(l)=ABS(a(l, ie, ic))
          i_pivot(l)=ie
          offset_pivot(l)=0
        END DO
!       In the next line (and a little later), we limit the
!       range to the size of the matrix because in the last
!       layer there is no further boundary to consider.
        DO j=ie+1, MIN(ie+n_block*i_phase-i, 2*n_step*n_block)
          DO l=1, n_matrix
            aabs=ABS(a(l, j, ic))
            IF (aabs >  pivot(l)) THEN
              pivot(l)=aabs
              i_pivot(l)=j
            END IF
          END DO
        END DO
!       In the first phase we also need to consider rows
!       which will have a different offset as they involve
!       conditions on the next boundary below.
        IF (i_phase == 1) THEN
          DO j=n_block*(2*i_step-1)+1                                   &
            , n_block*(MIN(2*i_step+1, 2*n_step))
            DO l=1, n_matrix
              aabs=ABS(a(l, j, ic-2*n_block))
              IF (aabs >  pivot(l)) THEN
                pivot(l)=aabs
                i_pivot(l)=j
                offset_pivot(l)=2*n_block
              END IF
            END DO
          END DO
        END IF

!       Swap the rows regardlessly to allow vectorization.
        DO j=ic, right
          DO l=1, n_matrix
            tmp=a(l, ie, j)
            a(l, ie, j)=a(l, i_pivot(l), j-offset_pivot(l))
            a(l, i_pivot(l), j-offset_pivot(l))=tmp
          END DO
        END DO
        DO l=1, n_matrix
          tmp=b(l, ie)
          b(l, ie)=b(l, i_pivot(l))
          b(l, i_pivot(l))=tmp
        END DO

!       Now eliminate. In both phases we have equations dealing
!       with the same boundary as the pivoting equation, but in
!       the first phase there is also an additional elimination
!       to be performed referring to the next boundary below.
        row_first=ie+1
        IF (i_step <  n_step) THEN
          row_last=ie+n_block*i_phase-i
        ELSE
          row_last=ie+n_block-i
        END IF
        DO j=row_first, row_last
          DO l=1, n_matrix
            rho(l)=a(l, j, ic)/a(l, ie, ic)
            b(l, j)=b(l, j)-rho(l)*b(l, ie)
          END DO
          DO k=ic+1, right
            DO l=1, n_matrix
              a(l, j, k)=a(l, j, k)-rho(l)*a(l, ie, k)
            END DO
          END DO
        END DO
!       This is the extra elimination required during the first
!       phase.
        IF (i_phase == 1) THEN
          row_first=n_block*(2*i_step-1)+1
          IF (i_step <  n_step) THEN
            row_last=n_block*(2*i_step+1)
          ELSE
            row_last=n_block*2*i_step
          END IF
          DO j=row_first, row_last
            DO l=1, n_matrix
              rho(l)=a(l, j, ic-2*n_block)/a(l, ie, ic)
              b(l, j)=b(l, j)-rho(l)*b(l, ie)
            END DO
            DO k=ic+1, right
              DO l=1, n_matrix
                a(l, j, k-2*n_block)                                    &
                  =a(l, j, k-2*n_block)-rho(l)*a(l, ie, k)
              END DO
            END DO
          END DO
        END IF

      END DO
    END DO

    i_step=i_step+1

  END DO



! Back-subsititution:
  i_step=n_step
  DO WHILE (i_step >= 1)


    DO i_phase=2, 1, -1

      IF (i_step <  n_step) THEN
        right=n_block*(8-2*i_phase)
      ELSE
        right=n_block*(6-2*i_phase)
      END IF

      DO i=n_block, 1, -1

        ie=n_block*(2*i_step+i_phase-3)+i
        ic=(3-i_phase)*n_block+i

        DO l=1, n_matrix
          x(l, ie)=b(l, ie)
        END DO
        DO j=1, right-ic
          DO l=1, n_matrix
            x(l, ie)=x(l, ie)-a(l, ie, j+ic)*x(l, ie+j)
          END DO
        END DO
        DO l=1, n_matrix
          x(l, ie)=x(l, ie)/a(l, ie, ic)
        END DO

      END DO
    END DO

    i_step=i_step-1


  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sph_matrix_solver
END MODULE sph_matrix_solver_mod

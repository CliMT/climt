! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set up and solve the eigensystem.
!
! Purpose:
!   For a given value of the azimuthal quantum number, MS, this
!   routine returns the positive eigenvalues imposed by the trunctaion
!   in one layer and the corresponsing eigenvectors.
!
! Method:
!   The sub-diagonal of the full matrix is calculated and then reduced
!   to the diagonal and subdiagonal of the reduced matrix. The
!   eigenvalues are then found by calling the QR-algorithm and the
!   eigenvectors are obtained from a recurrence relation.
!
!- ---------------------------------------------------------------------
MODULE eig_sys_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'EIG_SYS_MOD'
CONTAINS
SUBROUTINE eig_sys(n_profile, ls_trunc, ms, n_red_eigensystem           &
    , cg_coeff, sqs                                                     &
    , mu, eig_vec                                                       &
    , nd_profile, nd_red_eigensystem, nd_max_order                      &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE eigenvalue_tri_mod, ONLY: eigenvalue_tri

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_red_eigensystem                                                &
!       Size allocated for the reduced eigensystem
    , nd_max_order
!       Size allocated for the order of the calculation


! Dummy variables
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric profiles
    , ls_trunc                                                          &
!       Order of L for truncation
    , ms                                                                &
!       Azimuthal quantum number
    , n_red_eigensystem
!       Size of the reduced eigenproblem
  REAL (RealK), INTENT(IN) ::                                           &
      cg_coeff(ls_trunc+1-ms)                                           &
!       Clebsch-Gordan coefficients
    , sqs(nd_profile, 0: nd_max_order)
!       Square roots of S-coefficients
  REAL (RealK), INTENT(OUT) ::                                          &
      mu(nd_profile, nd_red_eigensystem)                                &
!       Eigenvalues
    , eig_vec(nd_profile, 2*nd_red_eigensystem, nd_red_eigensystem)
!       Eigenvectors


! Local variables
  INTEGER                                                               &
      j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , ls                                                                &
!       Order of spherical harmonic
    , n_max_qr_iteration
!       Maximum number of QR iterations
  REAL (RealK) ::                                                       &
      tol                                                               &
!       Tolerance for assessing convergence of the QR-algorithm
    , ec(nd_profile, 2*n_red_eigensystem)                               &
!       Sub-diagonal of full matrix
    , e(nd_profile, n_red_eigensystem)                                  &
!       Sub-diagonal of reduced matrix
    , normalization(nd_profile)                                         &
!       Normalization factor for the eigenvector
    , c_ratio(nd_profile)                                               &
!       Common ratio in the geometric progression used to rescale
!       the recurrence relation to avoid overflows
    , rescale(nd_profile)
!       Multiplier to convert the terms of the scaled recurrence
!       to the the final quantities

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='EIG_SYS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerance for convergence of the algorithm from the
! precision of the machine.
  tol=1.0e+02*EPSILON(tol)

! Calculate the reduced matrix to yield the eigenvalues. EC_...
! represent elements of the sub-diagonal of the full matrix:
! D and E are the diagonal and sub-diagonal of the reduced matrix.


! Calculate the sub-diagonal of the full matrix.
  DO j=2, ls_trunc+1-ms
    ls=ms-1+j
    DO l=1, n_profile
      ec(l, j)=cg_coeff(j-1)/(sqs(l, ls)*sqs(l, ls-1))
    END DO
  END DO


! Retain odd rows and columns of the square of the preceeding
! matrix. The diagonal terms are stored in MU as this will be
! reduced to the eigenvalues later.
  DO l=1, n_profile
    mu(l, 1)=ec(l, 2)**2
  END DO
  DO j=2, n_red_eigensystem
    DO l=1, n_profile
      mu(l, j)=ec(l, 2*j-1)**2+ec(l, 2*j)**2
      e(l, j)=ec(l, 2*j-2)*ec(l, 2*j-1)
    END DO
  END DO

! Determine the eigenvalues of the reduced matrix, which
! are the squares of the (positive) eigenvalues of the
! full matrix. If the eigensystem is of size 1 no calculation
! is required.
  IF (n_red_eigensystem >  1) THEN
!   The number of iterations required for convergence increases
!   as the order of truncation rises. A small allowance is made
!   for extra iterations.
    n_max_qr_iteration=ls_trunc+25
    CALL eigenvalue_tri(n_profile, n_red_eigensystem                    &
      , mu, e, tol, n_max_qr_iteration                                  &
      , nd_profile                                                      &
      )
  END IF
  DO k=1, n_red_eigensystem
    DO l=1, n_profile
      mu(l, k)=SQRT(mu(l, k))
      IF (mu(l, k) >  1.0e+00_RealK) THEN
        c_ratio(l)=5.0e-01_RealK/mu(l, k)
      ELSE
        c_ratio(l)=1.0e+00_RealK
      END IF
    END DO

!   Use the recurrence relation to find the eigenvectors of the
!   full matrix. For large values of MU there will be an
!   eigenvector like MU^J and one like MU^-J. The latter (minimal)
!   solution is required, but for |MU|>1 the recurrence is
!   unstable, so the growing solution will will swamp the required
!   solution. Conversely, with downward recurrence, the desired
!   solution grows and will dominate in the recurrence. When
!   |MU|<1 the recurrence is stable in either direction so downward
!   recurrence is used consistently. On further complication must
!   be taken into account: if MU is very large (as can happen with
!   almost conservative scattering) the elements of the eigenvector
!   may be of so large a range of magnitudes that the recurrence
!   overflows. A scaling factor, c, is therefore introduced so that
!   the j'th element of the eigenvector, e_j=c^j.e_j'. s may not
!   be less than 1 for small eigenvalues or the same problem will
!   be introduced with them; the vector e' has elements of order 1.

    j=2*n_red_eigensystem
    DO l=1, n_profile
      eig_vec(l, j, k)=1.0e+00_RealK
    END DO
    j=j-1
    DO l=1, n_profile
      eig_vec(l, j, k)=c_ratio(l)*mu(l, k)/ec(l, j+1)
    END DO
    DO WHILE(j >  1)
      j=j-1
      DO l=1, n_profile
        eig_vec(l, j, k)                                                &
          =(mu(l, k)*eig_vec(l, j+1, k)                                 &
          -c_ratio(l)*ec(l, j+2)*eig_vec(l, j+2, k))                    &
          *c_ratio(l)/ec(l, j+1)
      END DO
    END DO

!   Remove the scaling factor, renormalize the eigenvector
!   and rescale by the s-coefficients for later efficiency.
    DO l=1, n_profile
      rescale(l)=c_ratio(l)
      eig_vec(l, 1, k)=eig_vec(l, 1, k)*rescale(l)
      normalization(l)=eig_vec(l, 1, k)*eig_vec(l, 1, k)
    END DO
    DO j=2, 2*n_red_eigensystem
      DO l=1, n_profile
        rescale(l)=rescale(l)*c_ratio(l)
        eig_vec(l, j, k)=eig_vec(l, j, k)*rescale(l)
        normalization(l)=normalization(l)                               &
          +eig_vec(l, j, k)*eig_vec(l, j, k)
      END DO
    END DO
    DO l=1, n_profile
      normalization(l)=SQRT(1.0e+00_RealK/normalization(l))
    END DO
    DO j=1, 2*n_red_eigensystem
      DO l=1, n_profile
        eig_vec(l, j, k)=eig_vec(l, j, k)*normalization(l)              &
          /sqs(l, j+ms-1)
      END DO
    END DO

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eig_sys
END MODULE eig_sys_mod

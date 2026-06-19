! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to perfrom back-substitution in SVD.
!
! Method:
!   Straightforward back-substitution.
!
!- ---------------------------------------------------------------------
      SUBROUTINE back_substitute(u, w, v
     &  , n_row, n_column
     &  , row_dimension, column_dimension
     &  , b, x
     &  , tmp
     &  )
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
!     Dummy variables.
      INTEGER, Intent(IN) ::
     &    n_row
!           Number of rows
     &  , n_column
!           Number of columns
     &  , row_dimension
!           First dimension
     &  , column_dimension
!           Second dimension
      REAL  (RealK), Intent(IN) ::
     &    u(row_dimension, column_dimension)
!           U-matrix
     &  , w(column_dimension)
!           Diagonal matrix
     &  , v(column_dimension, column_dimension)
!           V-matrix
     &  , b(row_dimension)
!           Right hand side
!
      REAL  (RealK), Intent(OUT) ::
     &    x(column_dimension)
!           Solution
!
!     Local variables.
      INTEGER
     &    i
!           Loop variable
     &  , j
!           Loop variable
     &  , jj
!           Loop variable
      REAL  (RealK) ::
     &    s
!           Summing variable
     &  , tmp(column_dimension)
!           Temporary store
!
!     Variables related to the treatment of ill-conditioning
      REAL  (RealK) ::
     &    tol_svd
!           The tolerance for SVD decomposition
!
!
!
!     Set the tolerances used in avoiding ill-conditioning.
      tol_svd=1.0e+03_RealK*epsilon(tol_svd)
!
!
!     Calculate UT.B
      DO j=1, n_column
        s=0.0_RealK
        IF (abs(w(j)) > tol_svd) THEN
          DO i=1, n_row
            s=s+u(i, j)*b(i)
          ENDDO
            s=s/w(j)
        ENDIF
        tmp(j)=s
      ENDDO
      DO j=1, n_column
        s=0.0_RealK
        DO jj=1, n_column
          s=s+v(j, jj)*tmp(jj)
        ENDDO
        x(j)=s
      ENDDO
!
!
!
      RETURN
      END

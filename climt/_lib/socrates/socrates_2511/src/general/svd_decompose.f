! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to perform SVD decomposition of a matrix.
!
!- ---------------------------------------------------------------------
      SUBROUTINE svd_decompose(ierr
     &  , a, n_row, n_column
     &  , row_dimension, column_dimension
     &  , w, v, rv1
     &  )

      USE realtype_rd
      USE rad_pcf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n_row
!           Number of rows
     &  , n_column
!           Number of columns
     &  , row_dimension
!           First dimension
     &  , column_dimension
!           Second dimension
      REAL  (RealK), Intent(INOUT) ::
     &    a(row_dimension, column_dimension)
!           SVD matrix: overwritten
      REAL  (RealK), Intent(OUT) ::
     &    w(row_dimension)
!           Diagonal terms
     &  , v(column_dimension, column_dimension)
!           Second matrix
     &  , rv1(column_dimension)
!           Working array
!
!     Local arguments.
      INTEGER
     &    i
!           Loop variable
     &  , j
!             "
     &  , k
!             "
     &  , l
!             "
     &  , iteration
!           Iteration number
     &  , nm
!           Temporary scalar
     &  , jj
!             "
      REAL  (RealK) ::
     &    scale
!           These variables are just
     &  , anorm
!           Part of the algorithm
     &  , g
!           And cnnot meaningfully
     &  , f
!           Be describeed
     &  , h
     &  , y
     &  , z
     &  , s
     &  , c
     &  , x

      INTEGER, PARAMETER :: max_svd_iteration=100
!           Maximum number of svd iterations

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
!     Check that there are not too many points in the fit.
      IF (n_column > n_row) THEN
        WRITE(iu_err, '(a)')
     &    'svd attempted with too many unknowns.'
        ierr=i_abort_calculation
        RETURN
      ENDIF
!
!     Householder reduction to bidiagonal form.
      g=0.0_RealK
      scale=0.0_RealK
      anorm=0.0_RealK
      DO i=1, n_column
        l=i+1
        rv1(i)=scale*g
        g=0.0_RealK
        s=0.0_RealK
        scale=0.0_RealK
        IF (i <= n_row) THEN
          DO k=i, n_row
            scale=scale+abs(a(k,i))
          ENDDO
          IF (scale > tol_svd) THEN
            DO k=i, n_row
              a(k, i)=a(k, i)/scale
              s=s+a(k, i)*a(k, i)
            ENDDO
            f=a(i, i)
            g=-sign(sqrt(s), f)
            h=f*g-s
            a(i, i)=f-g
            IF (i /= n_column) THEN
              DO j=l, n_column
                s=0.0_RealK
                DO k=i, n_row
                  s=s+a(k, i)*a(k, j)
                ENDDO
                f=s/h
                DO k=i, n_row
                  a(k, j)=a(k, j)+f*a(k, i)
                ENDDO
              ENDDO
            ENDIF
            DO k=i, n_row
              a(k, i)=scale*a(k, i)
            ENDDO
          ENDIF
        ENDIF
!
        w(i)=scale*g
        g=0.0_RealK
        s=0.0_RealK
        scale=0.0_RealK
        IF ( (i <= n_row).AND.(i /= n_column) ) THEN
          DO k=l, n_column
            scale=scale+abs(a(i, k))
          ENDDO
          IF (scale > tol_svd) THEN
            DO k=l, n_column
              a(i, k)=a(i, k)/scale
              s=s+a(i, k)*a(i, k)
            ENDDO
            f=a(i, l)
            g=-sign(sqrt(s), f)
            h=f*g-s
            a(i, l)=f-g
            DO k=l, n_column
              rv1(k)=a(i, k)/h
            ENDDO
            IF (i /= n_row) THEN
              DO j=l, n_row
                s=0.0_RealK
                DO k=l, n_column
                  s=s+a(j, k)*a(i, k)
                ENDDO
                DO k=l, n_column
                  a(j, k)=a(j, k)+s*rv1(k)
                ENDDO
              ENDDO
            ENDIF
            DO k=l, n_column
              a(i, k)=scale*a(i, k)
            ENDDO
          ENDIF
        ENDIF
        anorm=max(anorm, (abs(w(i))+abs(rv1(i))))
      ENDDO
!
!     Accumulate the right-hand sides of the transformations.
      DO i=n_column, 1, -1
        IF (i < n_column) THEN
          IF (abs(g) > tol_svd) THEN
            DO j=l, n_column
              v(j, i)=(a(i, j)/a(i, l))/g
            ENDDO
            DO j=l, n_column
              s=0.0_RealK
              DO k=l, n_column
                s=s+a(i, k)*v(k, j)
              ENDDO
              DO k=l, n_column
                v(k, j)=v(k, j)+s*v(k, i)
              ENDDO
            ENDDO
          ENDIF
          DO j=l, n_column
            v(i, j)=0.0_RealK
            v(j, i)=0.0_RealK
          ENDDO
        ENDIF
        v(i, i)=1.0_RealK
        g=rv1(i)
        l=i
      ENDDO
!
!     Accumulate the left-hand sides of the transformations.
      DO i=n_column, 1, -1
        l=i+1
        g=w(i)
        IF (i < n_column) THEN
          DO j=l, n_column
            a(i, j)=0.0_RealK
          ENDDO
        ENDIF
        IF (abs(g) > tol_svd) THEN
          g=1.0_RealK/g
          IF (i /= n_column) THEN
            DO j=l, n_column
              s=0.0_RealK
              DO k=l, n_row
                s=s+a(k, i)*a(k, j)
              ENDDO
              f=(s/a(i, i))*g
              DO k=i, n_row
                a(k, j)=a(k, j)+f*a(k, i)
              ENDDO
            ENDDO
          ENDIF
          DO j=i, n_row
            a(j, i)=a(j, i)*g
          ENDDO
        ELSE
          DO j=i, n_row
            a(j, i)=0.0_RealK
          ENDDO
        ENDIF
        a(i, i)=a(i, i)+1.0_RealK
      ENDDO
!
!     Diagonalize the bidiagonal form.
      DO k=n_column, 1, -1
        DO iteration=1, max_svd_iteration
          DO l=k, 1, -1
            nm=l-1
            IF (abs(rv1(l)) < tol_svd) goto 2
            IF (abs(w(nm)) < tol_svd) goto 1
          ENDDO
!
1         c=0.0_RealK
          s=1.0_RealK
          DO i=l, k
            f=s*rv1(i)
            IF (abs(f) < tol_svd) THEN
              g=w(i)
              h=sqrt(f*f+g*g)
              w(i)=h
              h=1.0_RealK/h
              c=g*h
              s=-f*h
              DO j=1, n_row
                y=a(j, nm)
                z=a(j, i)
                a(j, nm)=(y*c)+(z*s)
                a(j, i)=-(y*s)+(z*c)
              ENDDO
            ENDIF
          ENDDO
!
2         z=w(k)
          IF (l == k) THEN
            IF (z < 0.0_RealK) THEN
              w(k)=-z
              DO j=1, n_column
                v(j, k)=-v(j, k)
              ENDDO
            ENDIF
            goto 3
          ENDIF
!
          IF (iteration == max_svd_iteration) THEN
            WRITE(iu_err, '(/a)')
     &        'convergence of svd not achieved: too many iterations.'
          ENDIF
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_RealK*h*y)
          g=sqrt(f*f+1.0_RealK)
          f=((x-z)*(x+z)+h*((y/(f+sign(g, f)))-h))/x
!
!         Make the qr-transformation.
          c=1.0_RealK
          s=1.0_RealK
          DO j=l, nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=sqrt(f*f+h*h)
            rv1(j)=z
            c=f/z
            s=h/z
            f=(x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            DO jj=1, n_column
              x=v(jj, j)
              z=v(jj, i)
              v(jj, j)=(x*c)+(z*s)
              v(jj, i)=-(x*s)+(z*c)
            ENDDO
            z=sqrt(f*f+h*h)
            w(j)=z
            IF (z > tol_svd) THEN
              z=1.0_RealK/z
              c=f*z
              s=h*z
            ENDIF
            f=(c*g)+(s*y)
            x=-(s*g)+(c*y)
            DO jj=1, n_row
              y=a(jj, j)
              z=a(jj, i)
              a(jj, j)=(y*c)+(z*s)
              a(jj, i)=-(y*s)+(z*c)
            ENDDO
          ENDDO
          rv1(l)=0.0_RealK
          rv1(k)=f
          w(k)=x
        ENDDO
3     continue
      ENDDO
!
!
!
      RETURN
      END

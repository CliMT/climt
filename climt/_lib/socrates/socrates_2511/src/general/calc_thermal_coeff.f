! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the matrices for a fit to the Plack function.
!
! Method:
!   A polynomial fit in the temperature to the Planck function 
!   integrated across a spectral band is required. This is obtained
!   from a least-squares fit. Minimization of the error in the fit
!   leads to a set of linear equations. This routine determines the
!   matrix and source term in these equations as A and B.
!
!   This routine does basically the same as 'calc_thermal_coeff' below.
!   However, because now we have to integrate the product of the Planckian
!   and the filter function the least square forcing terms have to be 
!   evaluated numerically. 
!   The least squares forcing terms are evaluated by integrating the
!   product of the Planckian and the filter function over wavelength and
!   temperature. Integration is done by using Simpson's rule.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_THERMAL_COEFF_NUM(IERR
     &  , N_DEG_FIT, N_PTS
     &  , LAMBDA_MIN, LAMBDA_MAX
     &  , THETA_LOW, THETA_HIGH
     &  , C1, C2
     &  , WEIGHTS
     &  , INNER
     &  , B)

      USE realtype_rd
      USE rad_pcf
      USE dimensions_pp_ucf
      USE dimensions_spec_ucf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::
     &    N_DEG_FIT
!           Degree of fit to function
     &  , N_PTS
!           Number of points in filter function     
      REAL  (RealK), INTENT(IN) ::      
     &    LAMBDA_MIN
!           Short wavelength of band
     &  , LAMBDA_MAX
!           Long wavelength of band     
     &  , THETA_LOW
!           Low scaled temperature
     &  , THETA_HIGH
!           High scaled temperature    
     &  , C1
!           Numerical constant
     &  , C2
!           Numerical constant
     &  , WEIGHTS(0:N_PTS-1)
!           Weights contained in the filter function

      REAL  (RealK), INTENT(OUT) ::    
     &    INNER(0: NPD_THERMAL_COEFF-1, 0: NPD_THERMAL_COEFF-1) 
!           Matrix of inner product
     &  , B(0: NPD_THERMAL_COEFF-1)
!           Least squares forcing term      
     
! Local Variables
  
      INTEGER
     &    I, J, K
!           Looping variables
     &  , N_POWER
!           Power of temperature

      REAL (Realk) ::
     &    THETA(0:NPD_FIT_TEMP)
     &  , INTEG
!           Planckian
     &  , COEFF
     &  , WEIGHT_PLANCKIAN(0:NPD_FIT_TEMP,0:NPD_THERMAL_COEFF-1) 
!           Product of Planckian and weights                      
     &  , EXPON
!          Exponential in Planckian
     &  , LAMBDA
!           Wavelength
     &  , D_LAMBDA
!            Wavelength increment
     &  ,  D_THETA
!           Temperature increment 
     &  ,  DUMMY(0:N_PTS-1)
!           Dummy variable     

! Calculate the terms of the matrix A first of all. These terms
! represent the inner products of the basis functions over temperature.

      DO I=0, N_DEG_FIT
        DO J=0, I
          N_POWER=I+J+1
          INNER(I, J)=(THETA_HIGH**N_POWER-THETA_LOW**N_POWER)
     &      /REAL(N_POWER, RealK)
          INNER(J, I)=INNER(I, J)
        ENDDO
      ENDDO


! The terms to the right are more complicated, involving the 
! integral of the planckian times weights over wavelength 
! and an inner product over temperature      
 
      D_LAMBDA=(LAMBDA_MAX-LAMBDA_MIN)/REAL(N_PTS-1, RealK)
      D_THETA=(THETA_HIGH-THETA_LOW)/REAL(NPD_FIT_TEMP-1, RealK)
      
!      print*,'D_LAMBDA=',D_LAMBDA
      
     
      DO K=0, N_DEG_FIT
         DO J=0,NPD_FIT_TEMP
       
            THETA(J)=THETA_LOW+J*D_THETA
      
            INTEG=0.0
        
            DO I=0,N_PTS-1
            
               LAMBDA=LAMBDA_MIN+I*D_LAMBDA
               EXPON =EXP(C2/(LAMBDA*THETA(J))) 
               DUMMY(I)=WEIGHTS(I)*C1/(LAMBDA**5*(EXPON-1.0))
            ENDDO

            CALL SIMPSONS_RULE(N_PTS-1,D_LAMBDA,DUMMY,INTEG)
            
!           print*,j,theta(j),integ
            
            WEIGHT_PLANCKIAN(J,K)=INTEG*THETA(J)**K
            
         ENDDO
         
         print*,''   
      ENDDO 
 
! Integrate over temperature using Simpson's Rule

      DO K=0, N_DEG_FIT
        B(K)=0.0
         DO J=0,NPD_FIT_TEMP
   
            CALL SIMPSONS_RULE(NPD_FIT_TEMP,D_THETA,
     &                         WEIGHT_PLANCKIAN(0,K),B(K))

         ENDDO
      ENDDO
       
            
      RETURN
      END
!-------------------------------------------------------------------------
!
!+ Subroutine to calculate the matrices for a fit to the Plack function.
!
! Method:
!   A polynomial fit in the temperature to the Planck function 
!   integrated across a spectral band is required. This is obtained
!   from a least-squares fit. Minimization of the error in the fit
!   leads to a set of linear equations. This routine determines the
!   matrix and source term in these equations as A and B.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_THERMAL_COEFF(IERR
     &  , N_DEG_FIT
     &  , WAVELENGTH_SHORT, WAVELENGTH_LONG
     &  , THETA_LOW, THETA_HIGH
     &  , ALPHA, Q, A, B)

      USE realtype_rd
      USE rad_pcf
      USE dimensions_pp_ucf
      USE dimensions_spec_ucf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::
     &    N_DEG_FIT
!           Degree of fit to function
      REAL  (RealK), INTENT(IN) ::
     &    WAVELENGTH_SHORT
!           Short wavelength of band
     &  , WAVELENGTH_LONG
!           Long wavelength of band
     &  , THETA_LOW
!           Low scaled temperature
     &  , THETA_HIGH
!           High scaled temperature
     &  , ALPHA
!           Numerical constant
     &  , Q
!           Numerical constant
      REAL  (RealK), INTENT(OUT) ::
     &    A(0: NPD_THERMAL_COEFF-1, 0: NPD_THERMAL_COEFF-1)
!           Matrix of inner products
     &  , B(0: NPD_THERMAL_COEFF-1)
!           Rhs for least squares eqn.
!
!     Local arguments.
      INTEGER
     &    N_POWER
!           Power of temperature
     &  , N_POINT
!           Number of integration points
     &  , I
!           Loop variable
     &  , J
!           Loop variable
      REAL  (RealK) ::
     &    U_1
!           Limit of integral
     &  , U_2
!           Limit of integral
     &  , U_3
!           Limit of integral
     &  , U_4
!           Limit of integral
     &  , INTEGRAL_OLD
!           Old estimate of integral
     &  , INTEGRAL_NEW
!           Old estimate of integral
     &  , U(NPD_THERMAL_ABSCISSA)
!           Abscissae for u-integration
     &  , INTEGRAND(NPD_THERMAL_ABSCISSA)
!           Integrands for this integration

      REAL (RealK), PARAMETER :: tol_thermal=1.0e-06_RealK
!           Tolerance on integration

!     Functions called:
      REAL  (RealK) ::
     &    TRAPEZOID
!           Function for trapezoidal integration
      EXTERNAL 
     &     TRAPEZOID
!
!
!
!
!     Calculate the terms of the matrix A first of all. These terms
!     represent the inner products of the basis functions 
!     over temperature.
      DO I=0, N_DEG_FIT
        DO J=0, I
          N_POWER=I+J+1
          A(I, J)=(THETA_HIGH**N_POWER-THETA_LOW**N_POWER)
     &      /REAL(N_POWER, RealK)
          A(J, I)=A(I, J)
        ENDDO
      ENDDO
!
!     The terms to the right are more complicated, involving the 
!     integral of the planckian function over wavelength 
!     and an inner product over temperature. The order of 
!     integration has been inverted and the temperature integration 
!     performed analytically, leaving an integration over the 
!     transformed variable
!              U=(H_PLANCK*NU)/(K_BLOTZMANN*T_REF_THERMAL*THETA)
!     to be performed numerically. There are two cases.
!
      IF (WAVELENGTH_LONG*THETA_LOW.GE.WAVELENGTH_SHORT*THETA_HIGH) THEN
        U_1=1.0E+00_RealK/(ALPHA*THETA_HIGH*WAVELENGTH_LONG)
        U_2=1.0E+00_RealK/(ALPHA*THETA_LOW*WAVELENGTH_LONG)
        U_3=1.0E+00_RealK/(ALPHA*THETA_HIGH*WAVELENGTH_SHORT)
        U_4=1.0E+00_RealK/(ALPHA*THETA_LOW*WAVELENGTH_SHORT)
!
        DO I=0, N_DEG_FIT
!
          B(I)=0.0E+00_RealK
!
!         First integral:
          U(1)=U_1
          U(2)=0.5E+00_RealK*(U_1+U_2)
          U(3)=U_2
          INTEGRAL_OLD=0.0E+00_RealK
          N_POINT=3
1         DO J=1, N_POINT
            INTEGRAND(J)=FNC_11(U(J), I)
          ENDDO
          INTEGRAL_NEW=TRAPEZOID(N_POINT, U, INTEGRAND)
          IF (ABS(INTEGRAL_NEW-INTEGRAL_OLD).GT.ABS(INTEGRAL_NEW)
     &      *TOL_THERMAL) THEN
            INTEGRAL_OLD=INTEGRAL_NEW
            IF (2*N_POINT.GT.NPD_THERMAL_ABSCISSA) THEN
              WRITE(IU_ERR, '(/A, /A)')
     &          '** Error: Too many points are required to '
     &          //'evaluate the source function. '
     &          , 'increase NPD_THERMAL_ABSCISSA and re-compile.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            DO J=N_POINT, 2, -1
              U(2*J-1)=U(J)
              U(2*J-2)=0.5E+00_RealK*(U(J)+U(J-1))
            ENDDO
            N_POINT=2*N_POINT-1
            GOTO 1
          ENDIF
!         Use Richardson extrapolation to remove 
!         the leading error term:
          INTEGRAL_NEW=(4.0E+00_RealK*INTEGRAL_NEW-INTEGRAL_OLD)
     &      /3.0E+00_RealK
          B(I)=B(I)+INTEGRAL_NEW
!
!         Second integral:
          U(1)=U_2
          U(2)=0.5E+00_RealK*(U_2+U_3)
          U(3)=U_3
          INTEGRAL_OLD=0.0E+00_RealK
          N_POINT=3
2         DO J=1, N_POINT
            INTEGRAND(J)=FNC_12(U(J), I)
          ENDDO
          INTEGRAL_NEW=TRAPEZOID(N_POINT, U, INTEGRAND)
          IF (ABS(INTEGRAL_NEW-INTEGRAL_OLD).GT.TOL_THERMAL
     &       *ABS(INTEGRAL_NEW)) THEN
            INTEGRAL_OLD=INTEGRAL_NEW
            IF (2*N_POINT.GT.NPD_THERMAL_ABSCISSA) THEN
              WRITE(IU_ERR, '(/A, /A)')
     &          '** Error: Too many points are required to '
     &          //'evaluate the source function. '
     &          , 'increase NPD_THERMAL_ABSCISSA and re-compile.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            DO J=N_POINT, 2, -1
              U(2*J-1)=U(J)
              U(2*J-2)=0.5E+00_RealK*(U(J)+U(J-1))
            ENDDO
            N_POINT=2*N_POINT-1
            GOTO 2
          ENDIF
!         Use Richardson extrapolation to remove 
!         the leading error term:
          INTEGRAL_NEW=(4.0E+00_RealK*INTEGRAL_NEW-INTEGRAL_OLD)
     &      /3.0E+00_RealK
          B(I)=B(I)+INTEGRAL_NEW
!
!         Third integral:
          U(1)=U_3
          U(2)=0.5E+00_RealK*(U_3+U_4)
          U(3)=U_4
          INTEGRAL_OLD=0.0E+00_RealK
          N_POINT=3
3         DO J=1, N_POINT
            INTEGRAND(J)=FNC_13(U(J), I)
          ENDDO
          INTEGRAL_NEW=TRAPEZOID(N_POINT, U, INTEGRAND)
          IF (ABS(INTEGRAL_NEW-INTEGRAL_OLD).GT.TOL_THERMAL
     &       *ABS(INTEGRAL_NEW)) THEN
            INTEGRAL_OLD=INTEGRAL_NEW
            IF (2*N_POINT.GT.NPD_THERMAL_ABSCISSA) THEN
              WRITE(IU_ERR, '(/A, /A)')
     &          '** Error: Too many points are required to '
     &          //'evaluate the source function. '
     &          , 'increase NPD_THERMAL_ABSCISSA and re-compile.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            DO J=N_POINT, 2, -1
              U(2*J-1)=U(J)
              U(2*J-2)=0.5E+00_RealK*(U(J)+U(J-1))
            ENDDO
            N_POINT=2*N_POINT-1
            GOTO 3
          ENDIF
!         Use Richardson extrapolation to remove 
!         the leading error term:
          INTEGRAL_NEW=(4.0E+00_RealK*INTEGRAL_NEW-INTEGRAL_OLD)
     &      /3.0E+00_RealK
          B(I)=B(I)+INTEGRAL_NEW
!
        ENDDO
!
      ELSE
        U_1=1.0E+00_RealK/(ALPHA*THETA_HIGH*WAVELENGTH_LONG)
        U_2=1.0E+00_RealK/(ALPHA*THETA_HIGH*WAVELENGTH_SHORT)
        U_3=1.0E+00_RealK/(ALPHA*THETA_LOW*WAVELENGTH_LONG)
        U_4=1.0E+00_RealK/(ALPHA*THETA_LOW*WAVELENGTH_SHORT)
!
        DO I=0, N_DEG_FIT
!
          B(I)=0.0E+00_RealK
!
!         First integral:
          U(1)=U_1
          U(2)=0.5E+00_RealK*(U_1+U_2)
          U(3)=U_2
          INTEGRAL_OLD=0.0E+00_RealK
          N_POINT=3
4         DO J=1, N_POINT
            INTEGRAND(J)=FNC_21(U(J), I)
          ENDDO
          INTEGRAL_NEW=TRAPEZOID(N_POINT, U, INTEGRAND)
          IF (ABS(INTEGRAL_NEW-INTEGRAL_OLD).GT.TOL_THERMAL
     &       *ABS(INTEGRAL_NEW)) THEN
            INTEGRAL_OLD=INTEGRAL_NEW
            IF (2*N_POINT.GT.NPD_THERMAL_ABSCISSA) THEN
              WRITE(IU_ERR, '(/A, /A)')
     &          '** Error: Too many points are required to '
     &          //'evaluate the source function. '
     &          , 'increase NPD_THERMAL_ABSCISSA and re-compile.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            DO J=N_POINT, 2, -1
              U(2*J-1)=U(J)
              U(2*J-2)=0.5E+00_RealK*(U(J)+U(J-1))
            ENDDO
            N_POINT=2*N_POINT-1
            GOTO 4
          ENDIF
!         Use Richardson extrapolation to remove 
!         the leading error term:
          INTEGRAL_NEW=(4.0E+00_RealK*INTEGRAL_NEW-INTEGRAL_OLD)
     &      /3.0E+00_RealK
          B(I)=B(I)+INTEGRAL_NEW
!
!         Second integral:
          U(1)=U_2
          U(2)=0.5E+00_RealK*(U_2+U_3)
          U(3)=U_3
          INTEGRAL_OLD=0.0E+00_RealK
          N_POINT=3
5         DO J=1, N_POINT
            INTEGRAND(J)=FNC_22(U(J), I)
          ENDDO
          INTEGRAL_NEW=TRAPEZOID(N_POINT, U, INTEGRAND)
          IF (ABS(INTEGRAL_NEW-INTEGRAL_OLD).GT.TOL_THERMAL
     &       *ABS(INTEGRAL_NEW)) THEN
            INTEGRAL_OLD=INTEGRAL_NEW
            IF (2*N_POINT.GT.NPD_THERMAL_ABSCISSA) THEN
              WRITE(IU_ERR, '(/A, /A)')
     &          '** Error: Too many points are required to '
     &          //'evaluate the source function. '
     &          , 'increase NPD_THERMAL_ABSCISSA and re-compile.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            DO J=N_POINT, 2, -1
              U(2*J-1)=U(J)
              U(2*J-2)=0.5E+00_RealK*(U(J)+U(J-1))
            ENDDO
            N_POINT=2*N_POINT-1
            GOTO 5
          ENDIF
!         Use Richardson extrapolation to remove 
!         the leading error term:
          INTEGRAL_NEW=(4.0E+00_RealK*INTEGRAL_NEW-INTEGRAL_OLD)
     &      /3.0E+00_RealK
          B(I)=B(I)+INTEGRAL_NEW
!
!         Third integral:
          U(1)=U_3
          U(2)=0.5E+00_RealK*(U_3+U_4)
          U(3)=U_4
          INTEGRAL_OLD=0.0E+00_RealK
          N_POINT=3
6         DO J=1, N_POINT
            INTEGRAND(J)=FNC_23(U(J), I)
          ENDDO
          INTEGRAL_NEW=TRAPEZOID(N_POINT, U, INTEGRAND)
          IF (ABS(INTEGRAL_NEW-INTEGRAL_OLD).GT.TOL_THERMAL
     &       *ABS(INTEGRAL_NEW)) THEN
            INTEGRAL_OLD=INTEGRAL_NEW
            IF (2*N_POINT.GT.NPD_THERMAL_ABSCISSA) THEN
              WRITE(IU_ERR, '(/A, /A)')
     &          '** Error: Too many points are required to '
     &          //'evaluate the source function. '
     &          , 'increase NPD_THERMAL_ABSCISSA and re-compile.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            DO J=N_POINT, 2, -1
              U(2*J-1)=U(J)
              U(2*J-2)=0.5E+00_RealK*(U(J)+U(J-1))
            ENDDO
            N_POINT=2*N_POINT-1
            GOTO 6
          ENDIF
!         Use Richardson extrapolation to remove 
!         the leading error term:
          INTEGRAL_NEW=(4.0E+00_RealK*INTEGRAL_NEW-INTEGRAL_OLD)
     &      /3.0E+00_RealK
          B(I)=B(I)+INTEGRAL_NEW
!
        ENDDO
!
      ENDIF
      
!
!
!
      RETURN
!
!
!
      CONTAINS
!
!     Internal functions for evaluation of the Planckian integrals
!     in the thermal source functions.
!
      FUNCTION FNC_11(X, N) RESULT(Y)
        REAL  (RealK) :: X
        REAL  (RealK) :: Y
        INTEGER :: N
        Y=Q*X**3*(THETA_HIGH**(N+5)
     &    -1.0E+00_RealK/(WAVELENGTH_LONG
     &    *ALPHA*X)**(N+5))/(REAL(N+5, RealK)*(EXP(X)-1.0E+00_RealK))
      END FUNCTION
!
      FUNCTION FNC_12(X, N) RESULT(Y)
        REAL  (RealK) :: X
        REAL  (RealK) :: Y
        INTEGER :: N
        Y=Q*X**3*(THETA_HIGH**(N+5)-THETA_LOW**(N+5))
     &    /(REAL(N+5, RealK)*(EXP(X)-1.0E+00_RealK))
      END FUNCTION
!
      FUNCTION FNC_13(X, N) RESULT(Y)
        REAL  (RealK) :: X
        REAL  (RealK) :: Y
        INTEGER :: N
        Y=Q*X**3*(1./(WAVELENGTH_SHORT*ALPHA*X)**(N+5)
     &    -THETA_LOW**(N+5))/(REAL(N+5, RealK)*(EXP(X)-1.0E+00_RealK))
      END FUNCTION
!
      FUNCTION FNC_21(X, N) RESULT(Y)
        REAL  (RealK) :: X
        REAL  (RealK) :: Y
        INTEGER :: N
        Y=Q*X**3*(THETA_HIGH**(N+5)
     &    -1.0E+00_RealK/(WAVELENGTH_LONG
     &    *ALPHA*X)**(N+5))/(REAL(N+5, RealK)*(EXP(X)-1.0E+00_RealK))
      END FUNCTION
!
      FUNCTION FNC_22(X, N) RESULT(Y)
        REAL  (RealK) :: X
        REAL  (RealK) :: Y
        INTEGER :: N
        Y=Q*X**3
     &    *(1.0E+00_RealK/(WAVELENGTH_SHORT*ALPHA*X)**(N+5)
     &    -1.0E+00_RealK/(WAVELENGTH_LONG*ALPHA*X)**(N+5))
     &    /(REAL(N+5, RealK)*(EXP(X)-1.0E+00_RealK))
      END FUNCTION
!
      FUNCTION FNC_23(X, N) RESULT(Y)
        REAL  (RealK) :: X
        REAL  (RealK) :: Y
        INTEGER :: N
        Y=Q*X**3
     &    *(1.0E+00_RealK/(WAVELENGTH_SHORT*ALPHA*X)**(N+5)
     &    -THETA_LOW**(N+5))/(REAL(N+5, RealK)*(EXP(X)-1.0E+00_RealK))
      END FUNCTION
!
!
!
      END

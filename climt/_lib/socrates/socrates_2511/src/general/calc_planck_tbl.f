! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the Planck function table.
!
! Method:
!   A table in the temperature to the Planck function 
!   integrated across a spectral band is required. This is obtained
!   by explicit integration using Simpson's Rule.
!
!   This routine does basically the same as 'calc_planck_tbl' below.
!   However, but we have to integrate the product of the Planckian
!   and the filter function.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_PLANCK_FILT_TBL(IERR
     &  , N_THETA_POINTS, N_PTS
     &  , LAMBDA_MIN, LAMBDA_MAX
     &  , THETA_PLANCK_TBL
     &  , C1, C2
     &  , WEIGHTS
     &  , PLANCK_FILT_TBL)

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
     &    N_THETA_POINTS
!           Number of temperature points in table
     &  , N_PTS
!           Number of points in filter function     
      REAL  (RealK), INTENT(IN) ::      
     &    LAMBDA_MIN
!           Short wavelength of band
     &  , LAMBDA_MAX
!           Long wavelength of band     
     &  , THETA_PLANCK_TBL(0: NPD_THERMAL_COEFF-1)
!           Scaled temperature array
     &  , C1
!           Numerical constant
     &  , C2
!           Numerical constant
     &  , WEIGHTS(0:N_PTS-1)
!           Weights contained in the filter function

      REAL  (RealK), INTENT(OUT) ::
     &    PLANCK_FILT_TBL(0: NPD_THERMAL_COEFF-1)
!           Table of values for the band-integrated Planck function
     
! Local Variables
  
      INTEGER
     &    I, J
!           Looping variables

      REAL (Realk) ::
     &    EXPON
!          Exponential in Planckian
     &  , LAMBDA
!           Wavelength
     &  , D_LAMBDA
!            Wavelength increment
     &  ,  DUMMY(0:N_PTS-1)
!           Dummy variable     


! Step in wavelength
 
      D_LAMBDA=(LAMBDA_MAX-LAMBDA_MIN)/REAL(N_PTS-1, RealK)
      
!      print*,'D_LAMBDA=',D_LAMBDA
      
     
      DO J=0,N_THETA_POINTS
        
        DO I=0,N_PTS-1
          
          LAMBDA=LAMBDA_MIN+I*D_LAMBDA
          EXPON =EXP(C2/(LAMBDA*THETA_PLANCK_TBL(J))) 
          DUMMY(I)=WEIGHTS(I)*C1/(LAMBDA**5*(EXPON-1.0))
        ENDDO

        CALL SIMPSONS_RULE(N_PTS-1,D_LAMBDA,DUMMY,PLANCK_FILT_TBL(J))
        
!        print*,j,theta(j),planck_filt_tbl(j)
        
      ENDDO
      
      print*,''
      
      
      RETURN
      END
!-------------------------------------------------------------------------
!
!+ Subroutine to calculate the Planck function table.
!
! Method:
!   A table in the temperature to the Planck function 
!   integrated across a spectral band is required. This is obtained
!   by explicit numerical integration using Simpson's Rule.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_PLANCK_TBL(IERR
     &  , N_THETA_POINTS
     &  , WAVELENGTH_SHORT, WAVELENGTH_LONG
     &  , THETA_PLANCK_TBL
     &  , ALPHA, Q, PLANCK_TBL)

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
     &    N_THETA_POINTS
!           Degree of fit to function
      REAL  (RealK), INTENT(IN) ::
     &    WAVELENGTH_SHORT
!           Short wavelength of band
     &  , WAVELENGTH_LONG
!           Long wavelength of band
     &  , THETA_PLANCK_TBL(0: NPD_THERMAL_COEFF-1)
!           Scaled temperature array
     &  , ALPHA
!           Numerical constant
     &  , Q
!           Numerical constant
      REAL  (RealK), Intent(OUT) ::
     &    PLANCK_TBL(0: NPD_THERMAL_COEFF-1)
!           Table of values for the band-integrated Planck function
!
!     Local arguments.
      INTEGER
     &    N_POINT
!           Number of integration points
     &  , I
!           Loop variable
     &  , J
!           Loop variable
      REAL  (RealK) ::
     &    U_MIN
!           Limit of integral
     &  , U_MAX
!           Limit of integral
     &  , U
!           Local value of U
!
      REAL(RealK) :: D_U=1e-4_RealK
!           Step in integration over U parameter. The default value of
!           1e-4 should be small enough for all situations.
      REAL(RealK) :: D_U_USE
!           The actual integration step used to ensure an integer number
!           of steps in the interval.
      REAL(RealK), ALLOCATABLE :: DUMMY(:)
!           Dummy variable
!
!
!
!
!
!     
!     Loop over temperatures
      DO J=0,N_THETA_POINTS
!
!       Find minimum and maximum values for U at given THETA
        U_MIN=1.0E+00_RealK/(ALPHA*THETA_PLANCK_TBL(J)*WAVELENGTH_LONG)
        U_MAX=1.0E+00_RealK/(ALPHA*THETA_PLANCK_TBL(J)*WAVELENGTH_SHORT)
        N_POINT = NINT((U_MAX - U_MIN)/D_U) + 1
        D_U_USE = (U_MAX - U_MIN)/(N_POINT-1)
        
!       Calculate integrand
        ALLOCATE(DUMMY(0:N_POINT-1))
        DO I=0,N_POINT-1
          U = U_MIN + D_U_USE*I
          DUMMY(I)=U**3*EXP(-U)/(1 - EXP(-U))
        ENDDO
        
!       Integration over the transformed variable
!           U=(H_PLANCK*NU)/(K_BLOTZMANN*T_REF_THERMAL*THETA)
!       using Simpson's rule.
        CALL SIMPSONS_RULE(N_POINT-1,D_U_USE,DUMMY,PLANCK_TBL(J))
        DEALLOCATE(DUMMY)
      ENDDO
      PLANCK_TBL = PLANCK_TBL*Q*THETA_PLANCK_TBL**4
!
!
!
      END

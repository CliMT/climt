! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the gamma-function.
!
! Method:
!   The gamma function for the fractional part of the argument is
!   evaluated using Weierstrass's formula. The integral part is
!   handled by recurrence. A correction is applied to Weierstrass's
!   formula to reduce the influence of the truncation: this is
!   derived by imposing the condition that Gamma(1)=1. The
!   coefficient in this truncation could be precalculated, but it
!   depends on N_LIMIT and is therefore calculated here so that
!   N_LIMIT may be changed easily and consistently. With
!   N_LIMIT=10,000, the error should be less than one part in 10^7.
!   This routine is intended only for positive arguments and is
!   not especially efficient.
!
!- ---------------------------------------------------------------------
      function gamma_fnc(ierr, z)
!
!
!
!     Modules used
      USE realtype_rd
      USE error_pcf
      USE def_std_io_icf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      REAL  (RealK), Intent(IN) ::
     &    z
!           Point of evaluation
      REAL  (RealK) ::
     &    gamma_fnc
!           Gamma function
!
!     Local variables
      INTEGER
     &    n_limit
!           Limiting number of terms
     &  , z_int
!           Integer part of z
     &  , i
!           Loop varaible
      REAL  (RealK) ::
     &    exclude_zero
!           Region to exclude near 0.0
     &  , euler
!           Euler's constant
     &  , z_reduced
!           Reduced part of z
     &  , u
!           Temporary variable
     &  , v
!           Temporary variable
     &  , u_1
!           Temporary variable
     &  , v_1
!           Temporary variable
     &  , correction
!           Correction to allow for
!           Truncation
!
!     Define euler's constant
      DATA euler/5.77215665e-01_RealK/
!     Define numerical constants.
      DATA exclude_zero/5.0e-02_RealK/
      DATA n_limit/10000/
!
!
!     Reduce the argument to a value in (0,1)
      IF (z <= 0.0_RealK) THEN
        WRITE(iu_err, '(/a)') 
     &    '***error: the gamma function cannot '
     &    //'be evaluated below 0.0.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      z_int=int(z)
      z_reduced=z-z_int
!     There will be numerical difficulties near 0, 
!     so we exclude such cases.
      IF (z_reduced < exclude_zero) THEN
        z_reduced=z_reduced+1.0_RealK
        z_int=z_int-1
      ENDIF
      u=0.0_RealK
      u_1=0.0_RealK
      v=1.0_RealK
      v_1=1.0_RealK
      DO i=1, n_limit
        u=u+z_reduced/REAL(i, RealK)
        v=v*(1.0_RealK+z_reduced/REAL(i, RealK))
        u_1=u_1+1.0_RealK/REAL(i, RealK)
        v_1=v_1*(1.0_RealK+1.0_RealK/REAL(i, RealK))
      ENDDO
      correction=1.0_RealK+z_reduced**2
     &  *(1.0_RealK-EXP(-euler)/(v_1*EXP(-u_1)))
      gamma_fnc=1.0_RealK/(z_reduced*v*EXP(euler*z_reduced-u))
      gamma_fnc=gamma_fnc*correction
!
!     Convert from the reduced function to the full function.
      IF (z_int >= 1) THEN
        DO i=1, z_int
          gamma_fnc=gamma_fnc*(z-REAL(i, RealK))
        ENDDO
      ELSE if(z_int == -1) THEN
!       When z is small z_reduced will have been raised to avoid
!       ill-conditioning, so the conversion goes the other way.
        gamma_fnc=gamma_fnc/z
      ENDIF
!
!
!
      RETURN
      END

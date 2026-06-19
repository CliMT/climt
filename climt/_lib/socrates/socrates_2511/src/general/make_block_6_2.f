! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 6.
!
! Method:
!   A range of temperatures and whether to calculate a polynomial fit
!   or table are determined. If a polynomial, the degree of the fitting
!   polynomial is determined. Using SVD, a least squares fit
!   to the Planck function in each band over this range of
!   temperatures is made to determine the coefficients. If a table, the
!   number of data points in the table is determined. Using explicit
!   numerical integration, the Planck function in each band over this
!   range of temperatures is calculated and stored.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MAKE_BLOCK_6_2(IERR
     &  , N_BAND, WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG
     &  , L_EXCLUDE, N_BAND_EXCLUDE, INDEX_EXCLUDE
     &  , N_DEG_FIT, T_REF_THERMAL, THERMAL_COEFFICIENT
     &  , THETA_PLANCK_TBL, L_PRESENT_6, L_PLANCK_TBL
     &  )

      USE realtype_rd
      USE def_inst_flt
      USE rad_pcf
      USE rad_ccf
      USE dimensions_spec_ucf
      USE def_std_io_icf

      IMPLICIT NONE


! Dummy arguments.

      INTEGER, INTENT(INOUT) ::
     &    IERR
!           Error flag
     &  , N_DEG_FIT
!           Number of terms in fit
      LOGICAL, INTENT(IN) ::
     &    L_EXCLUDE
!           Flag for exclusion of specific bands
      INTEGER, INTENT(IN) ::
     &    N_BAND
!           Number of spectral bands
     &  , N_BAND_EXCLUDE(N_BAND)
!           Number of bands excluded from a specific region
     &  , INDEX_EXCLUDE(NPD_EXCLUDE, N_BAND)
!           Indices of excluded bands
      LOGICAL, INTENT(OUT) ::
     &    L_PRESENT_6
!           Presence flag for block
     &  , L_PLANCK_TBL
!           Flag for using a look-up table instead of a polynomial
      REAL  (RealK), INTENT(IN) ::
     &    WAVE_LENGTH_SHORT(N_BAND)
!           Shorter wavelengths of bands
     &  , WAVE_LENGTH_LONG(N_BAND)
!           Longer wavelengths of bands
      REAL  (RealK), INTENT(OUT) ::
     &    T_REF_THERMAL
!           Planckian reference temp
     &  , THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, N_BAND)
!           Coefficients for fits
      REAL(RealK), Intent(OUT) ::
     &    THETA_PLANCK_TBL(0: NPD_THERMAL_COEFF-1)
!           Holds the temperatures for the table

! Local variables.

      LOGICAL
     &    L_LAMBDA
!           Flag is set true if wavelength limits is filter
!           and spectral agree
     &  , L_EXIST
!           Check whether file exists
!
      INTEGER
     &    IOS
!           I/O error flag
     &  , IO_STATUS
     &  , I, J, K
!           Loop variable
     &  , N_COEFFICIENT
!           Number of coefficients
     &  , I_EXCLUDE
!           Number of excluded band

      REAL (KIND=KIND(1.0D0)), ALLOCATABLE ::
     &     WEIGHTS(:)
!           Weights contained in the filter function
     &  , LAMBDA(:)
!           Wavelengths corresponding to these weights
     &  , DUMMY(:)
!           Dummy_variable

      REAL  (KIND=KIND(1.0D0)) ::
     &    LAMBDA_MIN
!           Smallest Wavelength
     &  , LAMBDA_MAX
!           Largest Wavelength
     &  , D_LAMBDA
!           Wavelength Increment
     &  , T_FIT_LOW
!           Lowest fitted temperature
     &  , T_FIT_HIGH
!           Lowest fitted temperature
     &  , T_TEMP
!           Temporary temperature
     &  , C1
!           Numerical constant
     &  , C2
!           Numerical constant
     &  , EPS1,EPS2
!           Small numbers
     &  , A(0: NPD_THERMAL_COEFF-1, 0: NPD_THERMAL_COEFF-1)
!           Least squares coefficient matrix
     &  , A_EXCLUDE(0: NPD_THERMAL_COEFF-1, 0: NPD_THERMAL_COEFF-1)
!           Least squares coefficient matrix for excluded bands
     &  , B(0: NPD_THERMAL_COEFF-1)
!           Least squares forcing term
     &  , B_EXCLUDE(0: NPD_THERMAL_COEFF-1)
!           Least squares forcing term for excluded bands
     &  , W(0: NPD_THERMAL_COEFF-1)
!           Diagonal of svd decomposition
     &  , V(0: NPD_THERMAL_COEFF-1, 0: NPD_THERMAL_COEFF-1)
!           Least squares coefficient matrix
     &  , WRK(0: NPD_THERMAL_COEFF-1)
!           Work array
     &  , WMAX
!           Maximum diagonal element
     &  , THRESHOLD
!           Threshold for zeroing

      CHARACTER(LEN=1) :: L_PLANCK_TBL_CHAR
!           Response to fit/table question
      REAL(RealK) :: D_LOG_THETA
!           Step in temperature array

!     Set precision in calculating the thermal source.
      REAL (RealK), PARAMETER :: p_min_thermal_range=1.0e+00_RealK
!           Minium range of temperature

      TYPE  (StrFiltResp) :: filter
!   Instrumental response function

! Subroutines called:

      EXTERNAL
     &    CALC_THERMAL_COEFF_NUM, SVD_DECOMPOSE, BACK_SUBSTITUTE


! Should a polynomial fit or a table be used?

      WRITE(*, '(/A)') 'SHOULD A TABLE OR POLYNOMIAL FIT BE USED FOR '
     &   //'THE BAND-INTEGRATED PLANCK FUNCTION (T/P)?'
      DO
        READ(*, *, IOSTAT=IOS) L_PLANCK_TBL_CHAR
        IF (IOS /= 0) THEN
          WRITE(*, '(A)') '***ERROR: UNRECOGNIZED RESPONSE'
          WRITE(*, '(A)') 'PLEASE RE-ENTER.'
        ELSE
          IF (L_PLANCK_TBL_CHAR=='T' .OR. L_PLANCK_TBL_CHAR=='T') THEN
            L_PLANCK_TBL=.TRUE.
          ELSE
            L_PLANCK_TBL=.FALSE.
          ENDIF
          EXIT
        ENDIF
      ENDDO

! Obtain limits over which fit is required.

      WRITE(IU_STDOUT, '(/A)') 'ENTER THE RANGE OF TEMPERATURES OVER '
     &   //'WHICH THE TABLE/FIT IS REQUIRED.'
2     READ(IU_STDIN, *, IOSTAT=IOS) T_FIT_LOW, T_FIT_HIGH
      IF (IOS.NE.0) THEN
        WRITE(IU_ERR, '(A)') '+++ UNRECOGNIZED RESPONSE:'
        WRITE(IU_STDOUT, '(A)') 'PLEASE RE-ENTER.'
        GOTO 2
      ENDIF
      IF (T_FIT_LOW.GT.T_FIT_HIGH) THEN
        T_TEMP=T_FIT_LOW
        T_FIT_LOW=T_FIT_HIGH
        T_FIT_HIGH=T_TEMP
      ENDIF
      IF (T_FIT_LOW.LT.0.0E+00_RealK) THEN
        WRITE(IU_ERR, '(A)') '+++ ILLEGAL RESPONSE:'
        WRITE(IU_ERR, '(A)') 'A NEGATIVE ABSOLUTE TEMPERATURE WAS '
        WRITE(IU_STDOUT, '(A)') 'PLEASE RE-SPECIFY.'
        GOTO 2
      ENDIF
      IF (T_FIT_LOW.GT.(T_FIT_HIGH-P_MIN_THERMAL_RANGE)) THEN
        WRITE(IU_ERR, '(A)')
     &    '+++ RANGE OF TEMPERATURES IS TOO SMALL:'
        WRITE(IU_STDOUT, '(A)') 'PLEASE RE-ENTER.'
        GOTO 2
      ENDIF

! Obtain number of points in table or degree of fitting polynomial

      IF (L_PLANCK_TBL) THEN
        WRITE(IU_STDOUT, '(/A)')
     &    'ENTER THE NUMBER OF TEMPERATURE POINTS FOR THE TABLE.'
      ELSE
        WRITE(IU_STDOUT, '(/A)')
     &    'ENTER THE DEGREE OF THE FITTING POLYNOMIAL.'
      END IF
3     READ(IU_STDIN, *, IOSTAT=IOS) N_DEG_FIT
      IF (IOS.NE.0) THEN
        WRITE(IU_ERR, '(A)') '+++ UNRECOGNIZED RESPONSE'
        WRITE(IU_STDOUT, '(A)') 'PLEASE RE-ENTER.'
        GOTO 3
      ENDIF
      IF (N_DEG_FIT.LT.0) THEN
        WRITE(IU_ERR, '(A)') 'A NEGATIVE DEGREE IS UNACCEPTABLE: '
        WRITE(IU_STDOUT, '(A)') 'PLEASE RESPECIFY.'
        GOTO 3
      ENDIF
      IF (N_DEG_FIT.GT.NPD_THERMAL_COEFF-1) THEN
        WRITE(IU_ERR, '(A)') 'THIS DEGREE IS TOO HIGH: '
        WRITE(IU_STDOUT, '(A)') 'PLEASE RESPECIFY.'
        GOTO 3
      ENDIF

! Determine the number of thermal coefficients required

      IF (L_PLANCK_TBL) THEN
        N_DEG_FIT=N_DEG_FIT-1
      ENDIF
      N_COEFFICIENT=N_DEG_FIT+1

! The reference temperature can conveniently be taken as the
! geometric mean temperature of the range.

      T_REF_THERMAL=SQRT(T_FIT_LOW*T_FIT_HIGH)

! Scale the temperatures.

      T_FIT_LOW=T_FIT_LOW/T_REF_THERMAL
      T_FIT_HIGH=T_FIT_HIGH/T_REF_THERMAL

      C1=2.0*PI*H_PLANCK*C_LIGHT**2
      C2=H_PLANCK*C_LIGHT/(K_BOLTZMANN*T_REF_THERMAL)

! Construct logarithmic temperature array

      IF (L_PLANCK_TBL) THEN
        D_LOG_THETA = (LOG10(T_FIT_HIGH) - LOG10(T_FIT_LOW))/
     &    (N_DEG_FIT)
        DO I=0, N_DEG_FIT
          THETA_PLANCK_TBL(I) = 10**(LOG10(T_FIT_LOW) +
     &      D_LOG_THETA*I)
        ENDDO
      ENDIF

      DO I=1,N_BAND

! Read in filter function from file

        CALL READ_INSTRUMENT_RESPONSE_90(FILTER,IERR)

! and transform wavenumber back into wavelength. Here we need to
! be slightly careful since READ_INSTRUMENT_RESPONSE_90 reads in
! the wavelength transforms them into wavenumbers (m-1) and then
! reorders them in order of increasing wavenumber!

! Since all the files are in wavelength (and I prefer working in
! wavelength) we change back to wavelength.

        ALLOCATE(LAMBDA(FILTER%N_PTS))
        ALLOCATE(WEIGHTS(FILTER%N_PTS))

        DO J=1,FILTER%N_PTS
          LAMBDA(FILTER%N_PTS-J+1)=1.0/FILTER%WAVENUMBER(J)
        ENDDO
!
! Compare wavelength speciefied in the spectral file to those specified
! in the filter file.

        LAMBDA_MIN=WAVE_LENGTH_SHORT(I)
        LAMBDA_MAX=WAVE_LENGTH_LONG(I)

        D_LAMBDA=(LAMBDA_MAX-LAMBDA_MIN)/REAL(FILTER%N_PTS-1)
!         print*,'D_LAMBDA=',D_LAMBDA


        L_LAMBDA=.FALSE.

        EPS1=ABS(LAMBDA_MIN-MINVAL(LAMBDA))
        EPS2=ABS(LAMBDA_MAX-MAXVAL(LAMBDA))
        IF ((EPS1.LE.1.0E-08).AND.(EPS2.LE.1.0E-08)) THEN
          L_LAMBDA=.TRUE.
        ENDIF

! If the limits agree then interpolate filter function onto regular
! grid. This step is required in order to do the integration. If the
! limits do not agree the weights are all set to 1.

! EPS1 is used as a dummy variable

        IF (L_LAMBDA) THEN

          DO J=1,FILTER%N_PTS
             EPS1=(LAMBDA_MIN+(J-1)*D_LAMBDA)*1.0E+06
             CALL LINEAR_INTERPOLATION(LAMBDA*1.0E+06,
     &              FILTER%RESPONSE,
     &              FILTER%N_PTS,EPS1,
     &              WEIGHTS(FILTER%N_PTS-J+1))
          ENDDO

!          DO J=1,FILTER%N_PTS
!             WRITE(20,*) WEIGHTS(J),
!     &                     FILTER%RESPONSE(FILTER%N_PTS-J+1)
!          ENDDO

        ELSE
          WRITE(IU_STDOUT, '(A)') ' WAVELENGTH LIMITS IN THE '
     &      // 'SPECTRAL AND IN THE FILTER FILE DO NOT AGREE!'
          WRITE(IU_STDOUT, '(A)') ' ALL WEIGHTS HAVE BEEN SET TO 1'
          WEIGHTS=1.0
        ENDIF

!      WRITE(21,*) c1,c2

        IF (L_PLANCK_TBL) THEN

          CALL CALC_PLANCK_FILT_TBL(IERR
     &        , N_DEG_FIT, FILTER%N_PTS
     &        , LAMBDA_MIN, LAMBDA_MAX
     &        , THETA_PLANCK_TBL
     &        , C1, C2
     &        , WEIGHTS, THERMAL_COEFFICIENT(:,I))

          IF (L_EXCLUDE) THEN
            DO j=1, N_BAND_EXCLUDE(I)
              I_EXCLUDE=INDEX_EXCLUDE(J, I)
              CALL CALC_PLANCK_FILT_TBL(IERR
     &            , N_DEG_FIT, FILTER%N_PTS
     &            , WAVE_LENGTH_SHORT(I_EXCLUDE)
     &            , WAVE_LENGTH_LONG(I_EXCLUDE)
     &            , THETA_PLANCK_TBL
     &            , C1, C2
     &            , WEIGHTS
     &            , B_EXCLUDE)
              IF (IERR /= I_NORMAL) RETURN
              THERMAL_COEFFICIENT(:,I)=
     &          THERMAL_COEFFICIENT(:,I)-B_EXCLUDE(:)
            ENDDO
          ENDIF

        ELSE

          CALL CALC_THERMAL_COEFF_NUM(IERR
     &        , N_DEG_FIT, FILTER%N_PTS
     &        , LAMBDA_MIN, LAMBDA_MAX
     &        , T_FIT_LOW, T_FIT_HIGH
     &        , C1, C2
     &        , WEIGHTS,A, B)
          IF (IERR.NE.I_NORMAL) RETURN

   ! If regions of the band are excluded, the corresponding regions
   ! must be removed from B.

          IF (L_EXCLUDE) THEN
            DO j=1, N_BAND_EXCLUDE(I)
              I_EXCLUDE=INDEX_EXCLUDE(J, I)
              CALL CALC_THERMAL_COEFF_NUM(IERR
     &            , N_DEG_FIT, FILTER%N_PTS
     &            , WAVE_LENGTH_SHORT(I_EXCLUDE)
     &            , WAVE_LENGTH_LONG(I_EXCLUDE)
     &            , T_FIT_LOW, T_FIT_HIGH
     &            , C1, C2
     &            , WEIGHTS
     &            , A_EXCLUDE, B_EXCLUDE)
              IF (IERR /= I_NORMAL) RETURN
              DO K=0, N_DEG_FIT
                B(K)=B(K)-B_EXCLUDE(K)
              ENDDO
            ENDDO
          ENDIF
        ENDIF


! Perform an SVD decomposition of the matrix.

        CALL SVD_DECOMPOSE(IERR
     &      , A, N_COEFFICIENT, N_COEFFICIENT
     &      , NPD_THERMAL_COEFF, NPD_THERMAL_COEFF
     &      , W, V, WRK
     &      )
        IF (IERR.NE.I_NORMAL) RETURN

! Zero very small terms of the diagonal matrix.

        WMAX=0.0E+00_RealK
        DO J=0, N_DEG_FIT
          IF (W(J).GT.WMAX) THEN
              WMAX=W(J)
          ENDIF
        ENDDO
        THRESHOLD=1.0E+03_RealK*EPSILON(WMAX)*WMAX
        DO J=0, N_DEG_FIT
           IF (W(J).LT.THRESHOLD) THEN
              W(J)=0.0E+00_RealK
           ENDIF
        ENDDO

! Solve the equations by back-substitution.

        CALL BACK_SUBSTITUTE(A, W, V
     &      , N_COEFFICIENT, N_COEFFICIENT
     &      , NPD_THERMAL_COEFF, NPD_THERMAL_COEFF
     &      , B, THERMAL_COEFFICIENT(0,i), WRK
     &      )

!       do k=0,npd_thermal_coeff-1
!          print*,i,k,THERMAL_COEFFICIENT(k,I)
!       enddo
      ENDDO

      DEALLOCATE(LAMBDA)
      DEALLOCATE(WEIGHTS)

!     The data are now correctly assembled so the presence flag is set.

      L_PRESENT_6=.TRUE.

      RETURN
      END

      SUBROUTINE LINEAR_INTERPOLATION(LAMBDA,WEIGHTS,N_FILTER,L,
     &                                               INTERPOLATED)

!********************************************
! Linear Interpolation: y=c1*y(i)+c2*y(i+1) *
!********************************************

      IMPLICIT NONE

      INTEGER
     &    N_FILTER
!            Number of weights in filter function
     &  , NS
!            Index
     &  , I
!            Loop variable

      REAL(KIND=KIND(1.0D0)), DIMENSION (0:N_FILTER-1) ::
     &    WEIGHTS
!           Filter Function
     &  , LAMBDA
!           Wavelength

      REAL(KIND=KIND(1.0D0))
     &    L
!            Wavelength at which we require the new weight
     &  , INTERPOLATED
!            Interpolated value
     &  , A1,A2
!            Interpolation Constants
     &  , DIF,DIFT
! Distance between 'l' and 'lambda'

      NS = 0
      DIF=ABS(L-LAMBDA(0))

! Find the index of 'lambda' which is closest to lambda_reg

      DO I=0,N_FILTER-1
         DIFT=ABS(L-LAMBDA(i))
         IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
         ENDIF
      ENDDO

! Interpolate and return value

      A1=(LAMBDA(NS+1)-L)/(LAMBDA(NS+1)-LAMBDA(NS))
      A2=(L-LAMBDA(NS)  )/(LAMBDA(NS+1)-LAMBDA(NS))

      INTERPOLATED=A1*WEIGHTS(NS)+A2*WEIGHTS(NS+1)

      IF (INTERPOLATED.LT.0.0) THEN
         INTERPOLATED=0.0
      ENDIF

      END

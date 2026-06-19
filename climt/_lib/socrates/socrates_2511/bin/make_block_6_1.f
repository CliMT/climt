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
      SUBROUTINE make_block_6_1(ierr
     &  , n_band, wave_length_short, wave_length_long
     &  , l_exclude, n_band_exclude, index_exclude
     &  , n_deg_fit, t_ref_thermal, thermal_coefficient
     &  , theta_planck_tbl, l_present_6, l_planck_tbl
     &  )

      USE realtype_rd
      USE rad_pcf
      USE rad_ccf
      USE dimensions_spec_ucf
      USE def_std_io_icf

      IMPLICIT NONE


!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , n_deg_fit
!           Number of terms in fit
      LOGICAL, Intent(IN) ::
     &    l_exclude
!           Flag for exclusion of specific bands
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of spectral bands
     &  , n_band_exclude(n_band)
!           Number of bands excluded from a specific region
     &  , index_exclude(npd_exclude, n_band)
!           Indices of excluded bands
      LOGICAL, Intent(OUT) ::
     &    l_present_6,
!           Presence flag for block
     &    l_planck_tbl
!           Flag for using a look-up table instead of a polynomial
      REAL  (RealK), Intent(IN) ::
     &    wave_length_short(n_band)
!           Shorter wavelengths of bands
     &  , wave_length_long(n_band)
!           Longer wavelengths of bands
      REAL  (RealK), Intent(OUT) ::
     &    t_ref_thermal
!           Planckian reference temp
     &  , thermal_coefficient(0: npd_thermal_coeff-1, n_band)
!           Coefficients for fits
      REAL(RealK), Intent(OUT) ::
     &    theta_planck_tbl(0: npd_thermal_coeff-1)
!           Holds the temperatures for the table
!
!     Local variables.
!
      INTEGER
     &    ios
!           I/O error flag
     &  , i
!           Loop variable
     &  , j
!           Loop variable
     &  , k
!           Loop variable
     &  , n_coefficient
!           Number of coefficients
     &  , i_exclude
!           Number of excluded band
      REAL  (RealK) ::
     &    t_fit_low
!           Lowest fitted temperature
     &  , t_fit_high
!           Lowest fitted temperature
     &  , t_temp
!           Temporary temperature
     &  , alpha
!           Numerical constant
     &  , q
!           Numerical constant
     &  , a(0: npd_thermal_coeff-1, 0: npd_thermal_coeff-1)
!           Least squares coefficient matrix
     &  , b(0: npd_thermal_coeff-1)
!           Least squares forcing term
     &  , a_exclude(0: npd_thermal_coeff-1, 0: npd_thermal_coeff-1)
!           Least squares coefficient matrix for excluded bands
     &  , b_exclude(0: npd_thermal_coeff-1)
!           Least squares forcing term for excluded bands
     &  , w(0: npd_thermal_coeff-1)
!           Diagonal of svd decomposition
     &  , v(0: npd_thermal_coeff-1, 0: npd_thermal_coeff-1)
!           Least squares coefficient matrix
     &  , wrk(0: npd_thermal_coeff-1)
!           Work array
     &  , wmax
!           Maximum diagonal element
     &  , threshold
!           Threshold for zeroing
!
      CHARACTER(LEN=1) :: l_planck_tbl_char
!           Response to fit/table question
      REAL(RealK) :: d_log_theta
!           Step in temperature array

!     Set precision in calculating the thermal source.
      REAL (RealK), PARAMETER :: p_min_thermal_range=1.0e+00_RealK
!           Minium range of temperature


!     Should a polynomial fit or a table be used?
      WRITE(*, '(/A)') 'Should a Table or Polynomial fit be used for '
     &   //'the band-integrated Planck function (T/P)?'
      DO
        READ(*, *, IOSTAT=ios) l_planck_tbl_char
        IF (ios /= 0) THEN
          WRITE(*, '(A)') '***error: unrecognized response'
          WRITE(*, '(A)') 'Please re-enter.'
        ELSE
          IF (l_planck_tbl_char=='T' .OR. l_planck_tbl_char=='t') THEN
            l_planck_tbl=.TRUE.
          ELSE
            l_planck_tbl=.FALSE.
          ENDIF
          EXIT
        ENDIF
      ENDDO
!
!     Obtain limits over which fit is required.
      WRITE(iu_stdout, '(/a)') 'enter the range of temperatures over '
     &   //'which the table/fit is required.'
1     read(iu_stdin, *, iostat=ios) t_fit_low, t_fit_high
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ unrecognized response:'
        WRITE(iu_stdout, '(a)') 'please re-enter.'
        goto 1
      ENDIF
      IF (t_fit_low > t_fit_high) THEN
        t_temp=t_fit_low
        t_fit_low=t_fit_high
        t_fit_high=t_temp
      ENDIF
      IF (t_fit_low < 0.0_RealK) THEN
        WRITE(iu_err, '(a)') '+++ illegal response:'
        WRITE(iu_err, '(a)') 'a negative absolute temperature was '
        WRITE(iu_stdout, '(a)') 'please re-specify.'
        goto 1
      ENDIF
      IF (t_fit_low > (t_fit_high-p_min_thermal_range)) THEN
        WRITE(iu_err, '(a)')
     &    '+++ range of temperatures is too small:'
        WRITE(iu_stdout, '(a)') 'please re-enter.'
        goto 1
      ENDIF
!
!     Obtain number of points in table or degree of fitting polynomial
      IF (l_planck_tbl) THEN
        WRITE(iu_stdout, '(/a)')
     &    'enter the number of temperature points for the table.'
      ELSE
        WRITE(iu_stdout, '(/a)')
     &    'enter the degree of the fitting polynomial.'
      END IF
!
2     read(iu_stdin, *, iostat=ios) n_deg_fit
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ unrecognized response'
        WRITE(iu_stdout, '(a)') 'please re-enter.'
        goto 2
      ENDIF
      IF (n_deg_fit < 0) THEN
        WRITE(iu_err, '(a)') 'a negative degree is unacceptable: '
        WRITE(iu_stdout, '(a)') 'please respecify.'
        goto 2
      ENDIF
      IF (n_deg_fit > npd_thermal_coeff-1) THEN
        WRITE(iu_err, '(a)') 'this degree is too high: '
        WRITE(iu_stdout, '(a)') 'please respecify.'
        goto 2
      ENDIF
!
      IF (l_planck_tbl) THEN
        n_deg_fit=n_deg_fit-1
      ENDIF
      n_coefficient=n_deg_fit+1
!
!     The reference temperature can conveniently be taken as the
!     geometric mean temperature of the range.
      t_ref_thermal=sqrt(t_fit_low*t_fit_high)
!
!     Scale the temperatures.
      t_fit_low=t_fit_low/t_ref_thermal
      t_fit_high=t_fit_high/t_ref_thermal
!
      alpha=k_boltzmann*t_ref_thermal/(h_planck*c_light)
      q=2.0_RealK*pi*h_planck*alpha**4*c_light**2
!
!     Construct logarithmic temperature array
      IF (l_planck_tbl) THEN
        d_log_theta = (LOG10(t_fit_high) - LOG10(t_fit_low))/
     &    (n_deg_fit)
        DO i=0, n_deg_fit
          theta_planck_tbl(i) = 10**(LOG10(t_fit_low) +
     &      d_log_theta*i)
        ENDDO
      ENDIF
!
!     For each band a table or a least squares fit in the scaled
!     temperature is performed.
      DO i=1, n_band
        IF (l_planck_tbl) THEN
          CALL calc_planck_tbl(ierr
     &      , n_deg_fit
     &      , wave_length_short(i), wave_length_long(i)
     &      , theta_planck_tbl
     &      , alpha, q
     &      , thermal_coefficient(:,i))
!
!         If regions of the band are excluded, the corresponding
!         regions must be removed from the Planck table.
          IF (l_exclude) THEN
            DO j=1, n_band_exclude(i)
              i_exclude=index_exclude(j, i)
              CALL calc_planck_tbl(ierr
     &          , n_deg_fit
     &          , wave_length_short(i_exclude)
     &          , wave_length_long(i_exclude)
     &          , theta_planck_tbl 
     &          , alpha, q
     &          , b_exclude)
              IF (ierr /= i_normal) RETURN
              thermal_coefficient(:,i)=
     &          thermal_coefficient(:,i)-b_exclude(:)
            ENDDO
          ENDIF
        ELSE
          CALL calc_thermal_coeff(ierr
     &      , n_deg_fit
     &      , wave_length_short(i), wave_length_long(i)
     &      , t_fit_low, t_fit_high
     &      , alpha, q
     &      , a, b)
          IF (ierr /= i_normal) RETURN
!
!         If regions of the band are excluded, the corresponding
!         regions must be removed from B.
          IF (l_exclude) THEN
            DO j=1, n_band_exclude(i)
              i_exclude=index_exclude(j, i)
              CALL calc_thermal_coeff(ierr
     &          , n_deg_fit
     &          , wave_length_short(i_exclude)
     &          , wave_length_long(i_exclude)
     &          , t_fit_low, t_fit_high
     &          , alpha, q
     &          , a_exclude, b_exclude)
              IF (ierr /= i_normal) RETURN
              DO k=0, n_deg_fit
                b(k)=b(k)-b_exclude(k)
              ENDDO
            ENDDO
          ENDIF
!
!       Perform an SVD decomposition of the matrix.
          CALL svd_decompose(ierr
     &      , a, n_coefficient, n_coefficient
     &      , npd_thermal_coeff, npd_thermal_coeff
     &      , w, v, wrk
     &    )
          IF (ierr /= i_normal) RETURN
  !
  !       Zero very small terms of the diagonal matrix.
          wmax=0.0_RealK
          DO j=0, n_deg_fit
            IF (w(j) > wmax) THEN
              wmax=w(j)
            ENDIF
          ENDDO
          threshold=1.0e+03_RealK*epsilon(wmax)*wmax
          DO j=0, n_deg_fit
            IF (w(j) < threshold) THEN
              w(j)=0.0_RealK
            ENDIF
          ENDDO      
!
!       Solve the equations by back-substitution.
          CALL back_substitute(a, w, v
     &      , n_coefficient, n_coefficient
     &      , npd_thermal_coeff, npd_thermal_coeff
     &      , b, thermal_coefficient(0, i), wrk
     &      )
!
  
!        do k=0,npd_thermal_coeff-1
!          print*,i,k,THERMAL_COEFFICIENT(k,I)
!        enddo
  
        ENDIF
!
      ENDDO
!
!     The data are now correctly assembled so the presence flag is set.
      l_present_6=.true.
!
!
!
      RETURN
      END

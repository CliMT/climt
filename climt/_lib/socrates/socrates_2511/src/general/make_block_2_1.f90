! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 2.
!
! Method:
!   A solar spectrum, giving the irradiance against wavelength,
!   is read in. The total flux in the spectrum is determined,
!   using a Rayleigh-Jeans distribution for the far infra-red.
!   The fraction of this flux in each band is calculated. The
!   band will probably not cover the whole spectrum. If 
!   simulating observations this is as required, since the
!   spectral cut-offs of the observing instruments must be 
!   matched. In atmospheric models the whole flux should be
!   calculated and there is therefore an option to enhance
!   the fluxes in the outside bands to include the full flux.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_2_1(Sp, Sol, filter, l_filter, l_enhance, l_verbose, ierr)

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_normal, i_err_range
  USE def_spectrum, ONLY: StrSpecData
  USE def_solarspec, ONLY: StrSolarSpec
  USE def_inst_flt, ONLY: StrFiltResp
  USE spline_evaluate_mod, ONLY: spline_evaluate

  IMPLICIT NONE


! Dummy arguments
  TYPE (StrSpecData), INTENT(INOUT) :: Sp
!   Spectral file to be assigned
  TYPE (StrSolarSpec), INTENT(IN) :: Sol
!   Solar spectrum
  TYPE (StrFiltResp), INTENT(IN) :: filter
!   Instrumental response function
  LOGICAL, INTENT(IN) :: l_filter
!   Weight with a filter function
  LOGICAL, INTENT(IN) :: l_enhance
!   Enhance outer bands
  LOGICAL, INTENT(IN) :: l_verbose
!   Print diagnostic output
  INTEGER, INTENT(INOUT) :: ierr
!   Error flag

! Local variables
  INTEGER :: i_short
!       Beginning of band
  INTEGER :: i_long
!       End of band
  INTEGER :: n_points_band
!       Number of points within band
  INTEGER :: i, j
!       Loop variable
  INTEGER :: i_first_band
!       Index of shortest wavelength band
  INTEGER :: i_last_band
!       Index of longest wavelength band
  REAL (RealK) :: total_solar_irradiance
!       Total irradiance
  REAL (RealK), ALLOCATABLE :: x(:)
!       Absicissae of wavelength
  REAL (RealK), ALLOCATABLE :: y(:)
!       Ordinates of spectrum
  REAL (RealK), ALLOCATABLE :: bw(:)
!       Bandwidth of spectral bands
  REAL (RealK) :: irradiance_tail
!       Irradiance in tail
  REAL (RealK) :: irradiance_tail_band
!       Irradiance in Rayleigh-Jeans tail in the band
  REAL (RealK) :: wave_length_begin
!       Beginning of range
  REAL (RealK) :: wave_length_end
!       End of range
  REAL (RealK) :: sol_wavelength_begin
!       Beginning of solar spectrum
  REAL (RealK) :: sol_wavelength_end
!       End of solar spectrum

  REAL (RealK) :: response_0
!       Value of instrument response at given wavenumber

! Functions called:
  REAL (RealK), EXTERNAL :: trapezoid
!       Integrating function
  REAL (RealK), EXTERNAL :: rayleigh_jeans_tail
!       Tail of rayleigh jeans function
  REAL (RealK), EXTERNAL :: solar_intensity
!       Solar intensity at a given wavelength


! Calculate the total integrated irradiance of this spectrum.
  SELECT CASE (Sol%l_binned)
  CASE (.FALSE.)
    total_solar_irradiance = &
      trapezoid(Sol%n_points, Sol%wavelength, Sol%irrad)
    sol_wavelength_begin = Sol%wavelength(1)
    sol_wavelength_end   = Sol%wavelength(Sol%n_points)
  CASE (.TRUE.)
    total_solar_irradiance = SUM(Sol%bandsize * Sol%irrad)
    sol_wavelength_begin = Sol%bandbnds(1,1)
    sol_wavelength_end   = Sol%bandbnds(2,Sol%n_points)
  END SELECT  
  
  IF (l_verbose) WRITE(*, '(a,f21.9)') &
    'Total irradiance of solar spectrum = ', total_solar_irradiance

! Add on the tail of the distribution using a rayleigh-jeans law.
  irradiance_tail = &
    rayleigh_jeans_tail(Sol, sol_wavelength_end)
  total_solar_irradiance = total_solar_irradiance+irradiance_tail
  IF (l_verbose) WRITE(*, '(a,f16.9)') &
    'Total irradiance of Rayleigh-Jeans tail = ', irradiance_tail

! Find the position of the shortest and longest wavelength bands
  IF (l_enhance) THEN
    i_first_band=1
    i_last_band=Sp%Basic%n_band
    DO i=1, Sp%Basic%n_band
      IF (Sp%Basic%wavelength_short(i) < &
          Sp%Basic%wavelength_short(i_first_band)) i_first_band=i
      IF (Sp%Basic%wavelength_long(i) > &
          Sp%Basic%wavelength_long(i_last_band)) i_last_band=i
    ENDDO
  ENDIF

! For each band find the arrays of wavelength and irradiance.
  DO i=1, Sp%Basic%n_band

!   Find the points of the spectrum just within the band.
    wave_length_begin = Sp%Basic%wavelength_short(i)
    wave_length_end   = Sp%Basic%wavelength_long(i)
    IF (l_enhance .AND. i == i_first_band) THEN
      wave_length_begin = min(sol_wavelength_begin, &
                              wave_length_begin)
    END IF
    IF (l_enhance .AND. i == i_last_band) THEN
      wave_length_end   = max(sol_wavelength_end, &
                              wave_length_end)
    END IF

    IF (wave_length_begin > sol_wavelength_end) THEN
!     The whole band is within the Rayleigh-Jeans tail so we set directly
      Sp%Solar%solar_flux_band(i)= &
        (rayleigh_jeans_tail(Sol, wave_length_begin) - &
         rayleigh_jeans_tail(Sol, wave_length_end)) &
        /total_solar_irradiance
    ELSE
      IF (wave_length_end > sol_wavelength_end) THEN
!       Part of the band is in the Rayleigh-Jeans tail
        irradiance_tail_band = &
          irradiance_tail - rayleigh_jeans_tail(Sol, wave_length_end)
        wave_length_end = sol_wavelength_end
      ELSE
        irradiance_tail_band = 0.0_RealK
      END IF


      SELECT CASE (Sol%l_binned)
      CASE (.FALSE.)
        CALL inner_bracket(ierr, wave_length_begin, wave_length_end, &
          Sol%n_points, Sol%wavelength, i_short, i_long)
        IF (ierr /= i_normal) THEN
          IF (l_verbose) WRITE(*, '(a)') 'Error in call to inner_bracket'
          RETURN
        END IF
        n_points_band = 3+i_long-i_short
        ALLOCATE(x(n_points_band))
        ALLOCATE(y(n_points_band))
        x(1) = wave_length_begin
        y(1) = solar_intensity(x(1), Sol)
        x(n_points_band) = wave_length_end
        y(n_points_band) = solar_intensity(x(n_points_band), Sol)
        DO j=i_short, i_long
          x(j-i_short+2)=Sol%wavelength(j)
          y(j-i_short+2)=Sol%irrad(j)
        ENDDO
        
      CASE (.TRUE.)
        DO i_short=1, Sol%n_points
           IF (wave_length_begin <  Sol%bandbnds(2,i_short)) EXIT
        END DO
        DO i_long=Sol%n_points,1,-1
           IF (wave_length_end   >  Sol%bandbnds(1,i_long)) EXIT
        END DO
        n_points_band = 1+i_long-i_short 
        ALLOCATE(x(n_points_band))
        ALLOCATE(y(n_points_band))
        ALLOCATE(bw(n_points_band))
        DO j=i_short, i_long
          y(j-i_short+1)=Sol%irrad(j)
          IF (i_short == i_long) THEN
            bw(1) = wave_length_end - wave_length_begin
            x(1)  = wave_length_begin + (bw(1)/2.0)
          ELSE IF (j == i_short) THEN
            bw(j-i_short+1) = Sol%bandbnds(2,j) - wave_length_begin
            x(j-i_short+1)  = wave_length_begin + (bw(j-i_short+1)/2.0)
          ELSE IF (j == i_long) THEN
            bw(j-i_short+1) = wave_length_end - Sol%bandbnds(1,j)
            x(j-i_short+1)  = wave_length_end - (bw(j-i_short+1)/2.0)
          ELSE
            bw(j-i_short+1)=Sol%bandsize(j)
            x(j-i_short+1)=Sol%wavelength(j)
          ENDIF
        ENDDO
      END SELECT

!     Weight the irradiance values with the filter function
      IF (l_filter) THEN
        DO j=1, n_points_band
          CALL spline_evaluate(ierr, filter%n_pts, &
            filter%wavenumber, filter%response, filter%d2_response, &
            1.0_RealK/x(j), response_0)
          IF (ierr == i_err_range) THEN
!           The filter function is taken to be 0 outside 
!           the explicit range. We therefore zero the response and
!           recover from the error.
            response_0 = 0.0_RealK
            ierr = i_normal
          ENDIF
          IF (ierr /= i_normal) THEN
            IF (l_verbose) WRITE(*, '(a)') 'Error in call to spline_evaluate'
            RETURN
          END IF
          y(j) = y(j) * response_0
        ENDDO
      ENDIF

      SELECT CASE (Sol%l_binned)
      CASE (.FALSE.)
!       Integrate across the band and normalize by the total irradiance.
        Sp%Solar%solar_flux_band(i)=(trapezoid(n_points_band, x, y) &
          +irradiance_tail_band)/total_solar_irradiance
        DEALLOCATE(x,y)
      CASE (.TRUE.)
!       Integrate across the band and normalize by the total irradiance.
        Sp%Solar%solar_flux_band(i)=(SUM(bw * y) &
          +irradiance_tail_band)/total_solar_irradiance
        DEALLOCATE(x,y,bw)
      END SELECT
      
    ENDIF
  ENDDO

! Add the tail to the last band if required. The setting of limits
! above dealt with the contributions within the region of solar data.
  IF (l_enhance) THEN
    IF (Sp%Basic%wavelength_long(i_last_band) > &
        Sol%wavelength(Sol%n_points)) THEN
      irradiance_tail= &
        rayleigh_jeans_tail(Sol, Sp%Basic%wavelength_long(i_last_band))
    END IF
    Sp%Solar%solar_flux_band(i_last_band) = &
      Sp%Solar%solar_flux_band(i_last_band) &
      + irradiance_tail/total_solar_irradiance
  ENDIF

! Reduce the solar flux in bands which contain exclusions if necessary.
  IF (Sp%Basic%l_present(14)) THEN
    DO i=1, Sp%Basic%n_band
      DO j=1, Sp%Basic%n_band_exclude(i)
        Sp%Solar%solar_flux_band(i) = Sp%Solar%solar_flux_band(i) &
          - Sp%Solar%solar_flux_band(Sp%Basic%index_exclude(j, i))
      ENDDO
    ENDDO
  ENDIF

  Sp%Basic%l_present(2)=.TRUE.

  IF (l_verbose) WRITE(*, '(a,f16.9)') 'Sum of normalised flux = ', &
    SUM(Sp%Solar%solar_flux_band(1:Sp%Basic%n_band))

END SUBROUTINE make_block_2_1

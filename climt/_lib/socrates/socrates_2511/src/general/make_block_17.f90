! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to make spectral blocks of type 17.
!
! Description:
!   This routine creates a look-up table of spectral variability data
!
!------------------------------------------------------------------------------
SUBROUTINE make_block_17(Sp, Sol, ierr)

  USE netcdf
  USE realtype_rd, ONLY: RealK
  USE def_spectrum, ONLY: StrSpecData, StrSpecVar
  USE def_solarspec, ONLY: StrSolarSpec
  USE def_refract, ONLY: StrRefract
  USE def_inst_flt, ONLY: StrFiltResp
  USE missing_data_mod, ONLY: imdi
  USE rad_pcf, ONLY: ip_rayleigh_total

  IMPLICIT NONE

  TYPE(StrSpecData), Intent(INOUT) :: Sp
!   Spectral file data
  TYPE (StrSolarSpec), Intent(INOUT) :: Sol
!   Mean Solar spectrum
  INTEGER, INTENT(INOUT) :: ierr
!   Error flag

  TYPE(StrSpecData) :: SubSp
!   Temporary spectral file data for sub-bands
  TYPE(StrSpecVar) :: SpVarTmp
!   Temporary spectral variability data
  TYPE (StrSolarSpec) :: VSol, VSolMean
!   Varying Solar spectrum
  TYPE (StrRefract) :: Refract
!   Refractive index of atmosphere
  TYPE (StrFiltResp) :: filter
!   Instrumental response function

  INTEGER :: ios
!   Reading error flag
  INTEGER :: i, j, k, jj
!   Loop variables
  INTEGER :: i_format
!   Format of solar spectrum data file
  INTEGER :: iu_solar
!   Unit number for the solar spectrum data file
  LOGICAL :: l_count
!   Flag for counting points
  REAL (RealK) :: scale_wv, scale_irr
!   Scaling for wavelength and irradiance to correct units

  CHARACTER (LEN=1) :: char
  CHARACTER (LEN=80) :: line
  CHARACTER (LEN=256) :: cmip6_file
  CHARACTER (LEN=4) :: dim_name
  CHARACTER (LEN=3) :: frequency
  INTEGER :: ncid, varid, dimid_time, dimid_wlen, time_len, wlen_len
  INTEGER :: band, sub_band, number_term, index_absorb
  INTEGER :: sub_bands(Sp%Dim%nd_band), band_sort(Sp%Dim%nd_band)
  INTEGER :: n_times, yearstart=imdi, monthstart, daystart
  INTEGER, ALLOCATABLE :: calyear(:), calmonth(:), calday(:), seconds(:)
  LOGICAL :: l_monthly
  REAL (RealK), ALLOCATABLE :: tsi(:), ssi(:,:), wbinsize(:), wbinbnds(:,:)
  REAL (RealK) :: wavelength(2, &
                             MAX(Sp%Dim%nd_k_term, Sp%Dim%nd_sub_band_gas), &
                             Sp%Dim%nd_band)
  REAL (RealK) :: wave_inc

  INTERFACE
    SUBROUTINE map_heap_func(a, map)
      USE realtype_rd, ONLY: RealK
      REAL(RealK), INTENT(IN), DIMENSION(:) :: a
      INTEGER, INTENT(OUT), DIMENSION(:) :: map
    END SUBROUTINE map_heap_func
  END INTERFACE

  IF ( .NOT. Sp%Basic%l_present(17) ) THEN
    sub_bands=0
    sub_bands(1:Sp%Basic%n_band)=1
    DO
      WRITE(*, '(/a)') &
        'Enter band to be sub-divided (0 to finish, -1 to use gas sub-bands): '
      READ(*, *, IOSTAT=ios) band
      IF (band == 0) EXIT
      IF (band == -1) THEN
        DO band=1, Sp%Basic%n_band
          IF (Sp%Gas%n_band_absorb(band) > 0) THEN
            index_absorb = Sp%Gas%index_absorb(1, band)
            IF (Sp%Gas%n_sub_band_gas(band, index_absorb) > 1) THEN
              sub_bands(band) = Sp%Gas%n_sub_band_gas(band, index_absorb)
              wavelength(:, 1:sub_bands(band), band) = &
                Sp%Gas%wavelength_sub_band(:, 1:sub_bands(band), &
                                           band, index_absorb)
            END IF
          END IF
        END DO
        EXIT
      END IF
      index_absorb = Sp%Gas%index_absorb(1, band)
      IF (Sp%Gas%n_sub_band_gas(band, index_absorb) > 1) THEN
        sub_bands(band) = Sp%Gas%n_sub_band_gas(band, index_absorb)
        wavelength(:, 1:sub_bands(band), band) = &
          Sp%Gas%wavelength_sub_band(:, 1:sub_bands(band), band, index_absorb)
      ELSE
        number_term = Sp%Gas%i_band_k(band, index_absorb)
        sub_bands(band)=number_term
        WRITE(*, '(a,i5,a)') &
          'There are ', number_term, ' major gas k-terms in this band.'
        WRITE(*, '(a)') 'Do you want to divide equally in wavelength (E),'
        WRITE(*, '(a)') 'provide band limits (L) for each k-term,'
        WRITE(*, '(a)') 'or use a single sub-band for the band (S)?'
        READ(*, '(a)') char
        IF ( (char.eq.'s').OR.(char.eq.'S') ) THEN
          sub_bands(band)=1
        ELSE IF ( (char.eq.'e').OR.(char.eq.'E') ) THEN
          wave_inc = ( Sp%Basic%wavelength_long(band) - &
                       Sp%Basic%wavelength_short(band) ) / number_term
          DO i=1, number_term
            wavelength(1, i, band) = Sp%Basic%wavelength_short(band) + &
                                     wave_inc*(i-1)
            wavelength(2, i, band) = Sp%Basic%wavelength_short(band) + &
                                     wave_inc*(i)
          END DO
        ELSE
          WRITE(*, '(a)') 'Enter band limits (metres): '
          DO i=1, number_term
            READ(*, *, IOSTAT=ios) wavelength(1, i, band), &
                                          wavelength(2, i, band)
          END DO 
        END IF
      END IF
    END DO
    
    Sp%Var%n_sub_band = SUM(sub_bands)
    Sp%Dim%nd_sub_band = Sp%Var%n_sub_band
    
    IF (ALLOCATED(Sp%Var%index_sub_band)) &
       DEALLOCATE(Sp%Var%index_sub_band)
    IF (ALLOCATED(Sp%Var%wavelength_sub_band)) &
       DEALLOCATE(Sp%Var%wavelength_sub_band)
    ALLOCATE(Sp%Var%index_sub_band( 2, Sp%Dim%nd_sub_band ))
    ALLOCATE(Sp%Var%wavelength_sub_band( 0:2, Sp%Dim%nd_sub_band ))

    CALL map_heap_func(Sp%Basic%wavelength_short(1:Sp%Basic%n_band), &
      band_sort(1:Sp%Basic%n_band))
    sub_band=1
    DO j=1, Sp%Basic%n_band
      band=band_sort(j)
      IF (sub_bands(band) == 1) THEN
        Sp%Var%index_sub_band(1, sub_band) = band
        Sp%Var%index_sub_band(2, sub_band) = 0
        Sp%Var%wavelength_sub_band(1, sub_band) = &
          Sp%Basic%wavelength_short(band)
        Sp%Var%wavelength_sub_band(2, sub_band) = &
          Sp%Basic%wavelength_long(band)
        sub_band = sub_band + 1
      ELSE
        DO i=1, sub_bands(band)
          Sp%Var%index_sub_band(1, sub_band) = band
          index_absorb = Sp%Gas%index_absorb(1, band)
          IF (Sp%Gas%n_sub_band_gas(band, index_absorb) > 1) THEN
            Sp%Var%index_sub_band(2, sub_band) = &
              Sp%Gas%sub_band_k(i, band, index_absorb)
          ELSE
            Sp%Var%index_sub_band(2, sub_band) = i
          END IF
          Sp%Var%wavelength_sub_band(1, sub_band) = wavelength(1, i, band)
          Sp%Var%wavelength_sub_band(2, sub_band) = wavelength(2, i, band)
          sub_band = sub_band + 1
        END DO
      END IF
    END DO

    Sp%Var%n_times  = 0
    Sp%Dim%nd_times = 0
    Sp%Var%n_repeat_times  = 0
    IF (ALLOCATED(Sp%Var%time)) &
       DEALLOCATE(Sp%Var%time)
    IF (ALLOCATED(Sp%Var%total_solar_flux)) &
       DEALLOCATE(Sp%Var%total_solar_flux)
    IF (ALLOCATED(Sp%Var%solar_flux_sub_band)) &
       DEALLOCATE(Sp%Var%solar_flux_sub_band)
    IF (ALLOCATED(Sp%Var%rayleigh_coeff)) &
       DEALLOCATE(Sp%Var%rayleigh_coeff)
    ALLOCATE(Sp%Var%time( 4, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%total_solar_flux( Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%solar_flux_sub_band( Sp%Dim%nd_sub_band, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%rayleigh_coeff( Sp%Dim%nd_sub_band, 0:Sp%Dim%nd_times ))
  END IF

! Fill temporary spectral type to hold sub-bands as full bands
  SubSp%Basic%n_band = Sp%Var%n_sub_band
  SubSp%Dim%nd_band  = Sp%Var%n_sub_band
  ALLOCATE(SubSp%Basic%wavelength_short(SubSp%Basic%n_band))
  SubSp%Basic%wavelength_short = Sp%Var%wavelength_sub_band(1, :)
  ALLOCATE(SubSp%Basic%wavelength_long(SubSp%Basic%n_band))
  SubSp%Basic%wavelength_long = Sp%Var%wavelength_sub_band(2, :)
  ALLOCATE(SubSp%Basic%l_present(0:14))
  SubSp%Basic%l_present(14) = .FALSE.
  ALLOCATE(SubSp%Solar%solar_flux_band(SubSp%Basic%n_band))
  ALLOCATE(SubSp%Rayleigh%rayleigh_coeff(SubSp%Basic%n_band))

  IF ( .NOT. Sp%Basic%l_present(17) ) THEN
    ! Fill sub-band Rayleigh coefficients for mean solar spectrum
    WRITE(*, '(/a)') 'A mean solar spectrum is needed for the mean'
    WRITE(*, '(a)')  'Rayleigh coefficients per sub-band:'
    CALL make_block_3(SubSp, Sol, ierr)
    Sp%Var%rayleigh_coeff(:,0) = SubSp%Rayleigh%rayleigh_coeff
  ELSE
    ! Currently only Rayleigh coefficients for air are supported here
    SubSp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_total
  END IF

  DO
    IF (Sp%Var%n_times == 0) THEN
      ! If there are currently no times in the look-up table we can
      ! choose the number of varying Rayleigh coefficients to set 
      WRITE(*, '(a)') &
        'How many sub-bands will require a varying Rayleigh coefficient:'
      READ(*, *, IOSTAT=ios) Sp%Var%n_rayleigh_coeff
    END IF

    WRITE(*, '(a)') 'Number of times / dates to add to spectral data:'
    WRITE(*, '(a)') '  0 to finish'
    WRITE(*, '(a)') ' -1 to clear current data'
    WRITE(*, '(a)') ' -2 to add all available data times'
    READ(*, *, IOSTAT=ios) n_times
    IF (n_times == 0 .OR. ios /= 0) EXIT
    IF (n_times == -1) THEN
      Sp%Var%n_times  = 0
      Sp%Dim%nd_times = 0
      CYCLE
    END IF

    WRITE(*, '(a)') 'Enter format of data file:'
    WRITE(*, '(a)') '  5 : CMIP5 (see http://solarisheppa.geomar.de/cmip5)'
    WRITE(*, '(a)') '  6 : CMIP6 (see http://solarisheppa.geomar.de/cmip6)'
    READ(*, *, IOSTAT=ios) i_format
    SELECT CASE (i_format)
    CASE (5)
      IF (n_times == -2) THEN
        WRITE(*, '(a)') 'Number of times must be specified for CMIP5 data'
        CYCLE
      END IF
    CASE (6)
      WRITE(*, '(a)') 'Enter location of data file:'
      READ(*, *, IOSTAT=ios) cmip6_file
      ! Open the file for reading
      CALL nf(nf90_open(TRIM(cmip6_file),NF90_NOWRITE,ncid))
      ! Get number of times
      dim_name = 'time'
      CALL nf(nf90_inq_dimid(ncid, dim_name, dimid_time))
      CALL nf(nf90_inquire_dimension(ncid, dimid_time, dim_name, time_len))
      IF (n_times > 0) THEN
        WRITE(*, '(a)') 'Provide year, month, and day of first data time:'
        READ(*, *, IOSTAT=ios) yearstart, monthstart, daystart
      END IF
      IF (n_times == -2) THEN
        WRITE(*, '(a,i0)') 'Setting number of times to ', time_len
        n_times = time_len
      END IF
    CASE DEFAULT
      WRITE(*, '(a)') 'Unknown format'
      CYCLE
    END SELECT

    Sp%Dim%nd_times = Sp%Dim%nd_times + n_times

!   Save current variability data
    SpVarTmp = Sp%Var

!   Reallocate arrays to extend times
    DEALLOCATE(Sp%Var%time)
    DEALLOCATE(Sp%Var%total_solar_flux)
    DEALLOCATE(Sp%Var%solar_flux_sub_band)
    DEALLOCATE(Sp%Var%rayleigh_coeff)
    ALLOCATE(Sp%Var%time( 4, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%total_solar_flux( Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%solar_flux_sub_band( Sp%Dim%nd_sub_band, Sp%Dim%nd_times ))
    ALLOCATE(Sp%Var%rayleigh_coeff( Sp%Dim%nd_sub_band, 0:Sp%Dim%nd_times ))
    Sp%Var%time(:, 1:Sp%Var%n_times) = &
      SpVarTmp%time(:, 1:Sp%Var%n_times)
    Sp%Var%total_solar_flux(1:Sp%Var%n_times) = &
      SpVarTmp%total_solar_flux(1:Sp%Var%n_times)
    Sp%Var%solar_flux_sub_band(:, 1:Sp%Var%n_times) = &
      SpVarTmp%solar_flux_sub_band(:, 1:Sp%Var%n_times)
    Sp%Var%rayleigh_coeff(:, 0:Sp%Var%n_times) = &
      SpVarTmp%rayleigh_coeff(:, 0:Sp%Var%n_times)

    SELECT CASE (i_format)
    CASE (5)
      ! Read data from CMIP5 format spectral variability file
      CALL get_free_unit(ios, iu_solar)
      CALL open_file_in(ios, iu_solar, 'Enter location of data file:')
      
      ! Read first to find the number of points in the spectrum.
      VSol%n_points = 0
      l_count=.FALSE.
      DO
        READ(iu_solar, '(A)', IOSTAT=ios) line
        IF (ios /= 0) THEN
          EXIT
        ELSE IF (INDEX(line, 'wavelength (nm) grid') > 1) THEN
          l_count=.TRUE.
        ELSE IF (INDEX(line, 'wavelength bands') > 1) THEN
          EXIT
        ELSE IF (l_count) THEN
          VSol%n_points=VSol%n_points+1
        ENDIF
      ENDDO
      VSol%n_points = VSol%n_points * 5 ! Five wavelengths per line
      
      ALLOCATE(VSol%wavelength(VSol%n_points))
      ALLOCATE(VSol%irrad(     VSol%n_points))
      
      scale_wv=1.0E-09_RealK
      scale_irr=0.9965E+06_RealK ! Wm-3 in TIM scale
      l_monthly = .FALSE.
      REWIND(iu_solar)
      DO
        READ(iu_solar, '(a)', IOSTAT=ios) line
        IF (ios /= 0) THEN
          print*, 'end of file'
          STOP
        END IF
        IF (INDEX(line, 'wavelength (nm) grid') > 1) THEN
          DO i = 1, VSol%n_points/5
            READ(iu_solar, *, IOSTAT=ios) VSol%wavelength(i*5-4:i*5)
          END DO
          ! Read down to the start of the irradaince data
          DO
            READ(iu_solar, '(a)', IOSTAT=ios) line
            IF (INDEX(line, 'MONTH') > 1) THEN
              l_monthly = .TRUE.
            END IF
            IF (INDEX(line, 'TSI') > 1) EXIT
            IF (ios /= 0) STOP
          END DO
          EXIT
        ENDIF
      END DO
      ! Scale the values to the correct units:
      VSol%wavelength = VSol%wavelength * scale_wv
      
      ! Read in each time / date and calculate normalised spectrum
      DO i = Sp%Var%n_times + 1, Sp%Var%n_times + n_times
        IF (l_monthly) THEN
          READ(iu_solar, '(2i6,f17.6)', IOSTAT=ios) &
            Sp%Var%time(1:2, i), Sp%Var%total_solar_flux(i)
        ELSE
          READ(iu_solar, '(i6,f17.6)', IOSTAT=ios) &
            Sp%Var%time(1, i), Sp%Var%total_solar_flux(i)
          Sp%Var%time(2, i) = 1 ! Hardwire month to 1
        END IF
      
        ! Multiply by 0.9965 for TIM scale    
        Sp%Var%total_solar_flux(i) = Sp%Var%total_solar_flux(i) * 0.9965
      
        Sp%Var%time(3, i) = 1 ! Hardwire day of month to 1
        Sp%Var%time(4, i) = 0 ! Hardwire seconds since midnight to 0
      
        READ(iu_solar, *, IOSTAT=ios) VSol%irrad
        VSol%irrad = VSol%irrad * scale_irr
      
        ! Calculate the normalised solar flux in each sub-band
        CALL make_block_2_1(SubSp, VSol, filter,.FALSE.,.TRUE.,.FALSE.,ierr)
        Sp%Var%solar_flux_sub_band(:,i) = SubSp%Solar%solar_flux_band
      
        ! Calculate Rayleigh scattering coefficients in each sub-band
        CALL make_block_3_1(SubSp, VSol, Refract, .FALSE.)
        Sp%Var%rayleigh_coeff(:,i) = SubSp%Rayleigh%rayleigh_coeff
      END DO
      Sp%Var%n_times = Sp%Var%n_times + n_times
      
      DEALLOCATE(VSol%irrad)
      DEALLOCATE(VSol%wavelength)
      CLOSE(iu_solar)

    CASE (6)
      VSol%l_binned = .TRUE.
      ! Find the number of points in the spectrum.
      dim_name = 'wlen'
      CALL nf(nf90_inq_dimid(ncid, dim_name, dimid_wlen))
      CALL nf(nf90_inquire_dimension(ncid, dimid_wlen, dim_name, wlen_len))
      VSol%n_points = wlen_len

      ALLOCATE(VSol%wavelength(VSol%n_points))
      ALLOCATE(VSol%irrad(     VSol%n_points))
      ALLOCATE(VSol%bandsize(  VSol%n_points))
      ALLOCATE(VSol%bandbnds(  2, VSol%n_points))

      scale_wv=1.0E-09_RealK ! nm to m
      scale_irr=1.0E+09_RealK ! Wm-2nm-1 to Wm-3

      ! Read the wavelength bin centres
      CALL nf(nf90_inq_varid(ncid, 'wlen', varid))
      CALL nf(nf90_get_var(ncid, varid, VSol%wavelength))
      ! Scale the values to the correct units:
      VSol%wavelength = VSol%wavelength * scale_wv

      ! Read the start time for each time bin
      CALL nf(nf90_get_att(ncid, NF90_GLOBAL, 'frequency', frequency))
      ALLOCATE(calyear(time_len))
      ios = nf90_inq_varid(ncid, 'calyear', varid)
      IF (ios == NF90_NOERR) ios = nf90_get_var(ncid, varid, calyear)
      IF (ios /= NF90_NOERR) calyear(:) = 1850
      ALLOCATE(calmonth(time_len))
      ios = nf90_inq_varid(ncid, 'calmonth', varid)
      IF (ios == NF90_NOERR) ios = nf90_get_var(ncid, varid, calmonth)
      IF (ios /= NF90_NOERR) calmonth(:) = 1
      ALLOCATE(calday(time_len))
      ios = nf90_inq_varid(ncid, 'calday', varid)
      IF (ios == NF90_NOERR) ios = nf90_get_var(ncid, varid, calday)
      IF (ios /= NF90_NOERR) calday(:) = 1
      ALLOCATE(seconds(time_len))
      ios = nf90_inq_varid(ncid, 'seconds', varid)
      IF (ios == NF90_NOERR) ios = nf90_get_var(ncid, varid, seconds)
      IF (ios /= NF90_NOERR) seconds(:) = 0 ! Defaults to 0 if it fails
      
      ! Read the tsi
      ALLOCATE(tsi(time_len))
      CALL nf(nf90_inq_varid(ncid, 'tsi', varid))
      CALL nf(nf90_get_var(ncid, varid, tsi))

      ! Read the wbinsize
      ALLOCATE(wbinsize(wlen_len))
      CALL nf(nf90_inq_varid(ncid, 'wlenbinsize', varid))
      CALL nf(nf90_get_var(ncid, varid, wbinsize))
      VSol%bandsize = wbinsize * scale_wv

      ! Read the wbinbnds
      ALLOCATE(wbinbnds(2,wlen_len))
      CALL nf(nf90_inq_varid(ncid, 'wlen_bnds', varid))
      CALL nf(nf90_get_var(ncid, varid, wbinbnds))
      VSol%bandbnds = wbinbnds * scale_wv
      

      ! Find the ssi variable id
      ALLOCATE(ssi(1, wlen_len))
      CALL nf(nf90_inq_varid(ncid, 'ssi', varid))

      IF (yearstart == imdi) THEN
        jj = 1
      ELSE
        DO jj = 1, time_len
          IF (calyear(jj) >  yearstart) EXIT
          IF (calyear(jj) == yearstart) THEN
            IF (calmonth(jj) >  monthstart) EXIT
            IF (calmonth(jj) == monthstart) THEN
              IF (calday(jj) >= daystart) EXIT
            END IF
          END IF
        END DO
      END IF
      DO k = 1, n_times
         i = k + Sp%Var%n_times
         j = jj + k - 1
        SELECT CASE (TRIM(frequency))
        CASE DEFAULT
          Sp%Var%time(1, i) = calyear(j)
          Sp%Var%time(2, i) = calmonth(j)
          Sp%Var%time(3, i) = calday(j)
          Sp%Var%time(4, i) = 0 ! Hardwire seconds since midnight to 0
        CASE ("mon")
          Sp%Var%time(1, i) = calyear(j)
          Sp%Var%time(2, i) = calmonth(j)
          Sp%Var%time(3, i) = 1 ! Hardwire day of month to 1
          Sp%Var%time(4, i) = 0 ! Hardwire seconds since midnight to 0
        CASE ("fx")
          Sp%Var%time(1, i) = 1850
          Sp%Var%time(2, i) = 1
          Sp%Var%time(3, i) = 1
          Sp%Var%time(4, i) = 0
        CASE ("sec")
          Sp%Var%time(1, i) = calyear(j)
          Sp%Var%time(2, i) = calmonth(j)
          Sp%Var%time(3, i) = calday(j)
          Sp%Var%time(4, i) = seconds(j)
        END SELECT

        Sp%Var%total_solar_flux(i) = tsi(j)

        ! Read the ssi for this time
        CALL nf(nf90_get_var(ncid, varid, ssi, &
          start=(/1, j/), count=(/wlen_len, 1/) ))
        ! Scale the values to the correct units:
        VSol%irrad = ssi(1,:) * scale_irr

        ! Uncomment to output a mean solar spectrum file
        ! IF (k==1) THEN
        !   VSolMean%n_points    = VSol%n_points
        !   VSolMean%t_effective = VSol%t_effective
        !   VSolMean%radius      = VSol%radius
        !   VSolMean%l_binned    = VSol%l_binned
        !   ALLOCATE(VSolMean%wavelength( VSolMean%n_points))
        !   ALLOCATE(VSolMean%irrad(      VSolMean%n_points))
        !   ALLOCATE(VSolMean%bandsize(   VSolMean%n_points))
        !   ALLOCATE(VSolMean%bandbnds(2, VSolMean%n_points))
        !   VSolMean%wavelength = VSol%wavelength
        !   VSolMean%irrad = 0.0_RealK
        !   VSolMean%bandsize = VSol%bandsize
        !   VSolMean%bandbnds = VSol%bandbnds
        ! END IF
        ! VSolMean%irrad = VSolMean%irrad + VSol%irrad
        ! IF (k==n_times) THEN
        !   VSolMean%irrad = VSolMean%irrad/REAL(n_times, RealK)
        !   CALL write_solar_spectrum('cmip6_solar_spectrum', VSolMean, ierr)
        !   DEALLOCATE(VSolMean%bandbnds)
        !   DEALLOCATE(VSolMean%bandsize)
        !   DEALLOCATE(VSolMean%irrad)
        !   DEALLOCATE(VSolMean%wavelength)
        ! END IF

        ! Calculate the normalised solar flux in each sub-band
        CALL make_block_2_1(SubSp, VSol, filter,.FALSE.,.TRUE.,.FALSE.,ierr)
        Sp%Var%solar_flux_sub_band(:,i) = SubSp%Solar%solar_flux_band
      
        ! Calculate Rayleigh scattering coefficients in each sub-band
        CALL make_block_3_1(SubSp, VSol, Refract, .FALSE.)
        Sp%Var%rayleigh_coeff(:,i) = SubSp%Rayleigh%rayleigh_coeff
      END DO
      Sp%Var%n_times = Sp%Var%n_times + n_times
      
      CALL nf(nf90_close(ncid))
      DEALLOCATE(ssi)
      DEALLOCATE(wbinbnds)
      DEALLOCATE(wbinsize)
      DEALLOCATE(tsi)
      DEALLOCATE(seconds)
      DEALLOCATE(calday)
      DEALLOCATE(calmonth)
      DEALLOCATE(calyear)
      DEALLOCATE(VSol%bandbnds)
      DEALLOCATE(VSol%bandsize)
      DEALLOCATE(VSol%irrad)
      DEALLOCATE(VSol%wavelength)

    END SELECT
  END DO

  WRITE(*, '(a)') 'How many of the final times / dates should be '
  WRITE(*, '(a)') 'periodically repeated into the future:'
  READ(*, *, IOSTAT=ios) Sp%Var%n_repeat_times

  Sp%Basic%l_present(17)=.TRUE.

CONTAINS

  SUBROUTINE nf(status)
    USE netcdf
    INTEGER, INTENT(IN):: status
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       STOP 'STOPPED!'
    END IF
  END SUBROUTINE nf

END SUBROUTINE make_block_17

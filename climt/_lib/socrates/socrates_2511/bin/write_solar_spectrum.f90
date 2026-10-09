! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to write a solar spectrum file.
!
SUBROUTINE write_solar_spectrum(filename, SolarSpec, ierr) 

  USE def_std_io_icf, ONLY: iu_err
  USE def_solarspec, ONLY: StrSolarSpec
  USE rad_pcf, ONLY: i_err_fatal

  IMPLICIT NONE


  CHARACTER (LEN=*), Intent(IN) :: filename
!   Name of solar spectrum file
  TYPE (StrSolarSpec), Intent(IN) :: SolarSpec
!   Solar spectrum
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables.
  INTEGER :: iu_solar
!   Unit number for the solar spectrum
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i
!   Loop variable


! Get a free unit number for the solar spectrum.
  CALL get_free_unit(ierr, iu_solar)

! Open the file for writing.
  OPEN(UNIT=iu_solar, FILE=filename, IOSTAT=ios, STATUS='UNKNOWN')
  IF (ios /= 0) THEN
    WRITE(iu_err, '(/a)') &
      '***Error: Cannot open file to write solar spectrum'
    ierr=i_err_fatal
    RETURN
  ENDIF

  WRITE(iu_solar, '(A)', IOSTAT=ios) '*RADIUS'
  WRITE(iu_solar, '(1pe16.9)', IOSTAT=ios) SolarSpec%radius
  WRITE(iu_solar, '(A)', IOSTAT=ios) '*TEMPERATURE'
  WRITE(iu_solar, '(1pe16.9)', IOSTAT=ios) SolarSpec%t_effective
  IF (SolarSpec%l_binned) THEN
    WRITE(iu_solar, '(A)', IOSTAT=ios) &
      '      Wavelength Bounds           Irradiance'
    WRITE(iu_solar, '(A)', IOSTAT=ios) &
      '             (m)                    (W/m3)'
    WRITE(iu_solar, '(A)', IOSTAT=ios) '*BEGIN_BINNED_DATA'
    DO i = 1, SolarSpec%n_points
      WRITE(iu_solar, '(3(1pe16.9))', IOSTAT=ios) &
        SolarSpec%bandbnds(:,i), &
        SolarSpec%irrad(i)
    ENDDO
    WRITE(iu_solar, '(A)', IOSTAT=ios) '*END'
  ELSE
    WRITE(iu_solar, '(A)', IOSTAT=ios) &
      '   Wavelength      Irradiance'
    WRITE(iu_solar, '(A)', IOSTAT=ios) &
      '      (m)            (W/m3)'
    WRITE(iu_solar, '(A)', IOSTAT=ios) '*BEGIN_DATA'
    DO i = 1, SolarSpec%n_points
      WRITE(iu_solar, '(2(1pe16.9))', IOSTAT=ios) &
        SolarSpec%wavelength(i), &
        SolarSpec%irrad(i)
    ENDDO
    WRITE(iu_solar, '(A)', IOSTAT=ios) '*END'
  END IF

  CLOSE(iu_solar)

END SUBROUTINE write_solar_spectrum

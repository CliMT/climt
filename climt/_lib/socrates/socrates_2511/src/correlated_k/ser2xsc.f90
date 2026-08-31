! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert Serdyuchenko ozone data to HITRAN xsc format.
PROGRAM ser2xsc

  USE realtype_rd
  USE def_hitran_record

  IMPLICIT NONE

  INTEGER :: ierr
!   Error flag
  INTEGER :: iu_dat, iu_xsc
!   File unit numbers
  INTEGER :: i, j, d1, d2
  INTEGER, PARAMETER :: header_length = 45
  INTEGER, PARAMETER :: record_length = 12
  INTEGER, PARAMETER :: data_length = 88668
  TYPE(StrXscHead) :: xsc
  REAL (RealK), PARAMETER :: xsc_res = 0.5_RealK

  REAL (RealK), PARAMETER :: temperature(record_length) = &
    (/ 0.0, 293.0, 283.0, 273.0, 263.0, 253.0,            &
     243.0, 233.0, 223.0, 213.0, 203.0, 193.0 /)

  REAL (RealK) :: ser_data(record_length, data_length)
  REAL (RealK) :: ser_wn(data_length), ser_wn_diff(data_length)
  REAL (RealK), ALLOCATABLE :: xsc_data(:), xsc_wn(:)

! Open the file 'serdyuchenkogorshelev5digits.dat'
  CALL get_free_unit(ierr, iu_dat)
  OPEN(UNIT=iu_dat, FILE='serdyuchenkogorshelev5digits.dat', &
       IOSTAT=ierr, STATUS='OLD')

! Skip over the header
  DO i = 1, header_length
    READ(iu_dat, *)
  END DO
! Read the cross-section data
  READ(iu_dat, *) ser_data
  CLOSE(iu_dat)

  ser_data = ABS(ser_data)
! Convert wavelength in nm to wavenumber in cm-1
  ser_wn = 1.0e7_RealK/ser_data(1, :)

! Construct common XSC header information
  xsc%chemical_symbol='                  O3'
!  xsc%wavenumber_min = REAL(FLOOR(ser_wn(data_length)), RealK)
  xsc%wavenumber_min = 9091.0_RealK ! 1099.99 nm
  xsc%wavenumber_max = REAL(CEILING(ser_wn(1)), RealK)
  xsc%no_pts=(xsc%wavenumber_max-xsc%wavenumber_min)/xsc_res + 1
  xsc%pressure = 0.0
  WRITE(xsc%resolution,'(f5.1)') xsc_res
  xsc%common_name = '          ozone'
  xsc%broadening_gas = ''
  xsc%re_no = 19

  ALLOCATE(xsc_data(xsc%no_pts))
  ALLOCATE(xsc_wn(xsc%no_pts))

! Construct a regular wavenumber grid to cover the range of wavelengths
  xsc_wn(1) = xsc%wavenumber_min
  DO j = 2, xsc%no_pts
    xsc_wn(j) = xsc_wn(j-1) + xsc_res
  END DO

! Open a file for the output HITRAN xsc database
  CALL get_free_unit(ierr, iu_xsc)
  OPEN(UNIT=iu_xsc, FILE='serdyuchenko_o3.xsc', &
       IOSTAT=ierr, STATUS='UNKNOWN')

! Loop over temperatures
  DO i = 2, record_length
!   Interpolate wavelength values onto the wavenumber grid
    DO j = 1, xsc%no_pts
      ser_wn_diff = ser_wn-xsc_wn(j)
      d1 = MINLOC(ser_wn_diff, 1, ser_wn_diff > 0.0_RealK)
      d2 = d1 + 1 
      xsc_data(j) = (ser_data(i, d2)*ser_wn_diff(d1)  &
                    -ser_data(i, d1)*ser_wn_diff(d2)) &
                  / (ser_wn_diff(d1) - ser_wn_diff(d2))
    END DO
    xsc%temperature = temperature(i)
    xsc%max_xsc = MAXVAL(xsc_data)
!   Write out the data in HITRAN xsc format
    WRITE(iu_xsc, xsc_header_frmt) xsc
    WRITE(iu_xsc, xsc_data_frmt) xsc_data
  END DO
  CLOSE(iu_xsc)

  DEALLOCATE(xsc_wn)
  DEALLOCATE(xsc_data)

END PROGRAM ser2xsc

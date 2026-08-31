! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert JPL-2010 O2 cross-section data to HITRAN xsc format.
PROGRAM jpl2xsc

  USE realtype_rd
  USE def_hitran_record

  IMPLICIT NONE

  INTEGER :: ierr
!   Error flag
  INTEGER :: iu_dat, iu_xsc
!   File unit numbers
  INTEGER :: i, j, d1, d2
  INTEGER, PARAMETER :: record_length = 1
  INTEGER, PARAMETER :: data_length = 41
  INTEGER, PARAMETER :: header_length(record_length) = &
    (/ 24 /)
  TYPE(StrXscHead) :: xsc
  REAL (RealK), PARAMETER :: xsc_res = 100.0_RealK

  REAL (RealK), PARAMETER :: temperature(record_length) = &
    (/ 298.0 /)
  CHARACTER (LEN=41) :: jplfile(record_length) = &
    (/ 'O2_JPL-2010(2011)_298K_205-245nm(rec).txt' /)

  REAL (RealK) :: jpl_dat(2, data_length)
  REAL (RealK) :: jpl_data(record_length, data_length)
  REAL (RealK) :: jpl_wn(data_length), jpl_wn_diff(data_length)
  REAL (RealK), ALLOCATABLE :: xsc_data(:), xsc_wn(:)

  DO i = 1, record_length
!   Open the cross-section data file
    CALL get_free_unit(ierr, iu_dat)
    OPEN(UNIT=iu_dat, FILE=jplfile(i), IOSTAT=ierr, STATUS='OLD')

!   Skip over the header
    DO j = 1, header_length(i)
      READ(iu_dat, *)
    END DO
!   Read the cross-section data
    READ(iu_dat, *) jpl_dat
    CLOSE(iu_dat)

    jpl_data(i,:) = ABS(jpl_dat(2,:))
    IF (i==1) THEN
      jpl_wn = jpl_dat(1,:)
    ELSE
      IF (ANY(jpl_wn /= jpl_dat(1,:))) THEN
        PRINT*, "Wavelengths in files don't match"
        DO j = 1, data_length
          IF (jpl_wn(j) /= jpl_dat(1,j)) THEN
            PRINT*, 'Record ',jplfile(1),' ',jplfile(i)
            PRINT*, i, jpl_wn(j), jpl_dat(1,j)
          END IF
        END DO
        STOP
      END IF
    END IF
  END DO

! Convert wavelength in nm to wavenumber in cm-1
  jpl_wn = 1.0e7_RealK/jpl_wn

! Construct common XSC header information
  xsc%chemical_symbol='                  O2'
!  xsc%wavenumber_min = REAL(FLOOR(jpl_wn(data_length)), RealK)
!  xsc%wavenumber_max = REAL(CEILING(jpl_wn(1)), RealK)
  xsc%wavenumber_min = 40880.0_RealK ! 244.62 nm
  xsc%wavenumber_max = 48780.0_RealK ! 205.00 nm
  xsc%no_pts=(xsc%wavenumber_max-xsc%wavenumber_min)/xsc_res + 1
  xsc%pressure = 0.0
  WRITE(xsc%resolution,'(f5.1)') xsc_res
  xsc%common_name = '         oxygen'
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
  OPEN(UNIT=iu_xsc, FILE='jpl_o2.xsc', &
       IOSTAT=ierr, STATUS='UNKNOWN')

! Loop over temperatures
  DO i = 1, record_length
!   Interpolate wavelength values onto the wavenumber grid
    DO j = 1, xsc%no_pts
      jpl_wn_diff = jpl_wn-xsc_wn(j)
      d1 = MINLOC(jpl_wn_diff, 1, jpl_wn_diff > 0.0_RealK)
      d2 = d1 + 1 
      xsc_data(j) = (jpl_data(i, d2)*jpl_wn_diff(d1)  &
                    -jpl_data(i, d1)*jpl_wn_diff(d2)) &
                  / (jpl_wn_diff(d1) - jpl_wn_diff(d2))
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

END PROGRAM jpl2xsc

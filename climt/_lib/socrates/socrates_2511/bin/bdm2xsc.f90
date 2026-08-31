! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert Brion-Daumont-Malicet ozone data to HITRAN xsc format.
PROGRAM bdm2xsc

  USE realtype_rd
  USE def_hitran_record

  IMPLICIT NONE

  INTEGER :: ierr
!   Error flag
  INTEGER :: iu_dat, iu_xsc
!   File unit numbers
  INTEGER :: i, j, d1, d2
  INTEGER, PARAMETER :: record_length = 4
  INTEGER, PARAMETER :: data_length = 32402 ! 519.01nm
  INTEGER, PARAMETER :: header_length(record_length) = &
    (/ 51, 51, 51, 1 /)
  TYPE(StrXscHead) :: xsc
  REAL (RealK), PARAMETER :: xsc_res = 1.0_RealK

  REAL (RealK), PARAMETER :: temperature(record_length) = &
    (/ 218.0, 228.0, 243.0, 295.0 /)
  CHARACTER (LEN=19) :: bdmfile(record_length) = &
    (/ 'O3_CRS_BDM_218K.dat', 'O3_CRS_BDM_228K.dat', &
       'O3_CRS_BDM_243K.dat', 'O3_CRS_BDM_295K.dat' /)

  REAL (RealK) :: bdm_dat(2, data_length)
  REAL (RealK) :: bdm_data(record_length, data_length)
  REAL (RealK) :: bdm_wn(data_length), bdm_wn_diff(data_length)
  REAL (RealK), ALLOCATABLE :: xsc_data(:), xsc_wn(:)

  DO i = 1, record_length
!   Open the cross-section data file
    CALL get_free_unit(ierr, iu_dat)
    OPEN(UNIT=iu_dat, FILE=bdmfile(i), IOSTAT=ierr, STATUS='OLD')

!   Skip over the header
    DO j = 1, header_length(i)
      READ(iu_dat, *)
    END DO
!   Read the cross-section data
    READ(iu_dat, *) bdm_dat
    CLOSE(iu_dat)

    bdm_data(i,:) = ABS(bdm_dat(2,:))
    IF (i==1) THEN
      bdm_wn = bdm_dat(1,:)
    ELSE
      IF (ANY(bdm_wn /= bdm_dat(1,:))) THEN
        PRINT*, "Wavelengths in files don't match"
        DO j = 1, data_length
          IF (bdm_wn(j) /= bdm_dat(1,j)) THEN
            PRINT*, 'Record ',bdmfile(1),' ',bdmfile(i)
            PRINT*, i, bdm_wn(j), bdm_dat(1,j)
          END IF
        END DO
        STOP
      END IF
    END IF
  END DO

! Convert wavelength in Angstrom to wavenumber in cm-1
  bdm_wn = 1.0e8_RealK/bdm_wn

! Construct common XSC header information
  xsc%chemical_symbol='                  O3'
!  xsc%wavenumber_min = REAL(FLOOR(bdm_wn(data_length)), RealK)
  xsc%wavenumber_min = 46875.0_RealK ! 213.33 nm
  xsc%wavenumber_max = REAL(CEILING(bdm_wn(1)), RealK)
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
  OPEN(UNIT=iu_xsc, FILE='bdm_o3.xsc', &
       IOSTAT=ierr, STATUS='UNKNOWN')

! Loop over temperatures
  DO i = 1, record_length
!   Interpolate wavelength values onto the wavenumber grid
    DO j = 1, xsc%no_pts
      bdm_wn_diff = bdm_wn-xsc_wn(j)
      d1 = MINLOC(bdm_wn_diff, 1, bdm_wn_diff > 0.0_RealK)
      d2 = d1 + 1 
      xsc_data(j) = (bdm_data(i, d2)*bdm_wn_diff(d1)  &
                    -bdm_data(i, d1)*bdm_wn_diff(d2)) &
                  / (bdm_wn_diff(d1) - bdm_wn_diff(d2))
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

END PROGRAM bdm2xsc

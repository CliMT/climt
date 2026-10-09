! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to convert cross-section data to HITRAN xsc format.
PROGRAM dat2xsc

  USE realtype_rd, ONLY: RealK
  USE def_std_io_icf, ONLY: iu_err
  USE def_hitran_record, ONLY: StrXscHead, &
    xsc_header_format, xsc_data_format, uvxsc_header_frmt, uvxsc_data_frmt

  IMPLICIT NONE

  INTEGER :: ierr, ios
!   Error flag
  INTEGER :: iu_dat, iu_xsc
!   File unit numbers
  INTEGER :: i, jj, k, l
  INTEGER :: data_length
  TYPE(StrXscHead) :: xsc
  INTEGER :: xsc_res

  CHARACTER (LEN=256) :: infile, outfile, arg
  LOGICAL :: l_infile, l_outfile, l_offset
  
  REAL (RealK), ALLOCATABLE :: in_dat(:, :)
  REAL (RealK), ALLOCATABLE :: in_data(:)
  REAL (RealK), ALLOCATABLE :: in_wn(:)
  REAL (RealK), ALLOCATABLE :: out_wn(:, :)
  REAL (RealK), ALLOCATABLE :: xsc_data(:), xsc_wn(:)
  REAL (RealK) :: temperature = 298.0

  l_infile = .FALSE.
  l_outfile = .FALSE.
  l_offset = .FALSE.
  i = 0
  DO
    i=i+1
    IF (i > command_argument_count()) EXIT
    CALL get_command_argument(i, arg)
    SELECT CASE (arg)
    CASE ('-t','--temperature')
      i=i+1
      CALL get_command_argument(i, arg)
      READ(arg, *) temperature
    CASE ('--offset','--eric')
      l_offset = .TRUE.
    CASE default
      IF (l_infile) THEN
        ! If the input file has been provided, the second file is
        ! the output .(uv)xsc file
        outfile = arg
        l_outfile = .TRUE.
      ELSE
        ! First file given is the input .dat file
        infile = arg
        l_infile = .TRUE.
      END IF
    END SELECT
  END DO
  IF (.NOT.l_outfile) THEN
    WRITE(iu_err, '(a)') &
      'Usage: dat2xsc [-t <temperature>] [--offset] infile.dat outfile.[uv]xsc'
    STOP
  END IF

! Open the cross-section data file
  CALL get_free_unit(ierr, iu_dat)
  IF (ierr > 0) THEN
    WRITE(iu_err, '(A, i5)') 'Error in get_free_unit: ', ierr
    STOP
  END IF
  OPEN(UNIT=iu_dat, FILE=infile, IOSTAT=ierr, STATUS='OLD')
  IF (ierr > 0) THEN
    WRITE(iu_err, '(A, i5, A, A)') 'Error', ierr, ' opening file ', infile
    STOP
  END IF

! Count the number of lines in the file
  data_length = 0
  DO
    READ(iu_dat, *, iostat=ios)
    IF (ios /= 0) EXIT
    data_length = data_length + 1
  END DO

  REWIND(iu_dat)

  ALLOCATE ( in_wn(data_length)     )
  ALLOCATE ( in_data(data_length)   )
  ALLOCATE ( in_dat(2, data_length) )

! Read the cross-section data
  READ(iu_dat, *) in_dat
  CLOSE(iu_dat)

  in_wn = in_dat(1,:)
  in_data = ABS(in_dat(2,:))

  DEALLOCATE ( in_dat )
  ALLOCATE ( out_wn(2, data_length) )

! Construct the wavelngth bin limits to be used for each XSC record
  IF (l_offset) THEN
    ! The wavelengths in the .dat file designate the start of the bin
    DO jj = 1, data_length-1
      out_wn(1,jj) = in_wn(jj+1)
      out_wn(2,jj) = in_wn(jj)
    END DO
    jj = data_length
    out_wn(1,jj) = 2.0*in_wn(jj) - in_wn(jj-1)
    out_wn(2,jj) = in_wn(jj)
  ELSE
    ! The wavelengths in the .dat file designate point values
    ! or the centre of the bin
    jj = 1
    out_wn(1,jj) = (in_wn(jj) + in_wn(jj+1))/2.0
    out_wn(2,jj) = 2.0*in_wn(jj) - out_wn(1,jj)
    DO jj = 2, data_length-1
      out_wn(1,jj) = (in_wn(jj) + in_wn(jj+1))/2.0
      out_wn(2,jj) = (in_wn(jj-1) + in_wn(jj))/2.0
    END DO
    jj = data_length
    out_wn(2,jj) = (in_wn(jj-1) + in_wn(jj))/2.0
    out_wn(1,jj) = 2.0*in_wn(jj) - out_wn(2,jj)
  END IF

! Convert wavelength in nm to wavenumber in cm-1
  out_wn = 1.0e7_RealK/out_wn

  ! Construct common XSC header information
  xsc%chemical_symbol = ''
  xsc%temperature = temperature
  xsc%common_name = ''
  xsc%no_pts = 1
  xsc%pressure = 0.0
  xsc%broadening_gas = ''
  xsc%re_no = 0

! Open a file for the output HITRAN xsc database
  CALL get_free_unit(ierr, iu_xsc)
  IF (ierr > 0) THEN
    WRITE(iu_err, '(A, i5)') 'Error in get_free_unit: ', ierr
    STOP
  END IF
  OPEN(UNIT=iu_xsc, FILE=outfile, IOSTAT=ierr, STATUS='UNKNOWN')
  IF (ierr > 0) THEN
    WRITE(iu_err, '(A, i5, A, A)') 'Error', ierr, ' opening file ', outfile
    STOP
  END IF

! Select the cross-section file format based on the filename extension
! The default is the HITRAN xsc format.
  l = LEN_TRIM(outfile)
  IF (outfile(MAX(l-5,1):l) == '.uvxsc') THEN
    ! A bespoke format is used for UV cross-section data
    xsc_header_format = uvxsc_header_frmt
    xsc_data_format = uvxsc_data_frmt
  END IF

  DO jj = 1, data_length
    k = data_length + 1 - jj
    xsc%wavenumber_min = out_wn(1,k)
    xsc%wavenumber_max = out_wn(2,k)
!   Resolution in mAngstrom
    xsc_res = NINT(1.0e11_RealK/xsc%wavenumber_min &
                 - 1.0e11_RealK/xsc%wavenumber_max)
    ALLOCATE(xsc_data(xsc%no_pts))
    xsc_data(:) = in_data(k)
    xsc%max_xsc = MAXVAL(xsc_data)
    WRITE(xsc%resolution,'(i5)') xsc_res
!   Write out the data in HITRAN xsc format
    WRITE(iu_xsc, xsc_header_format) xsc
    WRITE(iu_xsc, xsc_data_format) xsc_data
    DEALLOCATE(xsc_data)
  END DO

  CLOSE(iu_xsc)

  DEALLOCATE ( out_wn  )
  DEALLOCATE ( in_data )
  DEALLOCATE ( in_wn   )

END PROGRAM dat2xsc

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to convert between cross-section data formats.
program xsc2xsc

  use realtype_rd, only: realk
  use def_std_io_icf, only: iu_err
  use def_hitran_record, only: StrXscHead, &
    xsc_header_frmt, xsc_data_frmt, uvxsc_header_frmt, uvxsc_data_frmt

  implicit none

  integer :: l
  integer :: ierr, iu_input, iu_output
  character (len=256) :: infile
  character (len=256) :: outfile
  character (len=256) :: input_header_format, input_data_format
  character (len=256) :: output_header_format, output_data_format
  type(StrXscHead) :: xsc_header
  real (RealK), allocatable :: xsc_data(:)

  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)

! Select the cross-section file format based on the filename extension
  l = len_trim(infile)
  IF (infile(MAX(l-5,1):l) == '.uvxsc') THEN
    ! A bespoke format is used for UV cross-section data
    input_header_format = uvxsc_header_frmt
    input_data_format = uvxsc_data_frmt
  ELSE
    ! The default is the HITRAN xsc format
    input_header_format = xsc_header_frmt
    input_data_format = xsc_data_frmt
  END IF
  l = len_trim(outfile)
  IF (outfile(MAX(l-5,1):l) == '.uvxsc') THEN
    output_header_format = uvxsc_header_frmt
    output_data_format = uvxsc_data_frmt
  ELSE
    output_header_format = xsc_header_frmt
    output_data_format = xsc_data_frmt
  END IF

! Open the files
  call get_free_unit(ierr, iu_input)
  if (ierr > 0) then
    write(iu_err, '(a, i5)') 'Error in get_free_unit: ', ierr
    stop
  end if
  open(UNIT=iu_input, FILE=infile, IOSTAT=ierr, STATUS='OLD')
  if (ierr > 0) then
    write(iu_err, '(a, i5, a, a)') 'Error', ierr, ' opening file ', infile
    stop
  end if
  call get_free_unit(ierr, iu_output)
  if (ierr > 0) then
    write(iu_err, '(a, i5)') 'Error in get_free_unit: ', ierr
    stop
  end if
  open(UNIT=iu_output, FILE=outfile, IOSTAT=ierr, STATUS='UNKNOWN')
  if (ierr > 0) then
    write(iu_err, '(a, i5, a, a)') 'Error', ierr, ' opening file ', outfile
    stop
  end if

! Read the input data and write the output
  do
    read(iu_input, input_header_format, IOSTAT=ierr) xsc_header
    if (ierr < 0) exit ! EOF
    allocate(xsc_data(xsc_header%no_pts))
    read(iu_input, input_data_format, IOSTAT=ierr) xsc_data
    write(iu_output, output_header_format) xsc_header
    write(iu_output, output_data_format) xsc_data
    deallocate(xsc_data)
  end do

  close(iu_input)
  close(iu_output)

end program xsc2xsc

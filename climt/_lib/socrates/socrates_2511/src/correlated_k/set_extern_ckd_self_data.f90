! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine to read the external data for the CKD self continuum.
SUBROUTINE set_extern_ckd_self_data(ierr)

! Method:
!   The name of the file is requested. The file is opened
!   and read. Directives encoutered are processed as the
!   data are read.


! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE ckd_extern_data
  USE error_pcf
  USE netcdf

  IMPLICIT NONE


  INTEGER, Intent(INOUT) :: ierr
!   Error flag


! Local variables.
  CHARACTER (LEN=256) :: line
!   Line read from input file
  CHARACTER (LEN=256) :: file_in
!   Name of input file
  CHARACTER(LEN=11) :: dim_name='wavenumbers'
!   Name of wavenumber dimension
  LOGICAL :: l_exist
!   Existence flag for file
  INTEGER :: iu_ckd_data, ncid
!   Unit number for reading CKD data
  INTEGER :: ios
!   I/O error flag
  INTEGER :: k
!   Loop variable
  INTEGER :: dimid, varid
!   Dimension and variable ID
  INTEGER :: status
!   Output status of netcdf call


  WRITE(*, '(/a)') &
    'Enter the name of the file containing the self-broadened CKD data at 296K.'
  READ(iu_stdin, '(a)', iostat=ios) file_in
  IF (ios /= 0) THEN
    WRITE(iu_err, "(/a)") '*** Error reading self-broadened continuum file name'
    ierr = i_err_fatal
    RETURN
  END IF
  INQUIRE(FILE=TRIM(file_in), EXIST=l_exist)
  IF (.NOT.l_exist) THEN
    WRITE(iu_err, '(a)') '*** The self-broadened continuum file does not exist'
    ierr=i_err_exist
    RETURN
  END IF

  IF (INDEX(file_in, '.nc') > 0) THEN

    ! Open the file for reading
    CALL nf(nf90_open(TRIM(file_in), NF90_NOWRITE, ncid))

    ! Get length of wavenumber dimension
    CALL nf(nf90_inq_dimid(ncid, dim_name, dimid))
    CALL nf(nf90_inquire_dimension(ncid, dimid, dim_name, &
            c_self_h2o_296%n_freq))
    ALLOCATE(c_self_h2o_296%c(c_self_h2o_296%n_freq))

    ! Get wavenumber data assuming a regular grid
    CALL nf(nf90_inq_varid(ncid, dim_name, varid))
    CALL nf(nf90_get_var(ncid, varid, c_self_h2o_296%c))
    c_self_h2o_296%table_start = c_self_h2o_296%c(1)
    c_self_h2o_296%table_end = c_self_h2o_296%c(c_self_h2o_296%n_freq)
    c_self_h2o_296%table_inc = ( c_self_h2o_296%table_end &
                               - c_self_h2o_296%table_start ) &
                             / REAL( c_self_h2o_296%n_freq - 1, RealK )

    ! Get self-broadened continuum coefficients
    CALL nf(nf90_inq_varid(ncid, 'self_absco_ref', varid))
    CALL nf(nf90_get_var(ncid, varid, c_self_h2o_296%c))

    ! Scale coefficients by 1e20 as used in table format
    c_self_h2o_296%c = c_self_h2o_296%c * 1.0e20_RealK

    ! Uncomment to write out data in old file format
    ! WRITE(296, "(a, a)") 'Self continuum data at 296K from ', TRIM(file_in)
    ! WRITE(296, "(a/)") '*BEGIN_DATA'
    ! WRITE(296, "(a, 1pe12.5, a)") '     Start of table                = ', &
    !   c_self_h2o_296%table_start, ' cm-1'
    ! WRITE(296, "(a, 1pe12.5, a)") '     End of table                  = ', &
    !   c_self_h2o_296%table_end, ' cm-1'
    ! WRITE(296, "(a, 1pe12.5, a)") '     Increment in table            = ', &
    !   c_self_h2o_296%table_inc, ' cm-1'
    ! WRITE(296, "(a, i5)") '     Number of points in table     = ', &
    !   c_self_h2o_296%n_freq
    ! WRITE(296, "(/a)") '     Continuum coefficients (cm^3.molec-1 * 1.0E-20)'
    ! WRITE(296, "(4(1pe15.5))") c_self_h2o_296%c
    ! WRITE(296, "(a/)") '*END'

    CALL nf(nf90_close(ncid))

  ELSE

    CALL get_free_unit(ierr, iu_ckd_data)
    OPEN(UNIT=iu_ckd_data, FILE=TRIM(file_in), iostat=ios, status='old')
    DO
      READ(iu_ckd_data, "(a)", IOSTAT=ios) line
      IF (ios /= 0) THEN
        WRITE(iu_err, "(/a)") "*** Error reading self-broadened continuum."
        ierr = i_err_fatal
        RETURN
      END IF
      IF (line(1:11) == "*BEGIN_DATA") EXIT
    END DO

    ! Set up the size of the table
    READ(iu_ckd_data, "(3(/, T38, E12.5))") &
      c_self_h2o_296%table_start, &
      c_self_h2o_296%table_end, &
      c_self_h2o_296%table_inc
    READ(iu_ckd_data, "(T38, i5)") c_self_h2o_296%n_freq
    ALLOCATE(c_self_h2o_296%c(c_self_h2o_296%n_freq))
    READ(iu_ckd_data, "(//, 4(3X, 1E12.5))", IOSTAT=ios) &
      (c_self_h2o_296%c(k), k=1, c_self_h2o_296%n_freq)
    IF (ios /= 0) THEN
      WRITE(iu_err, "(/a)") "*** Error reading self-broadened continuum."
      ierr = i_err_fatal
      RETURN
    END IF

    CLOSE(iu_ckd_data)

  END IF


  WRITE(*, '(/a)') &
    'Enter the name of the file containing the self-broadened CKD data at 260K'
  WRITE(*, '(a)') 'or the netCDF file with the temperature exponent.'
  READ(iu_stdin, '(a)', iostat=ios) file_in
  IF (ios /= 0) THEN
    WRITE(iu_err, "(/a)") '*** Error reading self-broadened continuum file name'
    ierr = i_err_fatal
    RETURN
  END IF
  INQUIRE(FILE=TRIM(file_in), EXIST=l_exist)
  IF (.NOT.l_exist) THEN
    WRITE(iu_err, '(a)') '*** The self-broadened continuum file does not exist'
    ierr=i_err_exist
    RETURN
  END IF

  IF (INDEX(file_in, '.nc') > 0) THEN

    c_self_h2o_260%n_freq = c_self_h2o_296%n_freq
    c_self_h2o_260%table_start = c_self_h2o_296%table_start
    c_self_h2o_260%table_end = c_self_h2o_296%table_end
    c_self_h2o_260%table_inc = c_self_h2o_296%table_inc

    ! Get self-broadened continuum temperature exponent
    ALLOCATE(c_self_h2o_296%texp(c_self_h2o_296%n_freq))
    CALL nf(nf90_open(TRIM(file_in), NF90_NOWRITE, ncid))
    CALL nf(nf90_inq_varid(ncid, 'self_texp', varid))
    CALL nf(nf90_get_var(ncid, varid, c_self_h2o_296%texp))
    CALL nf(nf90_close(ncid))

    ! Scale the reference continuum to 260K
    ALLOCATE(c_self_h2o_260%c(c_self_h2o_260%n_freq))
    c_self_h2o_260%c = c_self_h2o_296%c &
      * (296.0_RealK/260.0_RealK)**c_self_h2o_296%texp

    ! Uncomment to write out data in old file format
    ! WRITE(260, "(a, a)") 'Self continuum data at 260K from ', TRIM(file_in)
    ! WRITE(260, "(a/)") '*BEGIN_DATA'
    ! WRITE(260, "(a, 1pe12.5, a)") '     Start of table                = ', &
    !   c_self_h2o_260%table_start, ' cm-1'
    ! WRITE(260, "(a, 1pe12.5, a)") '     End of table                  = ', &
    !   c_self_h2o_260%table_end, ' cm-1'
    ! WRITE(260, "(a, 1pe12.5, a)") '     Increment in table            = ', &
    !   c_self_h2o_260%table_inc, ' cm-1'
    ! WRITE(260, "(a, i5)") '     Number of points in table     = ', &
    !   c_self_h2o_260%n_freq
    ! WRITE(260, "(/a)") '     Continuum coefficients (cm^3.molec-1 * 1.0E-20)'
    ! WRITE(260, "(4(1pe15.5))") c_self_h2o_260%c
    ! WRITE(260, "(a/)") '*END'

  ELSE

    CALL get_free_unit(ierr, iu_ckd_data)
    OPEN(UNIT=iu_ckd_data, FILE=TRIM(file_in), iostat=ios, status='old')
    DO 
      READ(iu_ckd_data, "(a)", IOSTAT=ios) line
      IF (ios /= 0) THEN
        WRITE(iu_err, "(/a)") "*** Error reading self-broadened continuum."
        ierr = i_err_fatal
        RETURN
      END IF
      IF (line(1:11) == "*BEGIN_DATA") EXIT
    END DO

    ! Set up the size of the table
    READ(iu_ckd_data, "(3(/, T38, 1E12.5))") &
      c_self_h2o_260%table_start, &
      c_self_h2o_260%table_end, &
      c_self_h2o_260%table_inc
    READ(iu_ckd_data, "(T38, i5)") c_self_h2o_260%n_freq
    ALLOCATE(c_self_h2o_260%c(c_self_h2o_260%n_freq))
    READ(iu_ckd_data, "(//, 4(3X, 1E12.5))", IOSTAT=ios) &
      (c_self_h2o_260%c(k), k=1, c_self_h2o_260%n_freq)
    IF (ios /= 0) THEN
      WRITE(iu_err, "(/a)") "*** Error reading self-broadened continuum."
      ierr = i_err_fatal
      RETURN
    END IF

    CLOSE(iu_ckd_data)

  END IF

CONTAINS

  SUBROUTINE nf(status)
    USE netcdf
    INTEGER, INTENT(IN):: status
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       STOP 'STOPPED!'
    END IF
  END SUBROUTINE nf

END SUBROUTINE set_extern_ckd_self_data

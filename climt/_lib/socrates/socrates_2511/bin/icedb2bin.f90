! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert the ice database to a binary format
!
PROGRAM icedb2bin
!
! Description:
!    This program receives the ASCII database of scattering 
! properties and convertis it to a direct access unformatted
! file for use in calculating single scattering properties
! averaged over distributions.
!
! Method:
!    The file is opened and read to determine how many blocks
! of data for a specific size and wavelength are contained 
! within it. An unformatted direct access file is then
! opened and the data are written to it. 
!
! Note:
!    The initial data are expected to be supplied in the 
! following units:
!    Wavelength:                   Micron
!    Mean maximum dimension:       Micron
!    Extinction cross-section:     m^2
!    Albedo of single scattering:  m^2
!
!
! Modules used:
  USE realtype_rd
  USE def_sct_db
  USE rad_ccf, ONLY: pi
  USE error_pcf
  USE def_std_io_icf
!
!
  IMPLICIT NONE
!
!
  INTEGER  ::  ierr
!                Error flag
  INTEGER  ::  ios
!                I/O error flag
  INTEGER  ::  iunit_in
!                Unit number for reading the ASCII database
  INTEGER  ::  iunit_unf
!                Unit number for writing the unformatted compressed
!                database
  INTEGER  ::  n_block
!                Number of blocks of scattering data
  INTEGER  ::  n_rec_header
!                Number of records in the header of the direct
!                access file
  INTEGER  ::  n_r4_per_rec
!                Number of 4-byte reals in a record
  INTEGER  ::  n_angle
!                Number of angles at which the phase function is
!                specified
  INTEGER  ::  i_block
!                Loop variable
  INTEGER  ::  i_data_rec
!                Record of the data file where the phase function
!                relative to the header line is stored
  INTEGER  ::  n_full_rec
!                Number of full records to hold the phase function
  LOGICAL  ::  l_partial_rec
!                Logical to indicate the need for a partial record
  INTEGER  ::  i
!                Loop variable
  INTEGER  ::  k
!                Loop variable
!
  CHARACTER(LEN=80)  ::  file_unf
!                          Name of the output file
  CHARACTER(LEN=80)  ::  line
!                          Line of text read from the input file
!
  TYPE (str_sctdb_entry)  ::  ice_sct
!                               Structure holding the ice scattering
!                               data
!
! Subroutines called:
  EXTERNAL get_free_unit, open_file_in
!
!- End of header
!
!
!
! Set generic sizes.
  n_r4_per_rec=sct_db_recl/8
!
! Check the size of the record
  IF (n_r4_per_rec < 10) THEN
    WRITE(iu_err, '(/A)') &
      '*** Error: The predefined record length is too small.'
    STOP
  ENDIF
!
! Get a unit and read the formatted database.
  CALL get_free_unit(ierr, iunit_in)
  IF (ierr /= I_NORMAL) STOP
  CALL open_file_in(ierr, iunit_in, &
    'Enter the name of the ASCII database.')
  IF (ierr /= I_NORMAL) STOP

! Open an unformatted file for output.
  CALL get_free_unit(ierr, iunit_unf)
  IF (ierr /= I_NORMAL) STOP
  WRITE(iu_stdout, '(A)') 'Enter the name of the compressed output file.'
  READ(iu_stdin, '(A)') file_unf
  OPEN(UNIT=iunit_unf, FILE=file_unf, FORM='UNFORMATTED', &
    ACCESS='DIRECT', RECL=sct_db_recl)
!
! Read through the input to find the number of blocks of data.
! In the current format the first element of an entry is the wavelength.
  n_block=0
  count_block:  DO 
    READ(iunit_in, '(A)', IOSTAT=ios) line
!   Conventionally negative errors denote an end of the file.
    IF (ios < 0) EXIT
    IF (line(17:28) == "; Wavelength") n_block=n_block+1
  ENDDO count_block
  REWIND(iunit_in)
!
!
! Calculate the space required for the header: one record
! holds the number of blocks; this is followed by one record
! for each block of data.
  n_rec_header=n_block+1
! Records of the data begin after the header.
  i_data_rec=n_rec_header+2
!
! Write the total number of blocks to the unformatted file.
  WRITE(iunit_unf, REC=1) n_block
!
!
! Read each block of data and convert to the binary format.
  i_block=1
  process_block:  DO
!
    READ(iunit_in, '(A)') line
    READ(line(1:16), *) ice_sct%wavelength
!   Convert to SI units.
    ice_sct%wavelength=1.0E-06_Real4*ice_sct%wavelength
!
    READ(iunit_in, '(A)') line
    READ(line(1:16), *) ice_sct%dm
!   Convert to SI units.
    ice_sct%dm=1.0E-06_Real4*ice_sct%dm
!
    READ(iunit_in, '(A)') line
    READ(line(1:16), *) ice_sct%cext
    READ(iunit_in, '(A)') line
    READ(line(1:16), *) ice_sct%omega
!
    READ(iunit_in, '(A)') line
    READ(line(1:16), *) ice_sct%asymm
!
!   The scattering cross-section is redundant, but is kept for historical
!   consistency.
    ice_sct%csca = ice_sct%cext * ice_sct%omega
!
!
!   The phase function:
    n_angle=0
    READ(iunit_in, '()')
    read_phase:  DO
!     If the line contains the string  "; Wavelength" it signals 
!     the start of the next block, or the block may be at the end 
!     of the file.
      READ(iunit_in, '(A)', IOSTAT=ios) line
      IF ( (ios < 0).OR.(line(17:28) == '; Wavelength') ) THEN
        BACKSPACE(iunit_in)
        EXIT
      ENDIF
      n_angle=n_angle+1
      IF (n_angle > npd_sct_db_angle) THEN
        WRITE(iu_err, '(/A, I5)') &
          '*** Error: Too many angles in the phase function in '// &
          'block ', i_block
        STOP
      ENDIF
      BACKSPACE(iunit_in)
      READ(iunit_in, *) ice_sct%phf_angle(n_angle), ice_sct%phf(n_angle)
!     Convert to SI units.
      ice_sct%phf_angle(n_angle)=(pi/180.0)*ice_sct%phf_angle(n_angle)
      ice_sct%n_angle=n_angle
!
!     Determine whether the angles need to be stored explicitly by
!     testing whether they are uniformly distributed.
      IF (n_angle == 2) THEN
        ice_sct%l_uniform=.TRUE.
        ice_sct%d_angle=ice_sct%phf_angle(2)-ice_sct%phf_angle(1)
      ELSE IF (n_angle > 2) THEN
!       Check that the differences in the increment are sufficiently
!       small for uniformity.
        IF ( ABS(ice_sct%phf_angle(n_angle) - &
                 ice_sct%phf_angle(n_angle-1) - &
                 ice_sct%d_angle) > 1.0E-02*ice_sct%d_angle) THEN
        WRITE(*,*) n_angle, ABS(ice_sct%phf_angle(n_angle) - &
                                ice_sct%phf_angle(n_angle-1) - &
                                ice_sct%d_angle)
          ice_sct%l_uniform=.FALSE.
        ENDIF
      ENDIF
    ENDDO read_phase
!   At least two angles are required:
    IF (n_angle < 2) THEN
      WRITE(iu_err, '(/A, /A, I5)') &
        '*** Error: At least two angles are required to specify ', &
        'the phase function: error in block ', i_block
      STOP
    ENDIF
!
!
!   Write the data to the direct access file, first completing the
!   header.
    WRITE(iunit_unf, REC=i_block+1) &
      ice_sct%dm, ice_sct%wavelength, ice_sct%csca, ice_sct%cext, &
      ice_sct%omega, ice_sct%asymm, ice_sct%n_angle, &
      ice_sct%l_uniform, ice_sct%phf_angle(1), ice_sct%d_angle, &
      i_data_rec
!   Now write the phase function.
    n_full_rec=ice_sct%n_angle/n_r4_per_rec
    l_partial_rec=((ice_sct%n_angle-n_r4_per_rec*n_full_rec) > 0)
    IF (.NOT.ice_sct%l_uniform) THEN
!     We need to write the angles as well.
      DO i=1, n_full_rec
        WRITE(iunit_unf, REC=i_data_rec) &
          (ice_sct%phf_angle(k), k=1+(i-1)*n_r4_per_rec, i*n_r4_per_rec)
        i_data_rec=i_data_rec+1
      ENDDO
      IF (l_partial_rec) THEN
        WRITE(iunit_unf, REC=i_data_rec) &
          (ice_sct%phf_angle(k), &
           k=1+n_full_rec*n_r4_per_rec, ice_sct%n_angle)
        i_data_rec=i_data_rec+1
      ENDIF
    ENDIF
!   Write the phase function.
    DO i=1, n_full_rec
      WRITE(iunit_unf, REC=i_data_rec) &
        (ice_sct%phf(k), k=1+(i-1)*n_r4_per_rec, i*n_r4_per_rec)
      i_data_rec=i_data_rec+1
    ENDDO
    IF (l_partial_rec) THEN
      WRITE(iunit_unf, REC=i_data_rec) &
        (ice_sct%phf(k), &
         k=1+n_full_rec*n_r4_per_rec, ice_sct%n_angle)
      i_data_rec=i_data_rec+1
    ENDIF
!
!   Advance the count of the number of blocks.
    i_block=i_block+1
!   Stop at the last block.
    IF (i_block > n_block) EXIT
!
  ENDDO process_block
!
!
!
  STOP
  END



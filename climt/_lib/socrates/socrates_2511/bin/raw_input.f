! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to convert data to CDL format.
!
!
! Method:
!   An input file with embedded directives is read. This file
!   is processed to generate a number of profiles.
!
!- ---------------------------------------------------------------------
      PROGRAM raw_input
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
      USE dimensions_spec_ucf
      USE gas_list_pcf
      USE input_head_pcf
      USE def_data_in_icf
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables
      INTEGER
     &    ierr
!           Error flag
     &  , length_base
!           Length of base file name
     &  , n_profile
!           Number of profiles read in
     &  , i_zero
!           Position of 0 in collating seq.
     &  , ios
!           I/O error flag
     &  , n_level_profile(npd_in_profile)
!           Number of levels in profiles
     &  , n_column_profile(npd_in_profile)
!           Number of columns in profiles
     &  , i_data_group(npd_data_column, npd_in_profile)
!           Groups of types of data in columns
     &  , i_data_type(npd_data_column, npd_in_profile)
!           Types of data within groups
     &  , length_name_profile
!           Length of profile name
     &  , n_level_ref
!           Number of levels in reference
     &  , i
!           Loop variable
     &  , j
!           Loop variable
      LOGICAL
     &    l_reference
!           True if reference profile present
     &  , l_remove_missing
!           True to remove missing data
      REAL  (RealK) ::
     &    profile(npd_layer+1, npd_phys_type, npd_in_profile)
!           Initial profiles
     &  , p_ref(npd_layer+1)
!           Reference pressure
     &  , z_ref(npd_layer+1)
!           Reference height
     &  , missing_data_flag
!           Missing data flag
      CHARACTER
     &    base_name*80
!           Base name for output files
     &  , name_profile*80
!           Name of profile
     &  , line*80
!           Line of input
     &  , char_yn*1
!           Character response variable
!
      data ierr/i_normal/
!     Subroutines called:
      EXTERNAL
     &    open_file_in, read_raw_profile, sort_raw_profile
     &  , set_state, write_profile
!
!
!
!     Set the error flag to the normal value.
      ierr=i_normal
!
!     Find the location of 0 in the compiler's collating sequence.
      i_zero=ichar('0')
!
!     Open input file.
      CALL open_file_in(ierr, iu_raw_in
     &  , 'Enter the name of the raw input file.')
      IF (ierr /= i_normal) STOP
!
      WRITE(iu_stdout, '(a)')
     &  'Enter the base name of the output file.'
      READ(iu_stdin, '(a)') base_name
!     Find length of base name.
      j=0
      length_base=0
1     j=j+1
      IF (base_name(j:j) /= ' ') THEN
        length_base=length_base+1
        goto 1
      ENDIF
!
!     Set the treatment of missing data.
      WRITE(iu_stdout, '(/a)') 'Do you want to remove
     &                                missing data? (y/n)'
2     read(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
        l_remove_missing=.true.
        WRITE(iu_stdout, '(/a)') 'Enter the missing data flag.'
3       read(iu_stdin, *, iostat=ios) missing_data_flag
        IF (ios /= 0) THEN
          WRITE(iu_err, '(/a)')
     &      '+++ Unrecognized response: please re-enter.'
          goto 3
        ENDIF
      ELSE IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
        l_remove_missing=.false.
      ELSE
        WRITE(iu_err, '(/a)') '+++ Illegal response: please re-enter.'
        goto 2
      ENDIF
!
      n_profile=0
!     Read each line until a directive is found.
4     read(iu_raw_in, '(a)', end=5) line
      IF (line(1:8) == '*PROFILE') THEN
        n_profile=n_profile+1
        CALL read_raw_profile(ierr
     &    , l_remove_missing, missing_data_flag
     &    , n_level_profile(n_profile), n_column_profile(n_profile)
     &    , i_data_group(1, n_profile), i_data_type(1, n_profile)
     &    , profile(1, 1, n_profile)
     &    )
        IF (ierr /= i_normal) STOP
        CALL sort_raw_profile(ierr
     &    , n_level_profile(n_profile), n_column_profile(n_profile)
     &    , i_data_type(1, n_profile), i_data_group(1, n_profile)
     &    , profile(1, 1, n_profile)
     &    )
        IF (ierr /= i_normal) STOP
      ELSE IF (line(1:10) == '*ATM_STATE') THEN
        CALL set_state(ierr
     &    , n_level_profile(n_profile), n_column_profile(n_profile)
     &    , i_data_type(1, n_profile), profile(1, 1, n_profile)
     &    , l_reference, n_level_ref, p_ref, z_ref
     &    )
        IF (ierr /= i_normal) STOP
      ELSE IF (line(1:1) == '*') THEN
        WRITE(iu_err, '(/a, /a)') 
     &    '*** Error: Unrecognized *-directive in input.'
     &    , 'Program aborted.'
        STOP
      ENDIF
      goto 4
!
5     close(iu_raw_in)
!
!     Write out each profile in turn.
      DO i=1, n_profile
        name_profile(1: length_base+2)=base_name(1: length_base)
     &    //'_'//char(i+i_zero)
        length_name_profile=length_base+2
        CALL write_profile(ierr
     &    , n_level_profile(i), n_column_profile(i)
     &    , i_data_group(1, i), i_data_type(1, i), profile(1, 1, i)
     &    , l_reference, n_level_ref, p_ref, z_ref
     &    , name_profile, length_name_profile
     &    , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &    )
        IF (ierr /= i_normal) STOP
      ENDDO
!
!
!
      STOP
      END

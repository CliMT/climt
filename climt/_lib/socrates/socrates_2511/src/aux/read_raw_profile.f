! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a single profile of raw data.
!
! Method:
!   The header is read to determine the contents. These are
!   read and converted to S. I. units.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_raw_profile(ierr
     &  , l_remove_missing, missing_data_flag
     &  , n_level_profile
     &  , n_column_profile, i_data_group, i_data_type, profile)
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_data_in_icf
      USE def_std_io_icf
      USE gas_list_pcf
      USE dimensions_spec_ucf
      USE rad_pcf
      USE dimensions_field_ucf
      USE input_head_pcf
!
!
      IMPLICIT NONE
!
!
!     Include headers
!
!     Dummy arguments
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , i_data_group(npd_data_column)
!           Groups of data in profile
     &  , i_data_type(npd_data_column)
!           Types of data in profile
     &  , n_column_profile
!           Number of columns in profile
     &  , n_level_profile
!           Number of levels in profile
      LOGICAL, Intent(IN) ::
     &    l_remove_missing
!           Logical to remove missing data
      REAL  (RealK), Intent(IN) ::
     &    missing_data_flag
!           Missing data flag
      REAL  (RealK), Intent(OUT) ::
     &    profile(npd_layer+1, npd_data_column)
!           Data of profile
!
!     Local arguments.
      CHARACTER
     &    line*132
!           Input line
     &  , word(npd_data_column)*14
!           Header words
     &  , type(npd_data_column)*10
!           Type of data
     &  , tab*1
!           Tab character
      CHARACTER(LEN=len_col_header) :: unit(npd_data_column)
!           Unit of data
      LOGICAL
     &    l_discard
!           Flag to discard line
      INTEGER
     &    i
!           Loop variable
     &  , j
!           Loop variable
     &  , k
!           Loop variable
     &  , n_word
!           Number of words in header line
     &  , length_word(npd_data_column)
!           Length of each header word
     &  , length_type
!           Length of name of type
     &  , length_unit
!           Length of name of unit
     &  , ios
!           I/O error flag
      REAL  (RealK) ::
     &    conversion_factor(npd_data_column)
!           Conversions to S.I. units
     &  , conversion_offset(npd_data_column)
!           Conversions to S.I. units
!
!     Include data statements.
!
      data tab/'	'/
!
!
!     Initialize character variables.
      DO j=1, len(line)
        line(j:j)=' '
      ENDDO
      DO i=1, npd_data_column
        DO j=1, len(type(1))
          type(i)(j:j)=' '
        ENDDO
        DO j=1, len(unit(1))
          unit(i)(j:j)=' '
        ENDDO
      ENDDO
!
!     Read the header line and split it into words.
      READ(iu_raw_in, '(a)') line
      j=0
      n_word=0
      DO i=1, npd_data_column
        length_word(i)=0
      ENDDO
1     j=j+1
      IF (j == 1) THEN
        IF ( (line(j:j) /= ' ').AND.(line(j:j) /= tab) ) THEN
!         A new word is found.
          n_word=n_word+1
        END IF
      ELSE IF ( ( (line(j:j) /= ' ').AND.(line(j:j) /= tab) ).AND.
     &      ( (line(j-1:j-1) == ' ').OR. (line(j-1:j-1) == tab) ) ) THEN
!       A new word is found.
        n_word=n_word+1
      ENDIF
      IF ( (line(j:j) /= ' ').AND.(line(j:j) /= tab) ) THEN
        k=length_word(n_word)+1
        length_word(n_word)=k
        word(n_word)(k:k)=line(j:j)
      ENDIF
      IF (j < np_max_length_line) goto 1
!
      n_column_profile=n_word
!
!     Process each word to find the type and the unit used.
      DO i=1, n_word
        j=0
        length_type=0
2       j=j+1
        IF (word(i)(j:j) /= '(') THEN
          type(i)(j:j)=word(i)(j:j)
          length_type=length_type+1
          goto 2
        ENDIF
!       Now find the unit.
        k=0
        length_unit=0
3       j=j+1
        k=k+1
        IF (word(i)(j:j) /= ')' ) THEN
          unit(i)(k:k)=word(i)(j:j)
          length_unit=length_unit+1
          goto 3
        ENDIF
!
!       Check that each type is valid and find the conversion factor
!       to S.I. units.
        k=0
4       k=k+1
        IF (type(i) == header_phys(k)) THEN
!         Column contains physical data.
          i_data_group(i)=IP_physical_data
          i_data_type(i)=k
        ELSE IF (k < npd_phys_type) THEN
          goto 4
        ELSE
!         The column does not contain physical data so check gaseous data.
          k=0
5         k=k+1
          IF (type(i) == header_gas(k)) THEN
!           Column contains gaseous data.
            i_data_group(i)=IP_gaseous_data
            i_data_type(i)=k
          ELSE IF (k < npd_gases) THEN
            goto 5
          ELSE
!           The column does not contain gaseous data either
!           so we check aerosol data.
            k=0
6           k=k+1
            IF (type(i) == header_aerosol(k)) THEN
!             Column contains physical data.
              i_data_group(i)=IP_aerosol_data
              i_data_type(i)=k
            ELSE IF (k < npd_aerosol_component) THEN
              goto 6
            ELSE
              WRITE(iu_err, '(/, 2a)') 
     &          '*** Error: Undefined type: ', type(i)
              ierr=i_err_fatal
              RETURN
            ENDIF
          ENDIF
        ENDIF
!
!       Set the conversion factor and offset
        k=0
7       k=k+1
        IF (unit(i) == name_unit(k)) THEN
          IF (factor_unit(k) < 0.0_RealK) THEN
            conversion_factor(i)=-factor_unit(k)
     &        *molar_weight(i_data_type(i))
          ELSE
            conversion_factor(i)=factor_unit(k)
          END IF
          conversion_offset(i)=offset_unit(k)
        ELSE IF (k < npd_unit) THEN
          goto 7
        ELSE
          WRITE(iu_err, '(/, 2a)') 
     &      '*** Error: Unknown unit: ', unit(i)
          ierr=i_err_fatal
          RETURN
        ENDIF
!
      ENDDO
!
!     Read in the profile, converting the values to SI units.
      j=0
8     read(iu_raw_in, '(a)') line
      IF (line(1:4) == '*END') THEN
        n_level_profile=j
        RETURN
      ELSE IF (j < npd_layer+1) THEN
        j=j+1
        backspace(iu_raw_in)
        READ(iu_raw_in, *, iostat=ios) (profile(j, i), i=1, n_word)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: missing *END directive.'
          ierr=i_err_fatal
          RETURN
        ENDIF
!
!       Check for missing data and skip this line if it is found and
!       to be removed.
        IF (l_remove_missing) THEN
          l_discard=.false.
          DO i=1, n_word
          IF (abs(profile(j, i)-missing_data_flag) < 
     &      1.0e-02_RealK*abs(missing_data_flag)) l_discard=.true.
          ENDDO
          IF (l_discard) THEN
            j=j-1
            goto 8
          ENDIF
        ENDIF
!       Carry out any required conversions.
        DO i=1, n_word
          profile(j, i)=profile(j, i)*conversion_factor(i)
     &      +conversion_offset(i)
        ENDDO
        goto 8
      ELSE
        WRITE(iu_err, '(/a)')
     &    '*** Error: Profile is too large: increase npd_layer and '
     &    //'recompile.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!
!
      END

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to sort the elements of a profile in increasing pressure.
!
! Method:
!   The profile is searched to find the pressure or height on
!   which to sort. A shell sorting algorithm is called to find
!   the correct order. Using the key thus defined the columns
!   os the profile are re-ordered.
!
!- ---------------------------------------------------------------------
      SUBROUTINE sort_raw_profile(ierr
     &  , n_level_profile, n_column_profile, i_data_type, i_data_group
     &  , profile
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE input_head_pcf
      USE def_std_io_icf
      USE error_pcf
      USE shell_sort_mod, ONLY: shell_sort
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n_level_profile
!           Number of levels
     &  , n_column_profile
!           Number of columns
     &  , i_data_type(npd_data_column)
!           Types of data in columns
     &  , i_data_group(npd_data_column)
!           Groups of types of data in columns
      REAL  (RealK), Intent(INOUT) ::
     &    profile(npd_layer+1, npd_data_column)
!           Profile of data
!
!     Local variables.
      INTEGER
     &    i_column
!           Column being searched
     &  , i_column_pressure
!           Column containing pressure
     &  , i_column_height
!           Column containing height
     &  , index_level(npd_layer+1)
!           Index for sorting
     &  , i
!           Loop variable
     &  , j
!           Loop variable
      REAL  (RealK) ::
     &    sort_key(npd_layer+1)
!           Sorting key
     &  , x(npd_layer+1)
!           Temporary variable
!
!
!     Proceed through the profile to find the pressure or the height.
      i_column_pressure=0
      i_column=0
1     i_column=i_column+1
      IF ( (i_data_type(i_column) == IP_pressure).AND.
     &     (i_data_group(i_column) == IP_physical_data) ) THEN
        i_column_pressure=i_column
!       The sorting key is the pressure.
        DO i=1, n_level_profile
          sort_key(i)=profile(i, i_column_pressure)
          index_level(i)=i
        ENDDO
      ELSE IF (i_column < n_column_profile) THEN
        goto 1
      ENDIF
!
      IF (i_column_pressure == 0) THEN
!       There is no pressure in this profile: search for the height.
        i_column_height=0
        i_column=0
2       i_column=i_column+1
        IF ( (i_data_type(i_column) == IP_height).AND.
     &       (i_data_group(i_column) == IP_physical_data) ) THEN
          i_column_height=i_column
          DO i=1, n_level_profile
            sort_key(i)=-profile(i, i_column_height)
            index_level(i)=i
          ENDDO
        ELSE IF (i_column < n_column_profile) THEN
          goto 2
        ELSE
          WRITE(iu_err, '(/a)') 
     &      '*** Error: This profile contains neither '
     &      //'pressure nor height data.'
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
!
!     Sort the array index_level in order of ascending value of sort_key.
      IF (n_level_profile > 1) THEN
        CALL shell_sort(n_level_profile, index_level, sort_key)
      ELSE
        sort_key(1)=1
      ENDIF
!
!
!     Pass through the profile re-ordering each column as required by
!     the index.
      DO j=1, n_column_profile
        DO i=1, n_level_profile
          x(i)=profile(i, j)
        ENDDO
        DO i=1, n_level_profile
          profile(i, j)=x(index_level(i))
        ENDDO
      ENDDO
!
!
!
      RETURN
      END

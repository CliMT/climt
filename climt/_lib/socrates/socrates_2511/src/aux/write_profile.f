! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to split a profile into CDL-files.
!
! Method:
!   For each field in a profile is examined. Appropriate pressure
!   coordinates are derived and the field is written to a file
!   with a suitable suffix.
!
!- ---------------------------------------------------------------------
      SUBROUTINE write_profile(ierr
     &  , n_level, n_column_profile, i_data_group, i_data_type, profile
     &  , l_reference, n_level_ref, p_ref, z_ref
     &  , name_profile, length
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE dimensions_field_ucf
      USE gas_list_pcf
      USE rad_pcf
      USE input_head_pcf
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
!
!     Declaration of data structures
      INCLUDE 'cdl_struc.finc'
!
      INTEGER, Intent(IN) ::
     &    n_level
!           Number of levels in profile
     &  , n_column_profile
!           Number of columns in profile
     &  , i_data_group(npd_data_column)
!           Data types of columns
     &  , i_data_type(npd_data_column)
!           Data types of columns
     &  , length
!           Length of profile name
     &  , n_level_ref
!           Number of levels in reference atm.
      CHARACTER(LEN=80), Intent(IN) ::
     &    name_profile
!           Name of profile
      LOGICAL, Intent(IN) ::
     &    l_reference
!           True if reference profile present
      REAL  (RealK), Intent(IN) ::
     &    profile(npd_layer+1, npd_data_column)
!           Data for output
     &  , p_ref(npd_layer+1)
!           Reference atmospheric pressure
     &  , z_ref(npd_layer+1)
!           Reference atmospheric heights
!
!     Local variables.
      INTEGER
     &    i_column
!           Column index
     &  , i_type
!           Type index
     &  , i
!           Loop variable
     &  , j
!           Loop variable
      LOGICAL
     &    l_pressure
!           Pressure flag
     &  , l_height
!           Height flag
      CHARACTER(LEN=80) ::
     &    file_name
!           Name of output file
     &  , empty = REPEAT(' ',80)
!           Empty character string to fill names
      REAL  (RealK) ::
     &    field(npd_profile, npd_layer+1)
!           Output field
     &  , p(npd_profile, npd_layer+1)
!           Output pressure
!
!     Subroutines called:
      EXTERNAL
     &    write_cdl_field


!     Initialize the logicals.
      l_pressure=.false.
      l_height=.false.
!
!     Scan the list of data types to find the pressure, if present, or
!     the height if not.
      i_column=0
1     i_column=i_column+1
      IF ( (i_data_group(i_column) == IP_physical_data).AND.
     &  (i_data_type(i_column) == IP_pressure) ) THEN
        l_pressure=.true.
        DO i=1, n_level
          p(1, i)=profile(i, i_column)
        ENDDO
      ELSE IF (i_column < n_column_profile) THEN
        goto 1
      ENDIF
!
      IF (.NOT.l_pressure) THEN
!       The pressure must be interpolated from the height.
!       check that a reference profile is set.
        IF (.NOT.l_reference) THEN
          WRITE(iu_err, '(/a)')
     &      '** Error: No reference profile is set.'
          ierr=i_err_fatal
          RETURN
        ENDIF
        i_column=0
2       i_column=i_column+1
        IF ( (i_data_group(i_column) == IP_physical_data).AND.
     &       (i_data_type(i_column) == IP_height) ) THEN
          l_height=.true.
        ELSE IF (i_column < n_column_profile) THEN
          goto 2
        ENDIF
!
        IF (.NOT.l_height) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: Current profile contains '
     &      //'no height or pressure data.'
          ierr=i_err_fatal
          RETURN
        ENDIF
!
        DO i=1, n_level
!         find the pressure levels in the reference profile bracketing
!         the desired level. note that heights in the reference are in
!         ascending order for conformity with the normal output order.
          j=n_level_ref
3         if ( (z_ref(j-1).lt.profile(i, i_column)).and.(j.gt.2) ) then
            j=j-1
            goto 3
          ENDIF
          p(1, i)=p_ref(j)*exp((profile(i, i_column)-z_ref(j))
     &      *log(p_ref(j-1)/p_ref(j))/(z_ref(j-1)-z_ref(j)))
        ENDDO
      ENDIF
!
!     Write a separate file for each column.
      DO j=1, n_column_profile
        i_type=i_data_type(j)
!       Discard specified columns.
        IF (.NOT.(i_type == IP_discard .AND.
     &            i_data_group(j) == IP_physical_data)) THEN
          DO i=1, len(file_name)
            file_name(i:i)=' '
          ENDDO
          IF (i_data_group(j) == IP_physical_data) THEN
            file_name(1: length+1+len_file_suffix)
     &        =name_profile(1: length)//'.'//phys_suffix(i_type)
            var_name(1)(1:24)=phys_suffix(i_type)
     &        //empty(1:24-len(phys_suffix(i_type)))
            var_long(1)(1: 30)=phys_title(i_type)(1: 30)
          ELSE IF (i_data_group(j) == IP_gaseous_data) THEN
            file_name(1: length+1+len_file_suffix)
     &        =name_profile(1: length)//'.'//gas_suffix(i_type)
            var_name(1)(1:24)=gas_suffix(i_type)
     &        //empty(1:24-len(gas_suffix(i_type)))
            var_long(1)(1: 20)=name_absorb(i_type)(1: 20)
          ELSE IF (i_data_group(j) == IP_aerosol_data) THEN
            file_name(1: length+1+len_file_suffix)
     &        =name_profile(1: length)//'.'//aerosol_suffix(i_type)
            var_name(1)(1:24)=aerosol_suffix(i_type)
     &        //empty(1:24-len(aerosol_suffix(i_type)))
            var_long(1)(1: 30)=aerosol_title(i_type)(1: 30)
          ENDIF
          DO i=31, len(var_long)
            var_long(1)(i: i)=' '
          ENDDO
          DO i=1, n_level
            field(1, i)=profile(i, j)
          ENDDO
!         Set up common features of the dimensions.
          n_dimension=3
          dimension_name(1)='lat'
          dimension_name(2)='lon'
          dimension_name(3)='plev'
          dimension_type(1)='float'
          dimension_type(2)='float'
          dimension_type(3)='float'
          dimension_unit(1)='degree'
          dimension_unit(2)='degree'
          dimension_unit(3)='pa'
          dimension_long(1)='latitude'
          dimension_long(2)='longitude'
          dimension_long(3)='pressure'
          dimension_size(1)=1
          dimension_size(2)=1
          dimension_size(3)=n_level
!
          dimension_array_fl(1, 1)=0.0_RealK
          dimension_array_fl(1, 2)=0.0_RealK
          DO i=1, n_level
            dimension_array_fl(i, 3)=p(1, i)
          ENDDO
!
          n_var=1
          var_type(1)(1:6)='float '
          var_unit(1)='none'
          n_data(1)=n_level
          n_dimension_var(1)=3
          list_dimension_var(1, 1)=1
          list_dimension_var(2, 1)=2
          list_dimension_var(3, 1)=3
          DO i=1, n_level
            data_fl(i, 1)=field(1, i)
          ENDDO
!
          CALL write_cdl(ierr
     &      , file_name
     &      , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &      , n_dimension, dimension_name
     &      , dimension_type, dimension_unit
     &      , dimension_long, dimension_size
     &      , dimension_array_int, dimension_array_fl
     &      , n_var, var_name, var_type, var_unit, var_long
     &      , n_dimension_var, list_dimension_var
     &      , n_data, data_int, data_fl
     &      )
          IF (ierr /= i_normal) RETURN
        ENDIF
      ENDDO
!
!
!
      RETURN
      END

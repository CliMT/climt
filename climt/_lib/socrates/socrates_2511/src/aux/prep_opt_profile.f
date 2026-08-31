! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to make an observed profile of optical properties.
!
!
! Description:
!   This program combines files of averaged scattering data into
!   a single CDL-file.
!
! Method:
!   A File listing vertical coordinates and names of files of
!   scattering properties is read. The vertical coordinates are
!   sorted and a combined profile is written as a CDL file, ready
!   to be picked up by the radiation code.
!
!- ---------------------------------------------------------------------
      PROGRAM prep_opt_profile
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
      USE dimensions_cdl_ucf
      USE dimensions_field_ucf
      USE dimensions_spec_ucf
      USE shell_sort_mod, ONLY: shell_sort
!
!
      IMPLICIT NONE
!
!
!
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
!
!
      CHARACTER (LEN=80), ALLOCATABLE :: filename(:)
!           Names of the files holding scattering data
      LOGICAL
     &    l_end
!           Flag for the end of a file
      INTEGER
     &    iunit_list
!           Unit number for reading the list of input files
     &  , n_file
!           Number of input files
     &  , i
!           Loop variable
     &  , nd_input_file
!           Size allocated for input files.
      INTEGER, ALLOCATABLE :: pointer(:)
!           Array of pointers to sets of input data sorted by pressure
!
      REAL  (RealK), ALLOCATABLE :: p(:)
!           Pressures at which scattering data are provided
      REAL  (RealK) ::
     &    latitude
!           Latitude of data
     &  , longitude
!           Latitude of data
!
      EXTERNAL
     &    get_free_unit, open_file_in
     &  , output_opt_profile_cdl
      data ierr/i_normal/
!
!
!
!     Set a wodking limit on the number of input files.
      nd_input_file=npd_opt_level_cloud_prsc
     &  +npd_opt_level_aerosol_prsc+1
!
      ALLOCATE ( filename(nd_input_file) )
      ALLOCATE ( pointer (nd_input_file) )
      ALLOCATE ( p       (nd_input_file) )

!     Determine and open the file listing the scattering data.
      CALL get_free_unit(ierr, iunit_list)
      CALL open_file_in(ierr, iunit_list
     &  , 'give the name of the file listing scattering data.')
!
!     Read the pressures and the names of the corresponding files.
      n_file=0
      l_end=.false.
      DO WHILE(.NOT.l_end)
        n_file=n_file+1
        IF (n_file > nd_input_file) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: The number of input files is too big.'
          STOP
        ENDIF
        READ(iunit_list, *, iostat=ios) p(n_file), filename(n_file)
        IF (ios /= 0) l_end=.true.
      ENDDO
      n_file=n_file-1
!
!     To cover a range of heights at least 2 files are needed.
      IF (n_file < 2) THEN
        WRITE(iu_err, '(/a)') 
     &    '*** Error: At least two sets of scattering data are '
     &    //'required.'
        STOP
      ENDIF
!
!     Sort the pressures into an increasing order using an array of
!     pointers, which is first set up in the default order.
      DO i=1, n_file
        pointer(i)=i
      ENDDO
      CALL shell_sort(n_file, pointer, p)
!
!     Obtain the position of the point.
      WRITE(iu_stdout, '(a)') 
     &  'Enter the latitude and longitude of the point.'
      READ(iu_stdin, *) latitude, longitude
!
!     Read in the first file to get the size of the data.
      CALL output_opt_profile_cdl(ierr
     &  , n_file, p, filename, pointer, latitude, longitude
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  , nd_input_file, npd_phase_term, npd_band)
      IF (ierr /= i_normal) STOP
!
      DEALLOCATE ( filename )
      DEALLOCATE ( pointer  )
      DEALLOCATE ( p        )
!
!
      STOP
      END

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert a samson file to a netCDF file.
!
! Purpose:
!   Reads in a Field in SAMSON format and writes it out as a netCDF field.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      PROGRAM genln2cdf
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_cdl_ucf
      USE error_pcf
      USE def_std_io_icf
!
!
      IMPLICIT NONE
!
!
!
!
      INTEGER
     &    ierr
!           Error flag
!
!     Variables for input file:
      CHARACTER
     &    file_in*80
!           Name of input file
     &  , title_in*80
!           Title of input file
      INTEGER
     &    n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
      REAL  (RealK) ::
     &    p(npd_cdl_dimen_size+1)
!           Pressure levels
     &  , pstar(npd_cdl_dimen_size)
!           Surface pressures
     &  , field(npd_cdl_dimen_size, npd_cdl_dimen_size+1)
!           Field of input data
!
!
!     Variables for output file:
      INTEGER
     &    n_dimension
!           Number of dimensions
     &  , n_var
!           Number of netCDF variables
     &  , var_size(1)
!           Size of data array
     &  , n_dimension_var(1)
!           Number of dimensions for each array
     &  , list_dimension_var(npd_cdl_dimen, 1)
!           List of dimensions for each array
      CHARACTER
     &    filename_out*80
!           Name of output file
     &  , dimension_name(npd_cdl_dimen)*24
!           Names of dimensions
     &  , dimension_type(npd_cdl_dimen)*6
!           Types of dimensions
     &  , dimension_long(npd_cdl_dimen)*40
!           Long names of dimensions
     &  , dimension_unit(npd_cdl_dimen)*20
!           Units of dimensions
     &  , var_name(1)*24
!           Names of dimensions
     &  , var_type(1)*6
!           Types of dimensions
     &  , var_long(1)*40
!           Long names of dimensions
     &  , var_unit(1)*20
!           Unit of variable
      INTEGER
     &    dimension_array_int(npd_cdl_dimen_size, npd_cdl_dimen)
!           Values of dimensions
      REAL  (RealK) ::
     &    dimension_array_fl(npd_cdl_dimen_size, npd_cdl_dimen)
!           Values of dimensions
      INTEGER
     &    dimension_size(npd_cdl_dimen)
!           Sizes of dimensions
      INTEGER
     &    var_int(npd_cdl_data, 1)
!           Integer variable
      REAL  (RealK) ::
     &    var_fl(npd_cdl_data, 1)
!           Floating-point variable
!
!
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
!
      data ierr/i_normal/
!     Subroutines called
      EXTERNAL
     &    read_samson_p_field, write_cdf
!
!
!
      WRITE(iu_stdout, '(/a)')
     &     'enter the name of the input file.'
      READ(iu_stdin, '(a)') file_in
      CALL read_genln2_flux(ierr, file_in
     &   , n_profile, n_level, p, field
     &   , npd_cdl_dimen_size, npd_cdl_dimen_size)
      IF (ierr /= i_normal) STOP
!
      n_profile=1
!
!     Write the field out in netCDF format
!
!     Set up the appropriate sizes. We map all points to the origin
!     since SAMSON files have no geographical information.
!
      n_dimension=3
!
      dimension_name(1)='lat'
      dimension_type(1)='float'
      dimension_unit(1)='degree'
      dimension_long(1)='latitude'
      dimension_size(1)=1
      DO i=1, n_profile
        dimension_array_fl(i, 1)=0.0_RealK
      ENDDO
!
      dimension_name(2)='lon'
      dimension_type(2)='float'
      dimension_unit(2)='degree'
      dimension_long(2)='longitude'
      dimension_size(2)=n_profile
      DO i=1, n_profile
        dimension_array_fl(i, 2)=0.0_RealK
      ENDDO
!
      dimension_name(3)='plev'
      dimension_type(3)='float'
      dimension_unit(3)='pa'
      dimension_long(3)='pressure'
      dimension_size(3)=n_level
      DO i=1, n_level
        dimension_array_fl(i, 3)=p(i)
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'enter name of the output file.'
      READ(iu_stdin, '(a)') filename_out
!
!     There can only be one variable in this case.
!      WRITE(iu_stdout, '(/a)') 'enter name of field.'
!      READ(iu_stdin, '(a)') var_name(1)
!      WRITE(iu_stdout, '(/a)') 'enter unit for field.'
!      READ(iu_stdin, '(a)') var_unit(1)
!      WRITE(iu_stdout, '(/a)') 'enter long name for field.'
!      READ(iu_stdin, '(a)') var_long(1)

      n_var=1
      var_type(1)='float'
      var_unit(1)='Wm-2'
      var_size(1)=n_level
      n_dimension_var=3
      list_dimension_var(1, 1)=1
      list_dimension_var(2, 1)=2
      list_dimension_var(3, 1)=3


      var_name(1)='uflx'
      var_long(1)='Upflux'
      DO i=1, n_level
        var_fl(i, 1)=field(2, i)
      ENDDO
      CALL write_cdf(ierr
     &  , TRIM(filename_out)//'.uflx'
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, 1
     &  , n_dimension, dimension_name, dimension_type, dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , var_size, var_int, var_fl
     &  )

      var_name(1)='dflx'
      var_long(1)='Downflux'
      DO i=1, n_level
        var_fl(i, 1)=field(3, i)
      ENDDO
      CALL write_cdf(ierr
     &  , TRIM(filename_out)//'.dflx'
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, 1
     &  , n_dimension, dimension_name, dimension_type, dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , var_size, var_int, var_fl
     &  )

      var_name(1)='nflx'
      var_long(1)='Netflux'
      DO i=1, n_level
        var_fl(i, 1)=field(4, i)
      ENDDO
      CALL write_cdf(ierr
     &  , TRIM(filename_out)//'.nflx'
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, 1
     &  , n_dimension, dimension_name, dimension_type, dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , var_size, var_int, var_fl
     &  )

      var_name(1)='hrts'
      var_long(1)='Heating Rate'
      n_level=n_level-1
      dimension_size(3)=n_level
      DO i=1, n_level
        dimension_array_fl(i, 3)=field(1,i)
      ENDDO
      var_size(1)=n_level
      DO i=1, n_level
        var_fl(i, 1)=field(6, i)
      ENDDO
      CALL write_cdf(ierr
     &  , TRIM(filename_out)//'.hrts'
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, 1
     &  , n_dimension, dimension_name, dimension_type, dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , var_size, var_int, var_fl
     &  )


      END

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
program cdf2samson

  use realtype_rd
  use def_std_io_icf

  implicit none

  integer :: ierr = 0                ! Error flag
  character (LEN=80) :: samson_name, cdf_name

  integer :: nd_cdl_dimen  = 7       ! Array size for netCDF dimensions
  integer :: nd_cdl_dimen_size = 768 ! Array size for number of values
                                     !  in a netCDF-dimension
  integer :: nd_cdl_data = 768*768   ! Array size for netCDF data
  integer :: nd_cdl_var =7           ! Number of netCDF variables

  integer :: n_dimension       ! Number of netCDF dimesions
  integer :: n_var             ! Number of netCDF variables

  character :: &
       dimension_name(7)*24 ,& ! Names of dimensions
       dimension_type(7)*6  ,& ! Types of dimensions
       dimension_long(7)*40 ,& ! Long names of dimensions
       dimension_unit(7)*20 ,& ! Units of dimensions
       var_name(7)*24       ,& ! Name of variable
       var_type(7)*6        ,& ! Type of variable
       var_long(7)*40       ,& ! Long name of variable
       var_unit(7)*20          ! Unit of variable

  ! Values of integral dimensions
  integer :: dimension_array_int(768, 7)

  ! Values of floating point dimensions
  real(RealK) :: dimension_array_fl(768, 7)

  ! Sizes of dimensions
  integer :: dimension_size(7)

  ! Number of data elements
  integer :: n_data(7)

  ! Number of dimensions used by variable
  integer :: n_dimension_var(7)

  ! List of dimensions used in the variable
  integer :: list_dimension_var(7, 7)

  ! Integral data fields
  integer :: data_int(768*768, 7)

  ! Floating point data fields
  real(RealK) :: data_fl(768*768, 7)

  write(iu_stdout, '(a)') 'Enter the input netCDF filename:'
  read(iu_stdin, '(a)') cdf_name

  call read_cdf(ierr, cdf_name,                                       &
       nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,      &
       n_dimension, dimension_name, dimension_type, dimension_unit,   &
       dimension_long, dimension_size,                                &
       dimension_array_int, dimension_array_fl,                       &
       n_var, var_name, var_type, var_unit, var_long,                 &
       n_dimension_var, list_dimension_var,                           &
       n_data, data_int, data_fl )

  write(iu_stdout, '(a)') 'Enter the output SAMSON filename:'
  read(iu_stdin, '(a)') samson_name

  call write_samson(ierr, samson_name,                                &
       nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,      &
       n_dimension, dimension_name, dimension_type, dimension_unit,   &
       dimension_long, dimension_size,                                &
       dimension_array_int, dimension_array_fl,                       &
       n_var, var_name, var_type, var_unit, var_long,                 &
       n_dimension_var, list_dimension_var,                           &
       n_data, data_int, data_fl )

end program cdf2samson

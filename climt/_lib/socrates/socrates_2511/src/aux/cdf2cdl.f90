! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
program cdf2cdl

  use realtype_rd
  use def_std_io_icf
  use filenamelength_mod, ONLY: filenamelength

  implicit none

  integer :: ierr = 0                ! Error flag
  integer :: j                       ! Loop variable
  character (LEN=filenamelength) :: cdl_name, cdf_name

  ! Array size for netCDF dimensions
  integer, parameter :: nd_cdl_dimen = 7

  ! Array size for number of values in a netCDF-dimension
  integer, parameter :: nd_cdl_dimen_size = 3072

  ! Array size for netCDF data
  integer, parameter :: nd_cdl_data = 3072*3072

  ! Number of netCDF variables
  integer, parameter :: nd_cdl_var = 7

  integer :: n_dimension       ! Number of netCDF dimesions
  integer :: n_var             ! Number of netCDF variables

  character :: &
       dimension_name(nd_cdl_dimen)*24 ,& ! Names of dimensions
       dimension_type(nd_cdl_dimen)*6  ,& ! Types of dimensions
       dimension_long(nd_cdl_dimen)*40 ,& ! Long names of dimensions
       dimension_unit(nd_cdl_dimen)*20 ,& ! Units of dimensions
       var_name(nd_cdl_var)*24         ,& ! Name of variable
       var_type(nd_cdl_var)*6          ,& ! Type of variable
       var_long(nd_cdl_var)*40         ,& ! Long name of variable
       var_unit(nd_cdl_var)*20            ! Unit of variable

  ! Values of integral dimensions
  integer :: dimension_array_int(nd_cdl_dimen_size, nd_cdl_dimen)

  ! Values of floating point dimensions
  real(RealK) :: dimension_array_fl(nd_cdl_dimen_size, nd_cdl_dimen)

  ! Sizes of dimensions
  integer :: dimension_size(nd_cdl_dimen)

  ! Number of data elements
  integer :: n_data(nd_cdl_var)

  ! Number of dimensions used by variable
  integer :: n_dimension_var(nd_cdl_var)

  ! List of dimensions used in the variable
  integer :: list_dimension_var(nd_cdl_dimen, nd_cdl_var)

  ! Integral data fields
  integer :: data_int(nd_cdl_data, nd_cdl_var)

  ! Floating point data fields
  real(RealK) :: data_fl(nd_cdl_data, nd_cdl_var)

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

  ! The list of dimensions will be in C order and need to be reversed
  ! before passing to the writing routine:
  do j = 1, n_var
    list_dimension_var(1:n_dimension_var(j), j) =                     &
      list_dimension_var(n_dimension_var(j):1:-1, j)
  end do

  write(iu_stdout, '(a)') 'Enter the output CDL filename:'
  read(iu_stdin, '(a)') cdl_name

  call write_cdl(ierr, cdl_name,                                      &
       nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,      &
       n_dimension, dimension_name, dimension_type, dimension_unit,   &
       dimension_long, dimension_size,                                &
       dimension_array_int, dimension_array_fl,                       &
       n_var, var_name, var_type, var_unit, var_long,                 &
       n_dimension_var, list_dimension_var,                           &
       n_data, data_int, data_fl )

end program cdf2cdl

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
Subroutine write_samson(ierr, filename,                           &
     nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,    &
     n_dimension, dimension_name, dimension_type, dimension_unit, &
     dimension_long, dimension_size,                              &
     dimension_array_int, dimension_array_fl,                     &
     n_var, var_name, var_type, var_unit, var_long,               &
     n_dimension_var, list_dimension_var,                         &
     n_data, data_int, data_fl )

  ! Modules to set types of variables:
  Use realtype_rd
  Use error_pcf
  Use def_std_io_icf

  Implicit None

  ! Dummy arguments:
  Integer, Intent(INOUT) :: ierr            ! Error flag
  Character (LEN=*), Intent(IN) :: filename ! Name of input file
  Include 'cdf_struc.finc'

  ! Local variables:
  Integer :: i, j
  Logical :: l_exist                 ! Existence flag for file
  Integer :: iunit                   ! Unit number for output
  Integer :: ios                     ! I/O error flag
  Integer :: n_profile, n_layer      ! samson file dimensions
  Integer :: n_third                 ! third samson dimension
  Integer :: n_start                 ! Start point for line output

  ! Find a free unit and open the unit for output.
  Call get_free_unit(ierr, iunit)
  If (ierr /= i_normal) Return

  ! Check whether the file already exists
  Inquire(file=Trim(filename), exist=l_exist)
  If (l_exist) Then
     Write(iu_err, '(3a)') &
          'Error: The file "',Trim(filename), &
          '" already exists: it will not be overwritten.'
     ierr=i_err_io
     Return
  End If

  ! Open the file for output
  Open(unit=iunit, file=Trim(filename), iostat=ios)
  If (ios /= 0) Then
    Write(iu_err, '(/a, /a, /a, i3, a, i3)') &
      '*** Error: The file', Trim(filename), &
      'could not be opened on unit ', iunit, ': iostat = ', ios
    ierr=i_err_io
    Return
  End If

  ! Calculate n_profile
  n_profile = dimension_size(1)*dimension_size(2)

  ! Calculate n_layer
  If (n_dimension > 2) Then
    n_layer = dimension_size(3)
  Else
    n_layer = 1
  End If

  ! Calculate n_third
  n_third = 1
  Do i=4, n_dimension
    n_third = n_third*dimension_size(i)
  End Do

  ! Write header
  Write(iunit,*) 'SAMSON FILE ',Trim(var_name(1))
  Write(iunit,'(3i10)') n_profile, n_layer, n_third

  ! Write data
  Do j = 1, n_third
    n_start = 1 + (j-1)*n_profile*n_layer
    Do i = 1, n_layer
      Write(iunit,'(f10.1,f14.5,(32e16.6))') 0.0, dimension_array_fl(i,3), &
        data_fl(n_start:n_start+n_profile-1, 1)
      n_start = n_start + n_profile
    End Do
  End Do  

End Subroutine write_samson

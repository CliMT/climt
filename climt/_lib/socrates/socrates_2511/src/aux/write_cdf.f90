! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
Subroutine write_cdf(ierr, filename,                              &
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

  Use netcdf
  Implicit None

  ! Dummy arguments:
  Integer, Intent(INOUT) :: ierr            ! Error flag
  Character (LEN=*), Intent(IN) :: filename ! Name of input file
  Include 'cdf_struc.finc'

  ! Local variables:
  Integer :: i, j
  Logical :: l_exist                 ! Existence flag for file
  Integer :: ncid                    ! netCDF file ID
  Integer :: dimid                   ! dimension ID
  Integer :: varid                   ! variable ID
  Integer, Allocatable :: vardim(:)  ! variable dimensions

  ! Check whether the file already exists
  Inquire(file=Trim(filename), exist=l_exist)
  If (l_exist) Then
     Write(iu_err, '(3a)') &
          'Error: The file "',Trim(filename), &
          '" already exists: it will not be overwritten.'
     ierr=i_err_io
     Return
  End If

  ! Create the file and open for writing
  Call nf(nf90_create(Trim(filename),NF90_NOCLOBBER,ncid))

  ! Write dimensions
  Do i=1, n_dimension
     Call nf(nf90_def_dim(ncid, dimension_name(i), dimension_size(i), &
          dimid))
     If (dimension_type(i)(1:5) == 'float') Then
        Call nf(nf90_def_var(ncid, dimension_name(i), NF90_FLOAT, &
             dimid, varid))
        Call nf(nf90_enddef(ncid))
        Call nf(nf90_put_var(ncid, varid, &
             dimension_array_fl(1:dimension_size(i),i)))
     Else
        Call nf(nf90_def_var(ncid, dimension_name(i), NF90_INT, &
             dimid, varid))
        Call nf(nf90_enddef(ncid))
        Call nf(nf90_put_var(ncid, varid, &
             dimension_array_int(1:dimension_size(i),i)))
     End If
     Call nf(nf90_redef(ncid))
     Call nf(nf90_put_att(ncid, varid, 'title', dimension_long(i)))
     Call nf(nf90_put_att(ncid, varid, 'long_name', dimension_long(i)))
     Call nf(nf90_put_att(ncid, varid, 'units', dimension_unit(i)))
  End Do

  ! Write variables
  Do i=1,n_var
     Allocate(vardim(n_dimension_var(i)))
     vardim = (/ (dimension_size(list_dimension_var(j,i)), &
          j=1,n_dimension_var(i)) /)
     If (var_type(i)(1:5) == 'float') Then
        Call nf(nf90_def_var(ncid, var_name(i), NF90_FLOAT, &
             list_dimension_var(1:n_dimension_var(i),i), varid))
        Call nf(nf90_enddef(ncid))
        Select Case (n_dimension_var(i))
        Case (1)
           Call nf(nf90_put_var(ncid, varid, data_fl(1:n_data(i),i) ))
        Case (2)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_fl(1:n_data(i),i),(/vardim(1), vardim(2)/)) ))
        Case (3)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_fl(1:n_data(i),i),(/vardim(1), vardim(2), &
              vardim(3)/)) ))
        Case (4)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_fl(1:n_data(i),i),(/vardim(1), vardim(2), &
              vardim(3), vardim(4)/)) ))
        Case (5)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_fl(1:n_data(i),i),(/vardim(1), vardim(2), &
              vardim(3), vardim(4), vardim(5)/)) ))
        End Select
     Else
        Call nf(nf90_def_var(ncid, var_name(i), NF90_INT, &
             list_dimension_var(1:n_dimension_var(i),i), varid))
        Call nf(nf90_enddef(ncid))

        Select Case (n_dimension_var(i))
        Case (1)
           Call nf(nf90_put_var(ncid, varid, data_int(1:n_data(i),i) ))
        Case (2)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_int(1:n_data(i),i),(/vardim(1), vardim(2)/)) ))
        Case (3)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_int(1:n_data(i),i),(/vardim(1), vardim(2), &
              vardim(3)/)) ))
        Case (4)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_int(1:n_data(i),i),(/vardim(1), vardim(2), &
              vardim(3), vardim(4)/)) ))
        Case (5)
           Call nf(nf90_put_var(ncid, varid, &
              Reshape(data_int(1:n_data(i),i),(/vardim(1), vardim(2), &
              vardim(3), vardim(4), vardim(5)/)) ))
        End Select
     End If
     Deallocate(vardim)
     Call nf(nf90_redef(ncid))
     Call nf(nf90_put_att(ncid, varid, 'title', var_long(i)))
     Call nf(nf90_put_att(ncid, varid, 'long_name', var_long(i)))
     Call nf(nf90_put_att(ncid, varid, 'units', var_unit(i)))
  End Do

  Call nf(nf90_close(ncid))

Contains

  Subroutine nf(status)
    Integer, Intent(IN):: status
    If (status /= NF90_NOERR) Then
       Write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       Stop 'STOPPED!'
    End If
  End Subroutine nf

End Subroutine write_cdf

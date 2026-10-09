! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
Subroutine read_nc(ierr, filename,                                 &
     nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var,      &
     n_dimension, dimension_name, dimension_type, dimension_unit,   &
     dimension_long, dimension_size,                                &
     dimension_array_int, dimension_array_fl,                       &
     n_var, var_name, var_type, var_unit, var_long,                 &
     n_dimension_var, list_dimension_var,                           &  
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
  Integer :: status                  ! netCDF status
  Integer :: ncid                    ! netCDF file - ID
  Integer, Allocatable :: dimvid(:)  ! corresponding variable ID
  Integer :: dimtype                 ! type of corresponding coordinate var.
  Integer :: vartype                 ! type of variable
  Integer :: varid                   ! variable ID
  Integer :: nvar                    ! number of variables (incl.dims)
  Integer, Allocatable :: dimids(:)  ! dimension IDs of variable
  Real(RealK), Allocatable :: data2D(:,:)  ! real (float) data
  Real(RealK), Allocatable :: data3D(:,:,:)  ! real (float) data
  Real(RealK), Allocatable :: data4D(:,:,:,:)  ! real (float) data
  Real(RealK), Allocatable :: data5D(:,:,:,:,:)  ! real (float) data
  Integer, Allocatable :: idata2D(:,:)        ! integer data
  Integer, Allocatable :: idata3D(:,:,:)      ! integer data
  Integer, Allocatable :: idata4D(:,:,:,:)    ! integer data
  Integer, Allocatable :: idata5D(:,:,:,:,:)  ! integer data

  ! Initial setting of dimensions and variables for error checking.
  n_dimension=0
  n_var=0
  list_dimension_var=0

  ! Initialize the character strings
  Do i=1, nd_cdl_dimen
     Do j=1, Len(dimension_name(i))
        dimension_name(i)(j:j)=' '
     End Do
     Do j=1, Len(dimension_unit(i))
        dimension_unit(i)(j:j)=' '
     End Do
     Do j=1, Len(dimension_long(i))
        dimension_long(i)(j:j)=' '
     End Do
  End Do
  Do i=1, nd_cdl_var
     Do j=1, Len(var_name(i))
        var_name(i)(j:j)=' '
     End Do
     Do j=1, Len(var_unit(i))
        var_unit(i)(j:j)=' '
     End Do
     Do j=1, Len(var_long(i))
        var_long(i)(j:j)=' '
     End Do
  End Do

  ! Check whether the file exists
  Inquire(file=Trim(filename), exist=l_exist)
  If (.Not.l_exist) Then
     Write(iu_err, '(3a)') &
          'Error: The file "',Trim(filename),'" does not exist.'
     ierr=i_err_exist
     Return
  End If

  ! Open the file for reading
  Call nf(nf90_open(Trim(filename),NF90_NOWRITE,ncid))

  ! Get number of dimensions and number of variables
  Call nf(nf90_Inquire(ncid, nDimensions=n_dimension, nVariables=nvar))

  ! Get dimension lengths, names, id and type of corresponding variable
  Allocate(dimvid(n_dimension))
  Do i=1, n_dimension

     Call nf(nf90_Inquire_Dimension(ncid, i, name=dimension_name(i), &
          len=dimension_size(i)))
     Call nf(nf90_inq_varid(ncid, Trim(dimension_name(i)), varid))
     dimvid(i) = varid

     status = nf90_get_att(ncid, dimvid(i), 'title', &
          dimension_long(i))
     If (status /= NF90_NOERR) Then
        status = nf90_get_att(ncid, dimvid(i), 'long_name', &
             dimension_long(i))
        If (status /= NF90_NOERR) Then
           Write(iu_err, '(3a)') 'Dimension "', &
                Trim(dimension_name(i)),'" has no title or long_name'
        End If
     End If

     status = nf90_get_att(ncid, dimvid(i), 'units', &
          dimension_unit(i))
     If (status /= NF90_NOERR) Then
        Write(iu_err, '(3a)') &
             'Dimension "',Trim(dimension_name(i)),'" has no units'
     End If

     Call nf(nf90_Inquire_Variable(ncid, varid, xtype=dimtype))

     Select Case (dimtype)
     Case (NF90_FLOAT,NF90_DOUBLE,0)
        dimension_type(i)='float'
        Call nf(nf90_get_var(ncid, dimvid(i), &
             dimension_array_fl(1:dimension_size(i),i)))
     Case (NF90_SHORT, NF90_INT, NF90_BYTE)
        dimension_type(i)='int'
        Call nf(nf90_get_var(ncid, dimvid(i), &
             dimension_array_int(1:dimension_size(i),i)))
     Case DEFAULT
     End Select
  End Do

  ! Get variable name and dimensions
  j=0
  Do i=1, nvar
     If (All(dimvid /= i)) Then
        j=j+1

        status = nf90_get_att(ncid, i, 'title', var_long(j))
        If (status /= NF90_NOERR) Then
           status = nf90_get_att(ncid, i, 'long_name', var_long(j))
           If (status /= NF90_NOERR) Then
              Write(iu_err, '(3a)') 'Variable "', &
                   Trim(var_name(j)),'" has no title or long_name'
           End If
        End If

        status = nf90_get_att(ncid, i, 'units', var_unit(j))
        If (status /= NF90_NOERR) Then
           Write(iu_err, '(3a)') 'Variable "', &
                Trim(var_name(j)),'" has no units'
        End If

        Call nf(nf90_Inquire_Variable(ncid, i, &
             ndims=n_dimension_var(j), xtype=vartype))
        Allocate(dimids(n_dimension_var(j)))
        Call nf(nf90_Inquire_Variable(ncid, i, &
             name=var_name(j), dimids=dimids))
        n_data(j)= Product(dimension_size(dimids))
        list_dimension_var(1:n_dimension_var(j),j) &
             = dimids(n_dimension_var(j):1:-1)

        Select Case (vartype)
        Case (NF90_FLOAT,NF90_DOUBLE,0)
           var_type(j)='float'
           Select Case (n_dimension_var(j))
           Case (1)
              Call nf(nf90_get_var(ncid, i, data_fl(1:n_data(j),j) ))
           Case (2)
              Allocate(data2D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2))  ))
              Call nf(nf90_get_var(ncid, i, data2D))
              data_fl(1:n_data(j),j)=Reshape(data2D,(/n_data(j)/))
              Deallocate(data2D)
           Case (3)
              Allocate(data3D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2)), &
                   dimension_size(dimids(3))  ))
              Call nf(nf90_get_var(ncid, i, data3D))
              data_fl(1:n_data(j),j)=Reshape(data3D,(/n_data(j)/))
              Deallocate(data3D)
           Case (4)
              Allocate(data4D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2)), &
                   dimension_size(dimids(3)), &
                   dimension_size(dimids(4))  ))
              Call nf(nf90_get_var(ncid, i, data4D))
              data_fl(1:n_data(j),j)=Reshape(data4D,(/n_data(j)/))
              Deallocate(data4D)
           Case (5)
              Allocate(data5D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2)), &
                   dimension_size(dimids(3)), &
                   dimension_size(dimids(4)), &
                   dimension_size(dimids(5))  ))
              Call nf(nf90_get_var(ncid, i, data5D))
              data_fl(1:n_data(j),j)=Reshape(data5D,(/n_data(j)/))
              Deallocate(data5D)
           Case DEFAULT
           End Select
        Case (NF90_SHORT, NF90_INT, NF90_BYTE)
           var_type(j)='int'
           Select Case (n_dimension_var(j))
           Case (1)
              Call nf(nf90_get_var(ncid, i, data_int(1:n_data(j),j) ))
           Case (2)
              Allocate(idata2D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2))  ))
              Call nf(nf90_get_var(ncid, i, idata2D))
              data_int(1:n_data(j),j)=Reshape(idata2D,(/n_data(j)/))
              Deallocate(idata2D)
           Case (3)
              Allocate(idata3D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2)), &
                   dimension_size(dimids(3))  ))
              Call nf(nf90_get_var(ncid, i, idata3D))
              data_int(1:n_data(j),j)=Reshape(idata3D,(/n_data(j)/))
              Deallocate(idata3D)
           Case (4)
              Allocate(idata4D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2)), &
                   dimension_size(dimids(3)), &
                   dimension_size(dimids(4))  ))
              Call nf(nf90_get_var(ncid, i, idata4D))
              data_int(1:n_data(j),j)=Reshape(idata4D,(/n_data(j)/))
              Deallocate(idata4D)
           Case (5)
              Allocate(idata5D(dimension_size(dimids(1)), &
                   dimension_size(dimids(2)), &
                   dimension_size(dimids(3)), &
                   dimension_size(dimids(4)), &
                   dimension_size(dimids(5))  ))
              Call nf(nf90_get_var(ncid, i, idata5D))
              data_int(1:n_data(j),j)=Reshape(idata5D,(/n_data(j)/))
              Deallocate(idata5D)
           Case DEFAULT
           End Select
        Case DEFAULT
        End Select
        Deallocate(dimids)
     End If
  End Do

  Deallocate(dimvid)
  n_var = nvar - n_dimension

  Call nf(nf90_close(ncid))

!!$  Print *,'' '
!!$  Print *,''Number of dimensions:',n_dimension
!!$  Print *,''---------------------------------'
!!$  Do i=1, n_dimension
!!$     Print *,'' '
!!$     Print *,''Name:      ',dimension_name(i)
!!$     Print *,''Type:      ',dimension_type(i)
!!$     Print *,''Long name: ',dimension_long(i)
!!$     Print *,''Units:     ',dimension_unit(i)
!!$     Print *,''Number of values:',dimension_size(i)
!!$     Print *,''Values:'
!!$     If (dimension_type(i) == ''float') Then
!!$        Print *,dimension_array_fl(1:dimension_size(i),i)
!!$     Else
!!$        Print *,dimension_array_int(1:dimension_size(i),i)
!!$     End If
!!$  End Do
!!$
!!$  Print *,'' '
!!$  Print *,''Number of variables:',n_var
!!$  Print *,''---------------------------------'
!!$  Do j=1, n_var
!!$     Print *,'' '
!!$     Print *,''Name:      ',var_name(j)
!!$     Print *,''Type:      ',var_type(j)
!!$     Print *,''Long name: ',var_long(j)
!!$     Print *,''Units:     ',var_unit(j)
!!$     Print *,''Number of values:',n_data(j)
!!$     Print *,''Number of dimensions used:',n_dimension_var(j)
!!$     Print *,''List of dimensions used:', &
!!$          list_dimension_var(1:n_dimension_var(j),j)
!!$     Print *,''Values:'
!!$     If (var_type(j) == ''float') Then
!!$        Print *,data_fl(1:n_data(j),j)
!!$     Else
!!$        Print *,data_int(1:n_data(j),j)
!!$     End If
!!$  End Do

Contains

  Subroutine nf(status)
    Integer, Intent(IN):: status
    If (status /= NF90_NOERR) Then
       Write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       Stop 'STOPPED!'
    End If
  End Subroutine nf

End Subroutine read_nc

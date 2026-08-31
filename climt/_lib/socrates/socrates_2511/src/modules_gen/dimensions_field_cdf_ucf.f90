! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module setting the dimensions of physical fields

MODULE dimensions_field_cdf_ucf

  IMPLICIT NONE

  INTEGER :: npd_latitude = 131
!       Number of latitudes
  INTEGER :: npd_longitude = 131
!       Number of longitudes
  INTEGER :: npd_profile = 131
!       Number of atmospheric profiles
  INTEGER :: npd_layer = 170
!       Number of atmospheric layers
  INTEGER :: npd_channel = 2
!       Number of spectral channels permitted for output
  INTEGER :: npd_column = 24
!       Maximum number of cloudy subcolumns
  INTEGER :: npd_direction = 63
!       Maximum number of directions for radiances
  INTEGER :: npd_max_order = 101
!       Maximum order of spherical harmonics used
  INTEGER :: npd_brdf_basis_fnc = 2
!       Number of BRDF basis functions
  INTEGER :: npd_brdf_trunc = 5
!       Order of BRDF truncation
  INTEGER :: npd_profile_aerosol_prsc = 9
!       Size allocated for profiles of prescribed
!       cloudy optical properties
  INTEGER :: npd_profile_cloud_prsc = 9
!       Size allocated for profiles of prescribed
!       aerosol optical properties
  INTEGER :: npd_opt_level_aerosol_prsc = 170
!       Size allocated for levels of prescribed
!       cloudy optical properties
  INTEGER :: npd_opt_level_cloud_prsc = 170
!       Size allocated for levels of prescribed
!       aerosol optical properties

Contains

  Subroutine set_dimensions_field(ierr, filename)

    Use error_pcf
    Use def_std_io_icf

    Use netcdf
    Implicit None

    ! Dummy arguments:
    Integer, Intent(INOUT) :: ierr            ! Error flag
    Character (LEN=*), Intent(IN) :: filename ! Name of input file

    ! Local variables:
    Integer :: i
    Logical :: l_exist                 ! Existence flag for file
    Integer :: status                  ! netCDF status
    Integer :: ncid                    ! netCDF file - ID
    Integer :: ndim                    ! number of dimensions
    Character (LEN=16) :: dimname      ! name of dimension
    Integer :: dimlen                  ! size of dimension

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

    ! Get number of dimensions
    Call nf(nf90_Inquire(ncid, nDimensions=ndim))

    ! Get dimension lengths, names, id and type of corresponding variable
    Do i=1, ndim
       Call nf(nf90_Inquire_Dimension(ncid, i, name=dimname, &
            len=dimlen))
       Select Case (Trim(dimname))
       Case ('lon')
          npd_longitude = 2*(dimlen/2)+1
          print *,'Setting npd_longitude: ',npd_longitude
       Case ('lat')
          npd_latitude = 2*(dimlen/2)+1
          print *,'Setting npd_latitude: ',npd_latitude
       Case ('plev')
          npd_layer = 2*(dimlen/2)+1
          print *,'Setting npd_layer: ',npd_layer
       Case ('direction')
          npd_direction = 2*(dimlen/2)+1
          print *,'Setting npd_direction: ',npd_direction
       Case DEFAULT
       End Select
    End Do

    npd_profile = npd_latitude*npd_longitude
    print *,'Setting npd_profile: ',npd_profile

  Contains

    Subroutine nf(status)
      Integer, Intent(IN):: status
      If (status /= NF90_NOERR) Then
         Write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
         Stop 'STOPPED!'
      End If
    End Subroutine nf

  End Subroutine set_dimensions_field

END MODULE dimensions_field_cdf_ucf

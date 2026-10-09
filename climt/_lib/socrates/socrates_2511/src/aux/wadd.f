! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to construct a weighted average of two CDL-files.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      PROGRAM wadd
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
      CHARACTER
     &    file_1*80
!           Name of first input file
     &  , file_2*80
!           Name of second input file
     &  , file_out*80
!           Output file
     &  , name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile_1
!           Number of profiles in first file
     &  , n_level_1
!           Number of levels in first file
     &  , n_profile_2
!           Number of profiles in second file
     &  , n_level_2
!           Number of levels in second file
     &  , i
!           Loop variable
     &  , l
!           Loop variable
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , p_1(npd_profile, npd_layer+1)
!           Pressures in first file
     &  , p_2(npd_profile, npd_layer+1)
!           Pressures in second file
     &  , field_1(npd_profile, npd_layer+1)
!           First field read in
     &  , field_2(npd_profile, npd_layer+1)
!           Second field read in
     &  , field_out(npd_profile, npd_layer+1)
!           Output field
     &  , weight_1
!           First weighting factor
     &  , weight_2
!           First weighting factor
!
!     Numerical precision:
      REAL  (RealK) ::
     &    tol_p
!           Tolerance for equality of pressures
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, output_vert_cdl
!
!
      data
     &    l_vert_coord/.false./
      data ierr/i_normal/
      data n_latitude/0/
      data n_longitude/0/
!
!
!
!
!     Set the numerical precision
      tol_p=1.0e+03_RealK*epsilon(tol_p)
!
      WRITE(iu_stdout, '(/a)') 'Enter the name of 
     &                          the first input file.'
      READ(iu_stdin, '(a)') file_1
      CALL assign_input_vert_cdl(ierr
     &  , file_1, 'first field', l_vert_coord, name_vert_coord
     &  , .true., n_level_1, .true.
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile_1, n_level_1
     &  , p_1, field_1
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
1     WRITE(iu_stdout, '(/a)') 
     &  'Enter the weighting factor for the first file.'
      READ(iu_stdin, *, iostat=ios) weight_1
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Illegal response: please re-enter.'
        goto 1
      ENDIF
!
      WRITE(iu_stdout, '(/a)') 'Enter the name of the 
     &                                      second input file.'
      READ(iu_stdin, '(a)') file_2
      CALL assign_input_vert_cdl(ierr
     &  , file_2, 'first field', l_vert_coord, name_vert_coord
     &  , .true., n_level_1, .true.
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile_2, n_level_2
     &  , p_2, field_2
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
2     WRITE(iu_stdout, '(/a)') 
     &  'Enter the weighting factor for the second file.'
      READ(iu_stdin, *, iostat=ios) weight_2
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Illegal response: please re-enter.'
        goto 2
      ENDIF
!
!     Check for compatibility.
      IF ( (n_profile_1 /= n_profile_2).OR.
     &     (n_level_1 /= n_level_2) ) THEN
        WRITE(iu_err, '(/a)')
     &    '*** Error: these files are of incompatible sizes.'
        STOP
      ENDIF
      DO i=1, n_level_1
        IF ( abs(p_1(1, i)-p_2(1, i)) < tol_p ) THEN
          DO l=1, n_profile_1
            field_out(l, i)=weight_1*field_1(l, i)
     &        +weight_2*field_2(l, i)
          ENDDO
        ELSE
          WRITE(iu_err, '(/a)')
     &      '*** Error: These files have incompatible pressure levels.'
          STOP
        ENDIF
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'Enter the name of the output file.'
      READ(iu_stdin, '(a)') file_out
      CALL output_vert_cdl(ierr
     &  , file_out
     &  , n_latitude, latitude, n_longitude, longitude, n_profile_1
     &  , n_level_1, trim(name_vert_coord), len(trim(name_vert_coord))
     &  , p_1
     &  , 'weighted sum', 'float', 'unknown', 'weighted sum', field_out
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      STOP
      END

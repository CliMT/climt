! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to scale a CDL-field.
!
! Method:
!   A CDL-field is read in. A scaling factor and a range of
!   pressures are obtained. Within the range of pressures the
!   field is multiplied by the scaling factor and written out
!   with a new title.
!
!- ---------------------------------------------------------------------
      PROGRAM scale_field
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      CHARACTER
     &    file_in*80
!           Name of input file
     &  , char_yn*1
!           Reponse variable
     &  , file_out*80
!           Output file
     &  , name_out*24
!           Name of the output field
     &  , unit_out*20
!           Units of the output field
     &  , long_out*40
!           Long name of the output field
     &  , name_vert_coord*24
!           Name of vertical coordinate
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
     &  , i
!           Loop variable
     &  , l
!           Loop variable
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , p(npd_profile, npd_layer+1)
!           Pressure levels
     &  , field(npd_profile, npd_layer+1)
!           Field read in
     &  , p_low
!           Low scaling pressure
     &  , p_high
!           High scaling pressure
     &  , scale_factor
!           Scaling factor
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, output_vert_cdl
!
      data
     &    l_vert_coord/.false./
      data ierr/i_normal/
      data n_latitude/0/
      data n_longitude/0/
!
!
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the file for scaling.'
      READ(iu_stdin, '(a)') file_in
      CALL assign_input_vert_cdl(ierr
     &  , file_in, 'input field', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, field
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
1     write(iu_stdout, '(/a)') 
     &  'enter the range of pressures for scaling.'
      READ(iu_stdin, *, iostat=ios) p_low, p_high
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ illegal response: please re-enter.'
        goto 1
      ENDIF
2     write(iu_stdout, '(/a)') 'enter the scaling factor.'
      READ(iu_stdin, *, iostat=ios) scale_factor
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ illegal response: please re-enter.'
        goto 2
      ENDIF
      DO i=1, n_level
        DO l=1, n_profile
          IF ( (p(l, i) >= p_low).AND.(p(l, i) <= p_high) ) THEN
            field(l, i)=scale_factor*field(l, i)
          ENDIF
        ENDDO
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'scale a new range? (y/n)'
      READ(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) goto 1
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_out
      WRITE(iu_stdout, '(/a)') 'enter the name of the output field.'
      READ(iu_stdin, '(a)') name_out
      WRITE(iu_stdout, '(/a)') 'enter the units of the output field.'
      READ(iu_stdin, '(a)') unit_out
      WRITE(iu_stdout, '(/a)') 'enter the long-name 
     &                                   of the output field.'
      READ(iu_stdin, '(a)') long_out
!
      CALL output_vert_cdl(ierr
     &  , file_out
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level, trim(name_vert_coord), len(trim(name_vert_coord)), p
     &  , name_out, 'float', unit_out, long_out, field
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      STOP
      END

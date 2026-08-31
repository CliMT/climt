! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to form a CDL-file of mid-points of layers.
!
! Method:
!   A CDL-file is read in, giving the boundaries of layers.
!   A CDL-file with levels at the mid-points of these layers 
!   is written out with a zero field.
!
!- ---------------------------------------------------------------------
      PROGRAM mid_point
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
      INTEGER
     &    ierr
!           Error flag
!
      CHARACTER
     &    file_in*80
!           Name of input file
     &  , file_out*80
!           Name of output file
     &  , name_vert_coord*24
!           Name of vertical coordinate
!
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
!
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of input profile
     &  , n_level_in
!           Number of input levels
     &  , n_level_out
!           Number of levels in the output file
     &  , i
!           Loop variable
     &  , l
!           Loop variable
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , field_in(npd_profile, npd_layer+1)
!           Input field
     &  , field_out(npd_profile, npd_layer)
!           Output field
     &  , p_in(npd_profile, npd_layer+1)
!           Pressure for interpolation
     &  , p_out(npd_profile, npd_layer)
!           Pressure levels of output field
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
!     Read in the basic field.
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the input file.'
      READ(iu_stdin, '(a)') file_in
      CALL assign_input_vert_cdl(ierr
     &  , file_in, 'input field', l_vert_coord, name_vert_coord
     &  , .true., n_level_in, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level_in
     &  , p_in, field_in
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!     Read in the name of the output file.
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_out
!
!     The output field will have one layer fewer than the input.
      n_level_out=n_level_in-1
      DO l=1, n_profile
        DO i=1, n_level_out
          p_out(l, i)=0.5_RealK*(p_in(l, i+1)+p_in(l, i))
        ENDDO
        DO i=1, n_level_out
          field_out(l, i)=0.0_RealK
        ENDDO
      ENDDO
!
      CALL output_vert_cdl(ierr
     &  , file_out
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level_out, trim(name_vert_coord), len(trim(name_vert_coord))
     &  , p_out, 'null', 'float', 'none', 'null field', field_out
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!
!
      STOP
      END

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert from humidity mixing ratios to specific humidity.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      PROGRAM hmrsph
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
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      CHARACTER
     &    file_hmr*80
!           File of humidity mixing ratios
     &  , file_sph*80
!           File of specific humidities
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
     &  , q(npd_profile, npd_layer+1)
!           Specific humidity
     &  , r(npd_profile, npd_layer+1)
!           Humidity mixing ratio
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
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the file containing '
     &  //'the humidity mixing ratios.'
      READ(iu_stdin, '(a)') file_hmr
      CALL assign_input_vert_cdl(ierr
     &  , file_hmr, 'mixing ratios', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, r
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
      DO i=1, n_level
        DO l=1, n_profile
          q(l, i)=r(l, i)/(1.0_RealK+r(l, i))
        ENDDO
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_sph
      CALL output_vert_cdl(ierr
     &  , file_sph
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level, trim(name_vert_coord), len(trim(name_vert_coord)), p
     &  , 'specific humidity', 'float', 'none', 'specific humidity', q
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      STOP
      END

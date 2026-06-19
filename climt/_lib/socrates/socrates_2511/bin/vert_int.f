! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to integrate mixing ratios vertically through the atmosphere.
!
! Method:
!   A CDL-file is read in. Its contents are assumed to be mixing
!   ratios on pressure levels. These pressure levels give the
!   masses in the layers so defined. The mean mixing ratio in the
!   layer is taken as the average of the values at the boundaries.
!   This enables us to calculate the amount of the quantity
!   integrating downwards. The results are written to a CDL-file.
!
!- ---------------------------------------------------------------------
      PROGRAM vert_int
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE error_pcf
      USE rad_ccf, ONLY: grav_acc
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
      CHARACTER
     &    file_mix*80
!           Name of file with mixing ratio
     &  , file_out*80
!           Name of file for output
     &  , name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      INTEGER
     &    ierr
!           Error flag
     &  , n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of levels
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
     &  , mix_ratio(npd_profile, npd_layer+1)
!           Mixing ratio
     &  , amount(npd_profile, npd_layer+1)
!           Amounts of quantity
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
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the file containing the mixing ratios.'
      READ(iu_stdin, '(a)') file_mix
      CALL assign_input_vert_cdl(ierr
     &  , file_mix, 'mixing ratios', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, mix_ratio
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
      DO l=1, n_profile
        amount(l, 1)=0.0_RealK
      ENDDO
      DO i=2, n_level
        DO l=1, n_profile
          amount(l, i)=amount(l, i-1)
     &      +0.5_RealK*(mix_ratio(l, i-1)+mix_ratio(l, i))
     &      *(p(l, i)-p(l, i-1))/grav_acc
        ENDDO
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_out
      CALL output_vert_cdl(ierr
     &  , file_out
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level, trim(name_vert_coord), len(trim(name_vert_coord)), p
     &  , 'vint', 'float', 'm-2'
     &  , 'vertically integrated field', amount
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      STOP
      END

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to calculate heating rates from net fluxes.
!
! Method:
!   From the net fluxes, the net flow of energy into a layer is
!   determined. From the pressures the mass within the layer
!   can be found. The heating rates are now calculated assuming
!   that the specific heat may be taken as that of dry air.
!   The heating rates are allocated pressurea at the mid-points
!   of layers.
!
!- ---------------------------------------------------------------------
      PROGRAM nettohr
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE error_pcf
      USE rad_ccf, ONLY: grav_acc, cp_air_dry, seconds_per_day
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
     &    file_net*80
!           Name of file of net fluxes
     &  , file_hr*80
!           Name of output file
     &  , name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
     &  , n_layer
!           Number of layers
     &  , i
!           Loop variable
!
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , p(npd_profile, 0: npd_layer)
!           Pressure levels
     &  , net_flux(npd_profile, 0: npd_layer)
!           Net fluxes
     &  , p_layer(npd_profile, npd_layer)
!           Pressures in layers
     &  , heat_rate(npd_profile, npd_layer)
!           Heating rates
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, output_vert_cdl
!
!
      data 
     &   l_vert_coord/.false./
      data ierr/i_normal/
      data n_latitude/0/
      data n_longitude/0/
!
!
!
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the file containing the net fluxes.'
      READ(iu_stdin, '(a)') file_net
      CALL assign_input_vert_cdl(ierr
     &  , file_net, 'net fluxes', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 0
     &  , n_profile, n_level
     &  , p, net_flux
     &  , npd_profile, npd_latitude, npd_longitude, 0, npd_layer
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      n_layer=n_level-1
      IF (ierr /= i_normal) STOP
!
      DO i=1, n_layer
        p_layer(1, i)=0.5_RealK*(p(1, i-1)+p(1, i))
        heat_rate(1, i)=-grav_acc*(net_flux(1, i)-net_flux(1, i-1))
     &    /(cp_air_dry*(p(1, i)-p(1, i-1)))*seconds_per_day
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the file '
     &   //'to contain the heating rates.'
      READ(iu_stdin, '(a)') file_hr
      CALL output_vert_cdl(ierr
     &  , file_hr
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_layer, trim(name_vert_coord), len(trim(name_vert_coord))
     &  , p_layer
     &  , 'heatrate', 'float', 'k.day-1', 'heating rate', heat_rate
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      STOP
      END

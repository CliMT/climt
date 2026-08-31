! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to control the writing of CDL-files of fluxes.
!
! Purpose:
!   This subroutine receives fluxes and heating rates as input and
!   calls routines to write them to CDL-files.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE output_flux_cdl(ierr
     &  , control, sp
     &  , base_name, length_name
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_layer, name_vert_coord, len_vert_coord, p, p_level
     &  , n_channel
     &  , flux_down, flux_down_diffuse, flux_up, flux_net, flux_direct
     &  , heating_rate, contrib_funci, contrib_funcf
     &  , nd_profile, nd_latitude, nd_longitude, nd_layer, nd_channel
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE def_control,  ONLY: StrCtrl
      USE def_spectrum, ONLY: StrSpecData
      USE realtype_rd
      USE def_std_io_icf
      USE rad_pcf
      USE gas_list_pcf
      USE input_head_pcf
!
!
      IMPLICIT NONE


!     Declaration of variables:
!
      INTEGER, INTENT(INOUT) :: ierr
!            Error flag

!     Control options:
      TYPE(StrCtrl),     INTENT(IN) :: control

!     Spectral information:
      TYPE(StrSpecData), INTENT(IN) :: sp

      INCLUDE 'cdl_struc.finc'

!     Sizes of arrays
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Allowed size for profiles
     &  , nd_latitude
!           Allowed size for latitudes
     &  , nd_longitude
!           Allowed size for longitudes
     &  , nd_layer
!           Allowed size for layers
     &  , nd_channel
!           Allowed size for spectral channels
!
      INTEGER, Intent(IN) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_layer
!           Number of layers
     &  , n_channel
!           Number of channels used
!
      REAL  (RealK), Intent(IN) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
!
      CHARACTER !, Intent(IN)
     &    base_name*(*)
!           Base name of input file
     &  , name_vert_coord*(*)
!           Name of vertical coordinate
      INTEGER, Intent(IN) ::
     &    length_name
!           Length of basename
     &  , len_vert_coord
!           Length of name of vertical coordinate
!
      REAL  (RealK), Intent(IN) ::
     &    p(nd_profile, nd_layer)
!           Atmsopheric pressures at the middle of layers
     &  , p_level(nd_profile, 0: nd_layer)
!           Pressures on levels
      REAL  (RealK), Intent(IN) ::
     &    flux_down(nd_profile, 0: nd_layer, nd_channel)
!           Total downward fluxes
     &  , flux_down_diffuse(nd_profile, 0: nd_layer, nd_channel)
!           Diffuse downward fluxes
     &  , flux_up(nd_profile, 0: nd_layer, nd_channel)
!           Upward fluxes
     &  , flux_net(nd_profile, 0: nd_layer, nd_channel)
!           Net fluxes
     &  , flux_direct(nd_profile, 0: nd_layer, nd_channel)
!           Direct fluxes
     &  , heating_rate(nd_profile, nd_layer, nd_channel)
!           Heating rates in layers
     &  , contrib_funci(nd_profile, nd_layer, nd_channel)
!           Contribution function (intensity)
     &  , contrib_funcf(nd_profile, nd_layer, nd_channel)
!           Contribution function (flux)
!
!
!     Local Variables
!
      INTEGER
     &    i
!           Loop variable
     &  , k
!           Loop variable
     &  , l
!           Loop variable
     &  , point
!           Point addressed in CDL-array
!
      CHARACTER !, Intent(IN)
     &    file_name*80
!           Name of input file
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     write_cdl
!
!
!
!
!     Set up common features of dimensions of fluxes.
      n_dimension=3
      dimension_name(1)='lon'
      dimension_name(2)='lat'
      dimension_type(1)='float'
      dimension_type(2)='float'
      dimension_unit(1)='degree'
      dimension_unit(2)='degree'
      dimension_long(1)='longitude'
      dimension_long(2)='latitude'
      dimension_size(1)=n_longitude
      dimension_size(2)=n_latitude
      DO l=1, n_longitude
        dimension_array_fl(l, 1)=longitude(l)
      ENDDO
      DO l=1, n_latitude
        dimension_array_fl(l, 2)=latitude(l)
      ENDDO
!
!     If multispectral output is to be generated add a dimension
!     for channels.
      IF (n_channel > 1) THEN
        n_dimension=4
        dimension_name(4)='channel'
        dimension_type(4)='int'
        dimension_unit(4)='none'
        dimension_long(4)='spectral channel'
        dimension_size(4)=n_channel
        DO k=1, n_channel
          dimension_array_int(k, 4)=k
        ENDDO
      ENDIF
!
!     Vertical coordinate:
      IF (name_vert_coord(1:len_vert_coord) == 'plev') THEN
        dimension_name(3)='plev'
        dimension_type(3)='float'
        dimension_unit(3)='pa'
        dimension_long(3)='pressure'
        dimension_size(3)=n_layer+1
        DO i=1, n_layer+1
          dimension_array_fl(i, 3)=p_level(1, i-1)
        ENDDO
      ELSE IF (name_vert_coord(1:len_vert_coord) == 'level') THEN
        dimension_name(3)='level'
        dimension_type(3)='int'
        dimension_unit(3)='none'
        dimension_long(3)='level'
        dimension_size(3)=n_layer+1
        DO i=1, n_layer+1
          dimension_array_int(i, 3)=i-1
        ENDDO
      ELSE
        WRITE(iu_err, '(/a)')
     &    '*** Error: An illegal vertical coordinate was specified'
     &    , 'for the calculated fluxes.'
      ENDIF
!
!     Set common properties of the variables.
      IF (n_channel == 1) THEN
        n_var=1
      ELSE
        n_var=4
        n_dimension_var(1)=1
        list_dimension_var(1, 1)=4
        var_name(1)='wl_short'
        var_type(1)='float'
        var_unit(1)='m'
        var_long(1)='Wavelength lower bound'
        n_data(1)=n_channel
        data_fl(1:n_channel, 1)=MAXVAL(sp%basic%wavelength_short
     &    (control%first_band:control%last_band))
        DO i=control%first_band, control%last_band
           data_fl(control%map_channel(i), 1) =
     &       MIN( data_fl(control%map_channel(i), 1),
     &            sp%basic%wavelength_short(i) )
        END DO
        n_dimension_var(2)=1
        list_dimension_var(1, 2)=4
        var_name(2)='wl_long'
        var_type(2)='float'
        var_unit(2)='m'
        var_long(2)='Wavelength upper bound'
        n_data(2)=n_channel
        data_fl(1:n_channel, 2)=MINVAL(sp%basic%wavelength_long
     &    (control%first_band:control%last_band))
        DO i=control%first_band, control%last_band
           data_fl(control%map_channel(i), 2) =
     &       MAX( data_fl(control%map_channel(i), 2),
     &            sp%basic%wavelength_long(i) )
        END DO
        n_dimension_var(3)=1
        list_dimension_var(1, 3)=4
        var_name(3)='bandwidth'
        var_type(3)='float'
        var_unit(3)='m'
        var_long(3)='Channel width'
        n_data(3)=n_channel
        data_fl(1:n_channel, 3)=0.0_RealK
        DO i=control%first_band, control%last_band
          data_fl(control%map_channel(i), 3) =
     &      data_fl(control%map_channel(i), 3) +
     &      sp%basic%wavelength_long(i) - sp%basic%wavelength_short(i)
          DO l=1, sp%basic%n_band_exclude(i)
            k = Sp%Basic%index_exclude(l, i)
            data_fl(control%map_channel(i), 3) =
     &        data_fl(control%map_channel(i), 3) -
     &        sp%basic%wavelength_long(k) + sp%basic%wavelength_short(k)
          END DO
        END DO
      END IF  
      IF (trim(name_vert_coord) == 'level') THEN
        n_dimension_var(n_var)=3
        list_dimension_var(1, n_var)=1
        list_dimension_var(2, n_var)=2
        list_dimension_var(3, n_var)=3
        var_name(n_var)='plev'
        var_type(n_var)='float'
        var_unit(n_var)='pa'
        var_long(n_var)='pressure levels'
        n_data(n_var)=n_profile*(n_layer+1)
        DO i=1, n_layer+1
          DO l=1, n_profile
            point=l+(i-1)*n_profile
            data_fl(point, n_var)=p_level(l, i-1)
          ENDDO
        ENDDO
        n_var=n_var+1
      ENDIF
!
      IF (n_channel == 1) THEN
        n_dimension_var(n_var)=3
        list_dimension_var(1, n_var)=1
        list_dimension_var(2, n_var)=2
        list_dimension_var(3, n_var)=3
      ELSE
        n_dimension_var(n_var)=4
        list_dimension_var(1, n_var)=1
        list_dimension_var(2, n_var)=2
        list_dimension_var(3, n_var)=3
        list_dimension_var(4, n_var)=4
      ENDIF
!
!
!     Downward Flux:
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_flux_down)
      var_name(n_var)='dflx'
      var_type(n_var)='float'
      var_unit(n_var)='wm-2'
      IF (control%isolir == IP_solar) THEN
        var_long(n_var)='diffuse downward flux'
        n_data(n_var)=n_profile*(n_layer+1)*n_channel
        DO k=1, n_channel
          DO i=1, n_layer+1
            DO l=1, n_profile
              point=l+(i-1+(k-1)*(n_layer+1))*n_profile
              data_fl(point, n_var)=flux_down_diffuse(l, i-1, k)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        var_long(n_var)='downward flux'
        n_data(n_var)=n_profile*(n_layer+1)*n_channel
        DO k=1, n_channel
          DO i=1, n_layer+1
            DO l=1, n_profile
              point=l+(i-1+(k-1)*(n_layer+1))*n_profile
              data_fl(point, n_var)=flux_down(l, i-1, k)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      CALL write_cdl(ierr
     &  , file_name(1: length_name+1+len_file_suffix)
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Upward Flux:
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_flux_up)
      var_name(n_var)='uflx'
      var_type(n_var)='float'
      var_unit(n_var)='wm-2'
      var_long(n_var)='upward flux'
      n_data(n_var)=n_profile*(n_layer+1)*n_channel
      DO k=1, n_channel
        DO i=1, n_layer+1
          DO l=1, n_profile
            point=l+(i-1+(k-1)*(n_layer+1))*n_profile
            data_fl(point, n_var)=flux_up(l, i-1, k)
          ENDDO
        ENDDO
      ENDDO
      CALL write_cdl(ierr
     &  , file_name(1: length_name+1+len_file_suffix)
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Net Flux:
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_flux_net)
      var_name(n_var)='nflx'
      var_type(n_var)='float'
      var_unit(n_var)='wm-2'
      var_long(n_var)='net downward flux'
      n_data(n_var)=n_profile*(n_layer+1)*n_channel
      DO k=1, n_channel
        DO i=1, n_layer+1
          DO l=1, n_profile
            point=l+(i-1+(k-1)*(n_layer+1))*n_profile
            data_fl(point, n_var)=flux_net(l, i-1, k)
          ENDDO
        ENDDO
      ENDDO
      CALL write_cdl(ierr
     &  , file_name(1: length_name+1+len_file_suffix)
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
      IF (ierr /= i_normal) RETURN
!
      IF (control%isolir == IP_solar) THEN
!
!       Direct Flux:
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)//'.'//phys_suffix(IP_flux_direct)
        var_name(n_var)='sflx'
        var_type(n_var)='float'
        var_unit(n_var)='wm-2'
        var_long(n_var)='direct flux'
        n_data(n_var)=n_profile*(n_layer+1)*n_channel
        DO k=1, n_channel
          DO i=1, n_layer+1
            DO l=1, n_profile
              point=l+(i-1+(k-1)*(n_layer+1))*n_profile
              data_fl(point, n_var)=flux_direct(l, i-1, k)
            ENDDO
          ENDDO
        ENDDO
        CALL write_cdl(ierr
     &    , file_name(1: length_name+1+len_file_suffix)
     &    , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &    , n_dimension, dimension_name, dimension_type
     &    , dimension_unit
     &    , dimension_long, dimension_size
     &    , dimension_array_int, dimension_array_fl
     &    , n_var, var_name, var_type, var_unit, var_long
     &    , n_dimension_var, list_dimension_var
     &    , n_data, data_int, data_fl
     &    )
        IF (ierr /= i_normal) RETURN
!
!
!       Total Downward Flux:
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)//'.'
     &    //phys_suffix(IP_flux_total_down)
        var_name(n_var)='vflx'
        var_type(n_var)='float'
        var_unit(n_var)='wm-2'
        var_long(n_var)='total downward flux'
        n_data(n_var)=n_profile*(n_layer+1)*n_channel
        DO k=1, n_channel
          DO i=1, n_layer+1
            DO l=1, n_profile
              point=l+(i-1+(k-1)*(n_layer+1))*n_profile
              data_fl(point, n_var)=flux_down(l, i-1, k)
            ENDDO
          ENDDO
        ENDDO
        CALL write_cdl(ierr
     &    , file_name(1: length_name+1+len_file_suffix)
     &    , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &    , n_dimension, dimension_name, dimension_type
     &    , dimension_unit
     &    , dimension_long, dimension_size
     &    , dimension_array_int, dimension_array_fl
     &    , n_var, var_name, var_type, var_unit, var_long
     &    , n_dimension_var, list_dimension_var
     &    , n_data, data_int, data_fl
     &    )
        IF (ierr /= i_normal) RETURN
!
      ENDIF
!
!
!     Heating Rates:
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_heating_rate)
!     Reset the third dimension.
      dimension_size(3)=n_layer
      IF (trim(name_vert_coord) == 'plev') THEN
        DO i=1, n_layer
          dimension_array_fl(i, 3)=p(1, i)
        ENDDO
      ELSE IF (trim(name_vert_coord) == 'level') THEN
        n_data(1)=n_profile*n_layer
        DO i=1, n_layer
          dimension_array_int(i, 3)=i
        ENDDO
      ENDIF
      
      var_name(n_var)='hrts'
      var_type(n_var)='float'
      var_unit(n_var)='k.d-1'
      var_long(n_var)='heating rates'
      n_data(n_var)=n_profile*n_layer*n_channel
      DO k=1, n_channel
        DO i=1, n_layer
          DO l=1, n_profile
            point=l+(i-1+(k-1)*n_layer)*n_profile
            data_fl(point, n_var)=heating_rate(l, i, k)
          ENDDO
        ENDDO
      ENDDO
      CALL write_cdl(ierr
     &  , file_name(1: length_name+1+len_file_suffix)
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
        IF (ierr /= i_normal) RETURN

      IF (control%l_contrib_func) THEN
!     Contribution function (intensity):
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)//'.'//phys_suffix(IP_contrib_funci)

        var_name(n_var)='cfi'
        var_type(n_var)='float'
        var_unit(n_var)=''
        var_long(n_var)='contribution function (intensity)'
        n_data(n_var)=n_profile*(n_layer)*n_channel
        DO k=1, n_channel
          DO i=1, n_layer
            DO l=1, n_profile
              point=l+(i-1+(k-1)*(n_layer))*n_profile
              data_fl(point, n_var)=contrib_funci(l, i, k)
            ENDDO
          ENDDO
        ENDDO
        CALL write_cdl(ierr
     &    , file_name(1: length_name+1+len_file_suffix)
     &    , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &    , n_dimension, dimension_name, dimension_type
     &    , dimension_unit
     &    , dimension_long, dimension_size
     &    , dimension_array_int, dimension_array_fl
     &    , n_var, var_name, var_type, var_unit, var_long
     &    , n_dimension_var, list_dimension_var
     &    , n_data, data_int, data_fl
     &    )
          IF (ierr /= i_normal) RETURN

!     Contribution function (flux):
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)//'.'//phys_suffix(IP_contrib_funcf)

        var_name(n_var)='cff'
        var_type(n_var)='float'
        var_unit(n_var)=''
        var_long(n_var)='contribution function (flux)'
        n_data(n_var)=n_profile*(n_layer)*n_channel
        DO k=1, n_channel
          DO i=1, n_layer
            DO l=1, n_profile
              point=l+(i-1+(k-1)*(n_layer))*n_profile
              data_fl(point, n_var)=contrib_funcf(l, i, k)
            ENDDO
          ENDDO
        ENDDO
        CALL write_cdl(ierr
     &    , file_name(1: length_name+1+len_file_suffix)
     &    , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &    , n_dimension, dimension_name, dimension_type
     &    , dimension_unit
     &    , dimension_long, dimension_size
     &    , dimension_array_int, dimension_array_fl
     &    , n_var, var_name, var_type, var_unit, var_long
     &    , n_dimension_var, list_dimension_var
     &    , n_data, data_int, data_fl
     &    )
        IF (ierr /= i_normal) RETURN
      END IF

!
!
!
      RETURN
      END

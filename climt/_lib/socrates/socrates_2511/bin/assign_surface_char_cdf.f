! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign the surface chracteristics from a netCDF-file.
!
! Purpose:
!   This subroutine reads a single netCDF-file containing the 
!   weighting coefficients for the surface BRDFs or albedos.
!
!   Alternatively, enable the calculation of sea surface albedo if
!   a value of -1 is given for albedo in the .surf input file
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_surface_char_cdf(ierr
     &  , i_angular_integration, isolir, n_band
     &  , file_name, field_name
     &  , zen_0, wave_length_short, wave_length_long
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, n_surf_basis_fnc, ls_brdf_trunc, rho_alb, f_brdf
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_band, nd_brdf_basis_fnc, nd_brdf_trunc
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE rad_pcf
!
!
      IMPLICIT NONE
!
!
!     INCLUDE HEADER FILES.
!
!
!
!     Declaration of variables:
!
      INTEGER   !,Intent(OUT)
     &    ierr
!            Error flag
      INCLUDE 'cdl_struc.finc'
!
!     Sizes of arrays
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Allowed size for profiles
     &  , nd_latitude
!           Allowed size for latitudes
     &  , nd_longitude
!           Allowed size for longitudes
     &  , nd_band
!           Allowed size for spectral bands
     &  , nd_brdf_basis_fnc
!           Allowed size for BRDF basis functions
     &  , nd_brdf_trunc
!           Allowed order of truncation for BRDF basis functions
!
      INTEGER, Intent(IN) ::
     &    i_angular_integration
!           Type of angular integration
     &  , isolir
!           Spectral region
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
     &  , field_name*(*)
!           Name of input field
!
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of spectral bands
!
      REAL  (RealK), Intent(IN) ::
     &    zen_0(nd_profile)
!           Secants or cosines of solar zenith angles
     &  , wave_length_short(nd_band)
!           Shorter wavelength limits
     &  , wave_length_long(nd_band)
!           Longer wavelength limits

      INTEGER, Intent(INOUT) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
!
!
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
      INTEGER, Intent(OUT) ::
     &    n_surf_basis_fnc
!           Number of surface basis functions
     &  , ls_brdf_trunc
!           Order of truncation applied to BRDFs
      REAL  (RealK), Intent(OUT) ::
     &    rho_alb(nd_profile, nd_brdf_basis_fnc, nd_band)
!           Surface characteristics
     &  , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!           Array of BRDF basis terms
!
!
!     Local Variables
!
      INTEGER
     &    l
!           Loop variable
     &  , id_lat
!           CDL `index' for latitude
     &  , id_long
!           CDL `index' for longitude
     &  , id_band
!           CDL `index' for spectral band
     &  , id_basis
!           CDL `index' for basis function
     &  , i_lat
!           Loop variable
     &  , i_long
!           Loop variable
     &  , i_band
!           Loop variable
     &  , i_basis
!           Loop variable
     &  , stride_lat
!           Stride through latitudes
     &  , stride_long
!           Stride through longitudes
     &  , stride_band
!           Stride through spectral bands
     &  , stride_basis
!           Stride through basis functions
     &  , i_brdf_basis
!           Representation of the BRDF basis functions
     &  , rho_alb1
!           Used to check whether sea albedo should be calculated
      LOGICAL
     &    l_found_band
!           Flag for longitude
     &  , l_found_basis
!           Flag for basis functions
!
      REAL (RealK) ::
     &    sea_albedo(n_band)

!
!     Functions invoked:
      INTEGER
     &    calc_cdl_stride
!           Function to calculate the stride in a CDL-array
      EXTERNAL
     &    calc_cdl_stride
!
!     Subroutines called:
      EXTERNAL
     &     read_cdf, assign_horiz_cdl
!
!
!
!     Read in the file of surface characteristics.
      CALL read_cdf(ierr
     &  , file_name
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit, dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Check the contents of the file.
!
!     Only one variable is permitted here.
!
      IF (n_var > 1) THEN
        WRITE(iu_err, '(2(/a))')
     &    '*** Error: There are too many variables in the file '
     &    , file_name
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Check and assign the horizontal dimensions.
      CALL assign_horiz_cdl(ierr, file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_fl
     &  , id_lat, n_latitude, latitude, id_long, n_longitude, longitude
     &  , n_profile
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Spectral band:
!     This may not be present. Grey surfaces do not require the
!     bands to be given explicitly.
!
      l_found_band=.false.
      id_band=1
      DO WHILE ( (id_band <= n_dimension).AND.(.NOT.l_found_band) )
        IF (dimension_name(id_band)(1:4) == 'band') THEN
!
          l_found_band=.true.
          IF (dimension_type(id_band)(1:3) /= 'int') THEN
            WRITE(iu_err, '(/a)')
     &        '*** Error: Illegal type for spectral bands '
     &        //'while reading '//field_name//'.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!         Check the number of spectral bands.
          IF (n_band /= dimension_size(id_band)) THEN
            WRITE(iu_err, '(/a, /a)')
     &        '*** Error: Number of spectral bands is not '
     &        //'consistent with the spectral file after '
     &        , 'reading '//field_name//'.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!
        ELSE
!
          id_band=id_band+1
!
        ENDIF
      ENDDO
!
!
!
!     Basis functions:
!
      l_found_basis=.false.
      id_basis=1
      DO WHILE ( (id_basis <= n_dimension).AND.(.NOT.l_found_basis) )
        IF (dimension_name(id_basis)(1:5) == 'basis') THEN
!
          l_found_basis=.true.
          IF (dimension_type(id_basis)(1:3) /= 'int') THEN
            WRITE(iu_err, '(/a)')
     &        '*** Error: Illegal type for basis functions '
     &        //'while reading '//field_name//'.'
            ierr=i_err_fatal
            RETURN
          ENDIF
          n_surf_basis_fnc=dimension_size(id_basis)
!
        ELSE
!
          id_basis=id_basis+1
!
        ENDIF
      ENDDO
!
      IF (.NOT.l_found_basis) THEN
        WRITE(iu_err, '(/a)')
     &    '*** Error: No coordinate for basis functions found '
     &    //'while reading '//field_name//'.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!
!
!     The data in the CDL file will be stored with the last
!     `index' changing most rapidly (C-convention), but there
!     will be a common stride between values at the same
!     point.
!
      stride_lat=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_long
     &   , dimension_size, nd_cdl_dimen)
!
      IF (l_found_band) THEN
        stride_band=calc_cdl_stride(n_dimension_var(1)
     &     , list_dimension_var(1, 1), id_band
     &     , dimension_size, nd_cdl_dimen)
      ELSE
        stride_band=0
      ENDIF
!
      stride_basis=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_basis
     &   , dimension_size, nd_cdl_dimen)
!
!     Assign the data.
      DO i_basis=1, n_surf_basis_fnc
!       In the case of a grey field the stride over bands will be 0
!       so the code may be written generally.
        DO i_lat=1, n_latitude
          DO i_long=1, n_longitude
            l=i_long+(i_lat-1)*n_longitude
            rho_alb1 = int(data_fl(1+stride_lat*(i_lat-1)
     &        +stride_long*(i_long-1)+stride_basis*(i_basis-1), 1))

            IF (rho_alb1 .eq. -1) THEN
!             Calculate the sea surface albedo
              CALL seaalbedo_driver(zen_0(l), n_band,
     &          wave_length_short, wave_length_long, sea_albedo)

              DO i_band=1, n_band
                rho_alb(l,i_basis,i_band) = sea_albedo(i_band)
              ENDDO
            ELSE
              DO i_band=1, n_band
                rho_alb(l, i_basis, i_band)
     &            =data_fl(1+stride_lat*(i_lat-1)+stride_long*(i_long-1)
     &            +stride_band*(i_band-1)+stride_basis*(i_basis-1), 1)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!
!
!     Set the basis functions:
!
      IF (i_angular_integration == IP_two_stream) THEN
!
!       Explicit basis functions are not required, but if there is
!       no separate direct albedo the diffuse value must be copied
!       into the direct field.
!
        IF (isolir == IP_solar) THEN
!
          IF (n_surf_basis_fnc == 1) THEN
            DO i_band=1, n_band
              DO l=1, n_profile
                rho_alb(l, IP_surf_alb_dir, i_band)
     &            =rho_alb(l, IP_surf_alb_diff, i_band)
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF
!
      ELSE IF (i_angular_integration == IP_ir_gauss) THEN
!
!       No direct component is required here.
        continue
!
      ELSE IF (i_angular_integration == IP_spherical_harmonic) THEN
!
        WRITE(iu_stdout, '(/a)') 
     &    'Enter the type of brdf basis function'
        READ(iu_stdin, *) i_brdf_basis
!
        IF (i_brdf_basis == IP_surface_lambertian) THEN
!
!         The Lambertian surface has a single spherically symmetric
!         term: the BRDF is truncated at the zeroth order.
          IF (n_surf_basis_fnc /= 1) THEN
            WRITE(iu_err, '(/a, /a, /a)') 
     &        '*** Error: The file of surface CHARACTERistics '
     &        , 'contains too many weights and is not consistent'
     &        , 'with a lambertian surface.'
            ierr=i_err_fatal
            RETURN
          ENDIF
          ls_brdf_trunc=0
!         By defining F_BRDF as below the weighting becomes the albedo.
          f_brdf(1, 0, 0, 0)=4.0e+00_RealK
!
        ELSE IF (i_brdf_basis == IP_surface_roujean) THEN
!
          WRITE(iu_err, '(/a)') 
     &      '*** Error: Roujean''s scheme has not yet been implemented.'
          ierr=i_err_fatal
          RETURN
!
!        ELSE IF (I_BRDF_BASIS.EQ.IP_SURFACE_LOMMEL_SEELIGER) THEN
!!
!!         Crude provisional 4-stream Lommel-Seeliger function.
!          F_BRDF(1, 0, 0, 0)=13.48947919306655
!          F_BRDF(1, 1, 0, 0)=-3.0152211617935
!          F_BRDF(1, 0, 1, 0)=-3.0152211617935
!          F_BRDF(1, 1, 1, 0)=1.10089123515584
!
        ELSE
!
          WRITE(iu_err, '(/a)')
     &      '*** Error: An illegal surface scheme has been specified.'
          ierr=i_err_fatal
          RETURN
!
        ENDIF
!
      ENDIF
            
!
!
!
      RETURN
      END

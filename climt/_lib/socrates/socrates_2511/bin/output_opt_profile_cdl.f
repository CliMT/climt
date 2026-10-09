! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to write optical properties to a CDL-file.
!
! Method:
!   According to the pointer to the list of input files ordrered
!   by pressure spectrally averaged scattering propeties are read
!   from the input and written to a CDL-file.
!
!- ---------------------------------------------------------------------
      SUBROUTINE output_opt_profile_cdl(ierr 
     &  , n_file, p, filename, pointer, latitude, longitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , nd_input_file
     &  , nd_phase_term, nd_band)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays
      INTEGER, Intent(IN) ::
     &    nd_phase_term
!           Size allocated for terms in the phase function
     &  , nd_band
!           Size allocated for spectral bands
     &  , nd_input_file
!           Size allocated for input files
!
!
!
!
!     Dummy arguments:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
!
      CHARACTER  !, Intent(IN)
     &    filename(nd_input_file)*80
!           Names of the files holding scattering data
      INTEGER, Intent(IN) ::
     &    n_file
!           Number of input files
     &  , pointer(nd_input_file)
!           Array of pointers to sets of input data sorted by pressure
!
      REAL  (RealK), Intent(IN) ::
     &    p(nd_input_file)
!           Pressures at which scattering data are provided
     &  , latitude
!           Latitude of the point
     &  , longitude
!           Longitude of the point
!
!
!     Local variables:
!
!     CDL structures:
      INCLUDE 'cdl_struc.finc'
      CHARACTER
     &    file_cdl*80
!           Output file of optical properties
!
!     Optical properties
      INTEGER
     &    n_phf_term
!           Number of terms in the phase function
     &  , n_phf_term_in
!           Number of terms in the phase function 
!           in the current input file
     &  , n_band
!           Number of spectral bands
     &  , n_band_in
!           Number of spectral bands in the current input file
      REAL  (RealK) ::
     &    absorption(nd_band)
!           Absorption in the current input file
     &  , scattering(nd_band)
!           Scattering in the current input file
     &  , phase_fnc(nd_phase_term, nd_band)
!           Phase function in the current input file
!
      INTEGER
     &    i
!           Loop variable
     &  , j
!           Loop variable
     &  , k
!           Loop variable
!
      EXTERNAL
     &    read_averaged_scatter, write_cdl
!
!
!
!     Read in the first file to get the size of the data.
      CALL read_averaged_scatter(ierr
     &   , filename(1)
     &   , n_band, absorption
     &   , scattering, n_phf_term, phase_fnc
     &   , nd_band, nd_phase_term)
      IF (ierr /= i_normal) STOP
!
!     Set the characteristics of the field of data.
      n_dimension=5
      dimension_name(1)='lat'
      dimension_name(2)='lon'
      dimension_name(3)='plev'
      dimension_name(4)='mom'
      dimension_name(5)='band'
      dimension_type(1)='float'
      dimension_type(2)='float'
      dimension_type(3)='float'
      dimension_type(4)='int'
      dimension_type(5)='int'
      dimension_unit(1)='degree'
      dimension_unit(2)='degree'
      dimension_unit(3)='pa'
      dimension_unit(4)='none'
      dimension_unit(5)='none'
      dimension_long(1)='latitude'
      dimension_long(2)='longitude'
      dimension_long(3)='pressure level'
      dimension_long(4)='moment'
      dimension_long(5)='spectral band'
      dimension_size(1)=1
      dimension_size(2)=1
      dimension_size(3)=n_file
      dimension_size(4)=n_phf_term
      dimension_size(5)=n_band
!
      dimension_array_fl(1, 1)=latitude
      dimension_array_fl(1, 2)=longitude
      DO i=1, n_file
        dimension_array_fl(i, 3)=p(pointer(i))
      ENDDO
      DO i=1, n_phf_term
        dimension_array_int(i, 4)=i
      ENDDO
      DO i=1, n_band
        dimension_array_int(i, 5)=i
      ENDDO
!
      n_var=3
!
      var_name(1)='abs'
      var_type(1)='float'
      var_long(1)='absorption'
      var_unit(1)='m-1'
      n_dimension_var(1)=4
      list_dimension_var(1, 1)=1
      list_dimension_var(2, 1)=2
      list_dimension_var(3, 1)=3
      list_dimension_var(4, 1)=5
      n_data(1)=n_file*n_band
!
      var_name(2)='scat'
      var_type(2)='float'
      var_long(2)='scattering'
      var_unit(2)='m-1'
      n_dimension_var(2)=4
      list_dimension_var(1, 2)=1
      list_dimension_var(2, 2)=2
      list_dimension_var(3, 2)=3
      list_dimension_var(4, 2)=5
      n_data(2)=n_file*n_band
!
      var_name(3)='phf'
      var_type(3)='float'
      var_long(3)='phase function'
      var_unit(3)='none'
      n_dimension_var(3)=5
      list_dimension_var(1, 3)=1
      list_dimension_var(2, 3)=2
      list_dimension_var(3, 3)=3
      list_dimension_var(4, 3)=4
      list_dimension_var(5, 3)=5
      n_data(3)=n_file*n_band*n_phf_term
!
      j=1
      DO WHILE (j <= n_file)
!       Read the input files in the sorted order.
        CALL read_averaged_scatter(ierr, filename(pointer(j))
     &     , n_band_in, absorption
     &     , scattering, n_phf_term_in, phase_fnc
     &     , nd_band, nd_phase_term)
        IF (ierr /= i_normal) STOP
!
!       Check for consistency with the previous sizes.
        IF ( (n_band_in /= n_band).OR.
     &       (n_phf_term /= n_phf_term_in) ) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: The files of scattering properties are '
     &      //'mutually inconsistent.'
          ierr=i_err_fatal
          RETURN
        ENDIF
!
!       Convert the data to the CDL-format
        DO i=1, n_band
          data_fl(j+(i-1)*n_file, 1)
     &      =absorption(i)
          data_fl(j+(i-1)*n_file, 2)
     &      =scattering(i)
          DO k=1, n_phf_term
            data_fl(j+(k-1+(i-1)*n_phf_term)*n_file, 3)
     &        =phase_fnc(k, i)
          ENDDO
        ENDDO
!
        j=j+1
!
      ENDDO
      
!
!     Write the data out as a CDL-file.
      WRITE(iu_stdout, '(/a)')
     &   'Give the name of the file to be used for output.'
      READ(iu_stdin, '(a)') file_cdl
      CALL write_cdl(ierr
     &   , file_cdl
     &   , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &   , n_dimension, dimension_name, dimension_type, dimension_unit
     &   , dimension_long, dimension_size
     &   , dimension_array_int, dimension_array_fl
     &   , n_var, var_name, var_type, var_unit, var_long
     &   , n_dimension_var, list_dimension_var
     &   , n_data, data_int, data_fl
     &   )
!
!
!
      STOP
      END

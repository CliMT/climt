! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert densities to mass fractions.
!
! Method:
!   Fields are read in as CDL files.
!   The density of air is calculated from T and q whence
!   the mixing ratio is calculated.
!
!- ---------------------------------------------------------------------
      PROGRAM dentomix
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE error_pcf
      USE gas_list_pcf
      USE rad_ccf, ONLY: r_gas_dry, mol_weight_air
!
!
      IMPLICIT NONE
!
!
!    
!
!     Declaration of variables.
      CHARACTER
     &    file_density*80
!           Name of file with denisty
     &  , file_t*80
!           Name of temperature file
     &  , file_q*80
!           Name of humidity file
     &  , file_mix*80
!           Name of file with mixing ratio
     &  , name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      INTEGER
     &    ierr
!           Error flag
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
     &  , t(npd_profile, npd_layer+1)
!           Temperatures
     &  , q(npd_profile, npd_layer+1)
!           Humidity mixing ratios
     &  , density(npd_profile, npd_layer+1)
!           Mass density
     &  , mix_ratio(npd_profile, npd_layer+1)
!           Mixing ratio
     &  , density_air
!           Denisty of air
     &  , ratio_molar_weight
!           Molecular weight of dry air/ molecular weight of water
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, output_vert_cdl
!
!
      data l_vert_coord/.false./
      data ierr/i_normal/
      data n_latitude/0/
      data n_longitude/0/
!
      ratio_molar_weight
     &  = mol_weight_air/(molar_weight(ip_h2o)*1.0E-03_RealK)
!
      WRITE(iu_stdout, '(/a)')
     &   'enter the name of the file containing the mass densities.'
      READ(iu_stdin, '(a)') file_density
      CALL assign_input_vert_cdl(ierr
     &  , file_density, 'densities', l_vert_coord, name_vert_coord
     &  , .true., n_level, .true.
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, density
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
      WRITE(iu_stdout, '(/a)')
     &   'enter the name of the file containing the temperatures.'
      READ(iu_stdin, '(a)') file_t
      CALL assign_input_vert_cdl(ierr
     &  , file_t, 'temperatures', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, t
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the file containing the specific humidity.'
      READ(iu_stdin, '(a)') file_q
      CALL assign_input_vert_cdl(ierr
     &  , file_q, 'spcific humidities', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, q
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
      DO i=1, n_level
        DO l=1, n_profile
          density_air=p(l, i)/(r_gas_dry*t(l, i)
     &      *(1.+(ratio_molar_weight-1.)*q(l, i)))
          mix_ratio(l, i)=density(l, i)/density_air
        ENDDO
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_mix
      CALL output_vert_cdl(ierr
     &  , file_mix
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level, trim(name_vert_coord), len(trim(name_vert_coord)), p
     &  , 'mixrat', 'float', 'none', 'mixing ratio', mix_ratio
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      STOP
      END

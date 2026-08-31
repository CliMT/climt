! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to convert from dewpoints to q.
!
! Method:
!   The dewpoints are read in and an approprite routine
!   is called to determine the saturated mixing ratio.
!
!- ---------------------------------------------------------------------
      PROGRAM dewtoq
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE error_pcf
      USE qsat_alg_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
!
      CHARACTER
     &    file_dew*80
!           Name of file with dewpoint
     &  , file_mix*80
!           Name of file with mixing ratio
     &  , name_vert_coord*24
!           Name of vertical coordinate
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
     &  , n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
     &  , i_qsat_alg
!           Algorithm for calculating saturated mixing ratio
     &  , i
!           Loop variable
     &  , l
!           Loop variable
!
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , p(npd_profile, npd_layer+1)
!           Pressure levels
     &  , dewpoint(npd_profile, npd_layer+1)
!           Dewpoint
     &  , mix_ratio(npd_profile, npd_layer+1)
!           Mixing ratio
     &  , t_transition
!           Transition temperature for saturation w.r.t. water/ice
!           In the goff-gratsch formula.
!
      LOGICAL
     &    l_vert_coord
!           Flag asserting that the vertical coordinate has been set
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, output_vert_cdl
     &  , qsat, qsat_gg, qsat_gg_ice
!
!
      data
     &    l_vert_coord/.false./
      data ierr/i_normal/
!
!
!
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the file containing the dew points.'
      READ(iu_stdin, '(a)') file_dew
      CALL assign_input_vert_cdl(ierr
     &  , file_dew, 'dew points', l_vert_coord, name_vert_coord
     &  , .true., n_level, .true.
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, dewpoint
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!     Determine algorithm for conversion of dewpoints.
      WRITE(iu_stdout, '(/a)') 'enter the algorithm for conversion.'
1     read(iu_stdin, *, iostat=ios) i_qsat_alg
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') '+++ unrecognized response.'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(/a)') 'please re-enter.'
          goto 1
        ENDIF
      ENDIF
      IF ( (i_qsat_alg /= IP_um_look_up).AND.
     &     (i_qsat_alg /= IP_goff_gratsch) ) THEN
        WRITE(iu_err, '(/a)') '+++ invalid algorithm.'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(/a)') 'please re-enter.'
          goto 1
        ENDIF
      ENDIF
!
      IF (i_qsat_alg == IP_um_look_up) THEN
        DO i=1, n_level
          DO l=1, n_profile
            CALL qsat(mix_ratio(l, i), dewpoint(l, i), p(l, i), 1)
          ENDDO
        ENDDO
      ELSE IF (i_qsat_alg == IP_goff_gratsch) THEN
        WRITE(iu_stdout, '(a)') 'enter the transition temperature for '
     &    //'saturation with respect to water or ice.'
        READ(iu_stdin , *) t_transition
        DO i=1, n_level
          DO l=1, n_profile
            IF (dewpoint(l, i) >= t_transition) THEN
              CALL qsat_gg(mix_ratio(l, i)
     &          , dewpoint(l, i), p(l, i))
            ELSE
              CALL qsat_gg_ice(mix_ratio(l, i)
     &          , dewpoint(l, i), p(l, i))
            ENDIF
          ENDDO
        ENDDO
      ENDIF
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

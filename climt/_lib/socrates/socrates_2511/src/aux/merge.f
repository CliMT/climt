! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to merge a CDL-field on to a background.
!
! Method:
!   Two CDL-files are read, one giving the background field
!   and the other the field to be merged. Options for merging
!   are obtained and the merged field is written out.
!
!- ---------------------------------------------------------------------
      PROGRAM merge
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE error_pcf
      USE method_merge_pcf
      USE interpolate_p_mod
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
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
!
!     Common arrays:
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
!
!     Foreground field:
!
      CHARACTER
     &    file_front*80
!           File to be merged
      INTEGER
     &    n_latitude_front
!           Number of latitudes
     &  , n_longitude_front
!           Number of longitudes
     &  , n_profile_front
!           Number of profiles in merging field
     &  , n_level_front
!           Number of levels in merging field
      REAL  (RealK) ::
     &    p_front(npd_layer+1)
!           Pressure in merging field
     &  , field_front(npd_profile, npd_layer+1)
!           Field to be merged
     &  , x_front(npd_layer+1)
!           Abscissa for merged field
     &  , y_front(npd_layer+1)
!           Ordinate for merged field
     &  , y2_front(npd_layer+1)
!           2nd-derivative for merged field
     &  , a_front(npd_layer+1)
!           Data in merging field
!
!     Background field:
!
      CHARACTER
     &    file_back*80
!           File containing background field
      INTEGER
     &    n_latitude_back
!           Number of latitudes
     &  , n_longitude_back
!           Number of longitudes
     &  , n_profile_back
!           Number of profiles in the background
     &  , n_level_back
!           Number of levels in background
      REAL  (RealK) ::
     &    p_back(npd_layer+1)
!           Pressure in background field
     &  , field_back(npd_profile, npd_layer+1)
!           Background field
     &  , x_back(npd_layer+1)
!           Abscissa for background field
     &  , y_back(npd_layer+1)
!           Ordinate for background field
     &  , y2_back(npd_layer+1)
!             2nd-derivative for background field
     &  , a_back(npd_layer+1)
!           Data in background field
!
!     Composite Field:
!
      CHARACTER
     &    file_composite*80
!           File to be created
     &  , name_composite*24
!           Name of the composite field
     &  , unit_composite*24
!           Unit of the composite field
     &  , long_composite*24
!           Long name of the composite field
      INTEGER
     &    n_latitude_composite
!           Number of latitudes
     &  , n_longitude_composite
!           Number of longitudes
     &  , n_profile_composite
!           Number of profiles in the composite field
     &  , n_level_composite
!           Number of levels in composite
      REAL  (RealK) ::
     &    p_composite(npd_layer+1)
!           Pressure in composite field
     &  , field_composite(npd_profile, npd_layer+1)
!           Composite field
!
!
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
     &  , l_splined_front
!           Merging spline flag
     &  , l_splined_back
!           Background spline flag
     &  , l_top_only
!           Logical to merge only top
!
      INTEGER
     &    i_method
!           Method of merging
     &  , i_mode
!           Fitting mode
     &  , i
!           Loop variable
     &  , l
!           Loop variable
!
      CHARACTER
     &    char_yn*1
!           Character response variable
     &  , name_vert_coord*24
!           Name of vertical coordinate
!
      REAL  (RealK) ::
     &    p_front_low
!           Low pressure in merged field
     &  , p_front_high
!           High pressure in merged field
     &  , p_merge_low
!           Low merging pressure
     &  , p_merge_high
!           High merging pressure
     &  , field_front_int
!           Interpolated merging value
     &  , field_back_int
!           Interpolated background value
     &  , weight_front
!           Weight for merging field
     &  , weight_back
!           Weight for background field
     &  , pp
!           Temporary pressure variable
!
!     Numerical precision:
      REAL  (RealK) ::
     &    tol_p
!           Tolerance for equality of pressures
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, merge_pressure
     &  , output_vert_cdl
!
      data
     &    l_vert_coord/.false./
      data ierr/i_normal/
!
!
!
!
!     Set the numerical precision
      tol_p=1.0e+03_RealK*epsilon(tol_p)
!
      WRITE(iu_stdout, '(/a)')
     &  'Enter the name of the background field.'
      READ(iu_stdin, '(a)') file_back
      CALL assign_input_vert_cdl(ierr
     &  , file_back, 'background', l_vert_coord, name_vert_coord
     &  , .true., n_level_back, .NOT.l_vert_coord
     &  , n_latitude_back, latitude, n_longitude_back, longitude
     &  , 1
     &  , n_profile_back, n_level_back
     &  , p_back, field_back
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!     Read the field to be merged in.
      WRITE(iu_stdout, '(/a)') 'Enter the name of the file containing '
     &  //'the field to be merged.'
      READ(iu_stdin, '(a)') file_front
      CALL assign_input_vert_cdl(ierr
     &  , file_front, 'foreground', l_vert_coord, name_vert_coord
     &  , .false., n_level_front, .true.
     &  , n_latitude_front, latitude, n_longitude_front, longitude
     &  , 1
     &  , n_profile_front, n_level_front
     &  , p_front, field_front
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!     Carr out some checks.
      IF ( (n_profile_front /= n_profile_back).OR.
     &     (n_longitude_front /= n_longitude_back).OR.
     &     (n_latitude_front /= n_latitude_back) ) THEN
        WRITE(iu_err, '(/a)') 
     &    '*** Error: The horizontal structure of the field to be'
     &    , 'merged does not match that of the background.'
        ierr=i_err_fatal
        STOP
      ENDIF
!
      p_front_low=p_front(1)
      p_front_high=p_front(n_level_front)
      WRITE(iu_stdout, '(a)')
     &  'This profile contains data between pressure levels: '
      WRITE(iu_stdout, '(1pe12.5, a9, 1pe12.5, a4)')
     &  p_front_low, ' pa, and ', p_front_high
     &  , ' pa.'
!
!     Determine whether only the top should be merged.
      WRITE(iu_stdout, '(/a)')
     &  'Do you wish to merge only at the top? (y/n)'
1     read(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
        l_top_only=.false.
      ELSE IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
        l_top_only=.true.
      ELSE
        WRITE(iu_err, '(a)') '+++ Invalid response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
          goto 1
        ENDIF
      ENDIF
!
!     Determine the method of merging.
      WRITE(iu_stdout, '(/a)') 'Enter the method of merging data.'
2     read(iu_stdin, *, iostat=ios) i_method
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Illegal input:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'Please re-enter'
          goto 2
        ENDIF
      ENDIF
      IF ( (i_method < 1).OR.(i_method > p_method_merge) ) THEN
        WRITE(iu_err, '(a)') '+++ Input out of range:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
          goto 2
        ENDIF
      ENDIF
!
      IF (i_method == IP_merge_linear) THEN
        WRITE(iu_stdout, '(a)') 'Enter merging pressures.'
3       read(iu_stdin, *, iostat=ios) p_merge_low, p_merge_high
        IF (ios /= 0) THEN
          WRITE(iu_err, '(a)') '+++ Unrecognized response:'
          IF (lock_code(.true.)) THEN
            INCLUDE 'batch_error_main.finc'
          ELSE
            WRITE(iu_stdout, '(a)') 'Please re-type.'
            goto 3
          ENDIF
        ENDIF
        IF (p_merge_high < p_merge_low) THEN
          WRITE(iu_err, '(a)') '+++ Unsuitable merging pressures:'
          IF (lock_code(.true.)) THEN
            INCLUDE 'batch_error_main.finc'
          ELSE
            WRITE(iu_stdout, '(a)') 'Please re-enter.'
            goto 3
          ENDIF
        ENDIF
        IF (p_merge_low < p_front_low) THEN
          IF (p_front_low > p_back(1)) THEN
            WRITE(iu_err, '(a)') '+++ Unsuitable merging pressures:'
            IF (lock_code(.true.)) THEN
              INCLUDE 'batch_error_main.finc'
            ELSE
              WRITE(iu_stdout, '(a)') 'Please re-enter.'
              goto 3
            ENDIF
          ELSE
            p_merge_low=p_front_low
          ENDIF
        ENDIF
        IF (p_merge_high > p_front_high) THEN
          IF (p_front_high < p_back(n_level_back)) THEN
            WRITE(iu_err, '(a)') 'Unsuitable merging pressures:'
            IF (lock_code(.true.)) THEN
              INCLUDE 'batch_error_main.finc'
            ELSE
              WRITE(iu_stdout, '(a)') 'Please re-enter.'
              goto 3
            ENDIF
          ELSE
            p_merge_high=p_front_high
          ENDIF
        ENDIF
!
!       In this case an interpolation rule is required.
        WRITE(iu_stdout, '(/a)') 'Enter mode of interpolation.'
4       read(iu_stdin, *, iostat=ios) i_mode
        IF (ios /= 0) THEN
          WRITE(iu_err, '(a)') '+++ Unrecognized response:'
          IF (lock_code(.true.)) THEN
            INCLUDE 'batch_error_main.finc'
          ELSE
            WRITE(iu_stdout, '(a)') 'Please re-enter.'
            goto 4
          ENDIF
        ENDIF
        IF ( (i_mode /= IP_1_lin_lin).AND.
     &       (i_mode /= IP_1_log_lin).AND.
     &       (i_mode /= IP_1_lin_log).AND.
     &       (i_mode /= IP_1_log_log).AND.
     &       (i_mode /= IP_3_lin_lin).AND.
     &       (i_mode /= IP_3_log_lin).AND.
     &       (i_mode /= IP_3_lin_log).AND.
     &       (i_mode /= IP_3_log_log) ) THEN
          WRITE(iu_err, '(a)') '+++ Unsupported mode:'
          IF (lock_code(.true.)) THEN
            INCLUDE 'batch_error_main.finc'
          ELSE
            WRITE(iu_stdout, '(a)') 'Please re-enter.'
            goto 4
          ENDIF
        ENDIF
      ELSE
!       Set a default mode in these cases.
        i_mode=IP_1_lin_lin
      ENDIF
!
!     Construct the merged pressure field.
      CALL merge_pressure(ierr
     &  , n_level_front, p_front, n_level_back, p_back
     &  , n_level_composite, p_composite
     &  , p_front_low, p_front_high
     &  , l_top_only
     &  )
      IF (ierr /= i_normal) STOP
!
!     Work out the interpolated fields in the merged grid
      DO l=1, n_profile_front
        DO i=1, n_level_front
          a_front(i)=field_front(l, i)
        ENDDO
        DO i=1, n_level_back
          a_back(i)=field_back(l, i)
        ENDDO
      ENDDO
      DO i=1, n_level_composite
        pp=p_composite(i)
        IF (i_method == IP_merge_direct) THEN
          IF ( (pp < p_front_low).OR.(pp > p_front_high) ) THEN
            weight_front=0.0_RealK
          ELSE
            weight_front=1.0_RealK
          ENDIF
          weight_back=1.0_RealK-weight_front
        ELSE IF (i_method == IP_merge_linear) THEN
          IF ( (pp < p_front_low).OR.(pp > p_front_high) ) THEN
            weight_front=0.0_RealK
          ELSE IF ( (pp >= p_front_low).AND.
     &              (pp < p_merge_low) ) THEN
            IF (abs(p_merge_low-p_front_low) > tol_p) THEN
              weight_front=(pp-p_front_low)
     &          /(p_merge_low-p_front_low)
            ELSE
!             This situation arises only in the exceptional case
!             of non_overlapping data.
              weight_front=1.0_RealK
            ENDIF
          ELSE IF ( (pp >= p_merge_low).AND.
     &              (pp < p_merge_high) ) THEN
            weight_front=1.0_RealK
          ELSE IF ( (pp > p_merge_high).AND.
     &              (pp <= p_front_high) ) THEN
            IF (abs(p_merge_high-p_front_high) > tol_p) THEN
              weight_front=(p_front_high-pp)
     &          /(p_front_high-p_merge_high)
            ELSE
!             This situation arises only in the exceptional case
!             of non_overlapping data.
              weight_front=1.0_RealK
            ENDIF
          ENDIF
          weight_back=1.0_RealK-weight_front
        ELSE IF (i_method == IP_merge_zero) THEN
          IF ( (pp >= p_front_low).AND.(pp <= p_front_high) ) THEN
            weight_front=1.0_RealK
          ELSE
            weight_front=0.0_RealK
          ENDIF
          weight_back=0.0_RealK
        ENDIF
!
        DO l=1, n_profile_front
!         Interpolate the front and back fields.
          IF ( (pp < p_front_low).OR.(pp > p_front_high) ) THEN
            field_front_int=0.0_RealK
          ELSE
            CALL interpolate_p(n_level_front, p_front
     &        , a_front, x_front, y_front, y2_front
     &        , pp, field_front_int
     &        , i_mode, l_splined_front
     &        )
          ENDIF
          CALL interpolate_p(n_level_back, p_back
     &      , a_back, x_back, y_back, y2_back
     &      , pp, field_back_int
     &      , i_mode, l_splined_back
     &      )
            field_composite(l, i)
     &        =weight_front*field_front_int+weight_back*field_back_int
        ENDDO
      ENDDO
!
!
      WRITE(iu_stdout, '(/a)')
     &  'Enter the name of the output file.'
      READ(iu_stdin, '(a)') file_composite
      WRITE(iu_stdout, '(/a)')
     &  'Enter the name of the output field.'
      READ(iu_stdin, '(a)') name_composite
      WRITE(iu_stdout, '(/a)')
     &  'Enter the unit of the output field.'
      READ(iu_stdin, '(a)') unit_composite
      WRITE(iu_stdout, '(/a)')
     &  'Enter the long-name of the output field.'
      READ(iu_stdin, '(a)') long_composite
!
      n_profile_composite=n_profile_front
      n_longitude_composite=n_longitude_front
      n_latitude_composite=n_latitude_front
!
      CALL output_vert_cdl(ierr
     &  , file_composite
     &  , n_latitude_composite, latitude
     &  , n_longitude_composite, longitude
     &  , n_profile_front
     &  , n_level_composite
     &  , trim(name_vert_coord), len(trim(name_vert_coord)), p_composite
     &  , name_composite, 'float', unit_composite, long_composite
     &  , field_composite
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
      END

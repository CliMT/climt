! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to obtain information about size and shape of particles.
!
Subroutine particle_size_90 &
!
(l_interactive, nd_size, nd_mode, SizeDist, ierr)
!
! Description:
!   This subroutine sets the details of the size distribution 
!   assumed in scattering calculations and allocates the 
!   relevant arrays.
!
! Modules used
  USE realtype_rd
  USE error_pcf
  USE def_size_dist
  USE shape_particle_pcf
  USE measure_particle_pcf
  USE distribution_pcf
  USE def_std_io_icf
  USE spline_fit_mod, ONLY: spline_fit
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  LOGICAL, Intent(IN) :: l_interactive
!       Flag for interactive running
  INTEGER, Intent(OUT) :: nd_size
!       Size allocated for points in the size distribution
  INTEGER, Intent(OUT) :: nd_mode
!       Size allocated for number of modes
  TYPE (STR_size_dist), Intent(OUT) :: SizeDist
!   Size distribution
!      
  INTEGER, Intent(InOut) :: ierr
!       Error flag
!
!
! Local variables.
!
  CHARACTER (LEN=80) :: line
!       Line of input file
  INTEGER :: ios
!       I/O error flag
  INTEGER :: iu_input
!       Unit number for input of the size distribution
  INTEGER :: i
!       Loop variable
  INTEGER :: n_knot
!       Number of knot points
  INTEGER :: i_max
!       Index of mode of number
  LOGICAL :: l_data_region
!       True in regions of data
  REAL  (RealK) :: number_max
!       Peak value of number per unit size range
  REAL  (RealK) :: number_tot
!       Total number density
  REAL  (RealK) :: dim_tmp
!       Temporary variable holding the dimension of the particle
  REAL  (RealK) :: tol_zero
!       The tolerance for flagging zeros
  REAL  (RealK) :: kt
!       Parameter in Mitchell's specification of the 
!       size distribution
!
! Subroutines called:
  EXTERNAL open_file_in
!
! Functions called:
  REAL  (RealK) :: gamma_fnc
!       Gamma-function
  EXTERNAL gamma_fnc
!
!- End of Header
!
!
!
! Set the tolerances used in detecting zeros.
  tol_zero=1.6E+01_RealK*EPSILON(tol_zero)
!
  WRITE (iu_stdout, '(/A)') &
    'Specify the size distribution for averaging.'
!
  DO
!
    READ(iu_stdin, *, IOSTAT=ios) SizeDist%i_distribution
!
    IF (ios == 0) THEN
!
      IF ( (SizeDist%i_distribution <= 0) .OR. &
           (SizeDist%i_distribution >= npd_distribution) ) THEN
        WRITE(iu_err, '(A)') '+++ Illegal response: ' &
          //'This number is outside the supported range.'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-type.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
!
        WRITE (iu_stdout, '(/A)') &
        'Specify the measure used for the particles.'
!
        DO
!
          READ(iu_stdin, *, IOSTAT=ios) SizeDist%i_measure
!
          IF (ios == 0) THEN
!
            IF ( (SizeDist%i_measure /= ip_measure_radius)    .AND. &
                 (SizeDist%i_measure /= ip_measure_max_dimen) .AND. &
                 (SizeDist%i_measure /= ip_measure_proj_radius) ) THEN
              WRITE(iu_err, '(A)') '+++ Illegal response: ' &
                 //'This number does not represent a valid quantity.'
              IF (l_interactive) THEN
                WRITE(iu_stdout, '(A)') 'Please re-type.'
              ELSE
                ierr = i_err_fatal
                RETURN
              ENDIF
            ELSE
              EXIT
            ENDIF
!
          ELSE
!
            WRITE(iu_err, '(/A)') '+++ Unrecognized response: '
            IF (l_interactive) THEN
              WRITE(iu_stdout, '(A)') 'Please re-type.'
            ELSE
              ierr = i_err_fatal
              RETURN
            ENDIF
!
          ENDIF
!
        ENDDO
!
        EXIT
!
      ENDIF
!
!
    ELSE
!
      WRITE(iu_err, '(/A)') '+++ Unrecognized response: '
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(A)') 'Please re-type.'
      ELSE
        ierr = i_err_fatal
        RETURN
      ENDIF
!
    ENDIF
!
  ENDDO
!
!
  IF (SizeDist%i_distribution == ip_external) THEN
!
!   Input read from file.
    CALL get_free_unit(ierr, iu_input)
    CALL open_file_in(IERR, iu_input, &
      'Enter name of file of distributions.')
    IF (ierr /= i_normal) RETURN
!
!   Initialize distribution and read it once to determine 
!   the size to allocate.
    l_data_region             = .FALSE.
    SizeDist%i_shape_particle = ip_shape_null
!
!   Determine the number of points in the size distribution to
!   allocate arrays.
    nd_size=0
    DO
      READ(iu_input, '(A)', IOSTAT=ios) line
      IF (ios < 0) EXIT
      IF (l_data_region) THEN
        IF (line(1:4) /= "*END") THEN
          nd_size = nd_size+1
        ENDIF
      ELSE IF (line(1:11) == "*BEGIN_DATA") THEN
        l_data_region = .TRUE.
      ENDIF
    ENDDO
!   Rewind the file for the actual reading
    REWIND(iu_input)
!

    IF (nd_size > 0) THEN
      ALLOCATE(SizeDist%dimen(0:nd_size))
      ALLOCATE(SizeDist%number(0:nd_size))
      ALLOCATE(SizeDist%aspect(0:nd_size))
    ELSE
      WRITE(iu_err, '(/A)') &
        '*** Error: There are no points in the size distribution.'
    ENDIF
    SizeDist%n_size      = 0
    SizeDist%dimen(0)    = 0.0E+00_RealK
    SizeDist%number(0)   = 0.0E+00_RealK
    SizeDist%aspect(0)   = 1.0E+00_RealK
!
!   Find dimension where number is a maximum.
    l_data_region        = .FALSE.
    number_max           = 0.0E+00_RealK
    i_max                = 1
    DO 
      READ(iu_input, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
!
      IF (l_data_region) THEN
!
        BACKSPACE(iu_input)
        DO
          READ(iu_input, '(A)', IOSTAT=ios) line
          IF (line(1:4) == '*END') EXIT
          BACKSPACE(iu_input)
!         Read the dimension into a temporary variable and decide
!         whether it is sufficiently distinct from the previous
!         value not to cause problems with interpolation.
          READ(iu_input, *, IOSTAT=ios) dim_tmp
          IF (ios /= 0) THEN
            WRITE(iu_err, '(/A)') &
              '*** Error: Size distribution is corrupt'
            ierr = i_err_io
            RETURN
          ENDIF
          IF (dim_tmp > SizeDist%dimen(SizeDist%n_size) + tol_zero) &
               THEN
            SizeDist%n_size = SizeDist%n_size + 1
            BACKSPACE(iu_input)
            IF (SizeDist%i_shape_particle == ip_shape_sphere) THEN
              READ(iu_input, *, IOSTAT=ios) &
                SizeDist%dimen(SizeDist%n_size), &
                SizeDist%number(SizeDist%n_size)
              IF (ios /= 0) THEN
                WRITE(iu_err, '(/A)') &
                  '*** Error: Size distribution is corrupt'
                ierr = i_err_io
                RETURN
              ENDIF
            ELSE IF (SizeDist%i_shape_particle == ip_shape_hexcyl) THEN
              READ(iu_input, *, IOSTAT=ios) &
                SizeDist%dimen(SizeDist%n_size), &
                SizeDist%number(SizeDist%n_size), &
                SizeDist%aspect(SizeDist%n_size)
              IF (ios /= 0) THEN
                WRITE(iu_err, '(/A)') &
                  '*** Error: Size distribution is corrupt'
                ierr=i_err_io
                RETURN
              ENDIF
            ELSE IF (SizeDist%i_shape_particle == &
                     ip_shape_db_undefined) THEN
              READ(iu_input, *, IOSTAT=ios) &
                SizeDist%dimen(SizeDist%n_size), &
                SizeDist%number(SizeDist%n_size)
              IF (ios /= 0) THEN
                WRITE(iu_err, '(/A)') &
                  '*** Error: Size distribution is corrupt'
                ierr = i_err_io
                RETURN
              ENDIF
            ELSE IF (SizeDist%i_shape_particle == ip_shape_null) THEN
              WRITE(iu_err, '(A)') &
                '*** Error: The shape of the particle has not been set.'
              ierr=i_err_io
              RETURN
            ENDIF
!
            IF (SizeDist%number(SizeDist%n_size) > number_max) THEN
              i_max      = SizeDist%n_size
              number_max = SizeDist%number(i_max)
            ENDIF
!
          ENDIF
       
        ENDDO
        l_data_region = .FALSE.
!
      ELSE
!
        IF (line(1:11) == '*BEGIN_DATA') THEN
          l_data_region             = .TRUE.
        ELSE IF (line(1:7) == '*SPHERE') THEN
          SizeDist%i_shape_particle = ip_shape_sphere
        ELSE IF (line(1:7) == '*HEXCYL') THEN
          SizeDist%i_shape_particle = ip_shape_hexcyl
        ELSE IF (line(1:10) == '*UNDEFINED') THEN
!         Note: How best to deal with observed size distributions
!         is not clear, as shapes are not specified. An undefined
!         shape is assumed at the moment. Whether this should be 
!         made more explicit is for future thought.
          SizeDist%i_shape_particle = ip_shape_db_undefined
        ENDIF
!
      ENDIF
!
    ENDDO 
!
    CLOSE(iu_input)
!
!   n_knot is the number of spline knot points. n_size is the
!   number of points in the distribution as read in, excluding
!   the zero point.
    n_knot           = SizeDist%n_size+1
    SizeDist%dimen_0 = SizeDist%dimen(i_max)
!

!   Establish a spline fit to the data.
!   IF (SizeDist%n_size > 1) THEN
!     CALL spline_fit(n_knot, SizeDist%dimen, SizeDist%number, &
!       SizeDist%d2_number)
!     CALL spline_fit(n_knot, SizeDist%dimen, SizeDist%aspect, &
!       SizeDist%d2_aspect)
!   ENDIF

!
!
!
  ELSE IF (SizeDist%i_distribution == ip_log_normal) THEN
!
!   Log-normal distribution.
!
!   Spheres are assumed here.
    SizeDist%i_shape_particle = ip_shape_sphere
!
    WRITE(iu_stdout, '(/A)') &
      'Enter the overall number density.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) SizeDist%n0_ln
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (L_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!
!
    WRITE(iu_stdout, '(/A)') 'Enter the number of modes.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) SizeDist%n_mode
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (L_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!
!   Allocate modal arrays using this size.
    nd_mode = SizeDist%n_mode
    ALLOCATE(SizeDist%weight_mode(nd_mode))
    ALLOCATE(SizeDist%ln_r0_ln(nd_mode))
    ALLOCATE(SizeDist%ln_sigma_ln(nd_mode))
!
    SizeDist%dimen_0 = 0.0E+00_RealK
    DO i=1, SizeDist%n_mode
      WRITE(iu_stdout, '(/A, I3)') &
        'Enter the values of the weight, r_0 and sigma for mode', i
      DO
        READ(iu_stdin, *, IOSTAT=ios) SizeDist%weight_mode(i), &
                                      SizeDist%ln_r0_ln(i),    &
                                      SizeDist%ln_sigma_ln(i)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(A)') '+++ Erroneous input:'
          IF (l_interactive) THEN
            WRITE(iu_stdout, '(A)') 'Please Re-enter.'
          ELSE
            ierr = i_err_fatal
            RETURN
          ENDIF
        ELSE
          EXIT
        ENDIF
      ENDDO
!     SizeDist%ln_r0_ln(i) holds a size here.
      SizeDist%dimen_0 = SizeDist%dimen_0 &
        + SizeDist%weight_mode(i) * SizeDist%ln_r0_ln(i)
!     Since logs of r0 and sigma are needed we form these.
      SizeDist%ln_r0_ln(i)    = LOG(SizeDist%ln_r0_ln(i))
      SizeDist%ln_sigma_ln(i) = LOG(SizeDist%ln_sigma_ln(i))
    ENDDO
!
!
!
  ELSE IF (SizeDist%i_distribution == ip_modified_gamma) THEN
!
!   Spheres are assumed here.
    SizeDist%i_shape_particle=ip_shape_sphere
!
    WRITE(iu_stdout, '(/A)') 'Enter the total number density.'
    DO
      READ(iu_stdin, *, IOSTAT=ios)  number_tot
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO

    WRITE(iu_stdout, '(/A)') 'Enter the number of modes.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) SizeDist%n_mode
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!
!   Allocate modal arrays using this size.
    nd_mode = SizeDist%n_mode
    ALLOCATE(SizeDist%weight_mode(nd_mode))
    ALLOCATE(SizeDist%coeff_mg(nd_mode))
    ALLOCATE(SizeDist%alpha_mg(nd_mode))
    ALLOCATE(SizeDist%rm_mg(nd_mode))
    ALLOCATE(SizeDist%beta_mg(nd_mode))
!
    SizeDist%dimen_0=0.0E+00_RealK
    DO I=1, SizeDist%n_mode
      WRITE(iu_stdout, '(/A, /A, I3)') &
        'Enter the values of the weighting, the shape parameter, ' &
        //'the characteristic radius and the second power ' &
        , 'for mode ', i
      DO
        READ(iu_stdin, *, IOSTAT=ios) SizeDist%weight_mode(i), &
                                      SizeDist%alpha_mg(i), &
                                      SizeDist%rm_mg(i), &
                                      SizeDist%beta_mg(i)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(A)') '+++ Erroneous input:'
          IF (l_interactive) THEN
            WRITE(iu_stdout, '(A)') 'Please re-enter.'
          ELSE
            ierr = i_err_fatal
            RETURN
          ENDIF
        ELSE
          EXIT
        ENDIF
      ENDDO
!     The typical radius is set to the mode of the distribution.
      SizeDist%dimen_0 = SizeDist%dimen_0 &
                       + SizeDist%weight_mode(i) &
                       * SizeDist%rm_mg(i)
!     Precalculate the coefficient for future use.
      SizeDist%coeff_mg(i) = number_tot * SizeDist%beta_mg(i) / &
        (SizeDist%rm_mg(i) * &
         gamma_fnc(ierr, (SizeDist%alpha_mg(i)/SizeDist%beta_mg(i))) &
        )
    ENDDO
!
!
  ELSE IF (SizeDist%i_distribution == ip_heymsfield_platt) THEN
!
!   Hexagonal columns are assumed here.
    SizeDist%i_shape_particle = ip_shape_hexcyl
!
    WRITE(iu_stdout, '(/A)') &
      'Enter the smallest value of the maximum ' &
      //'dimension of the crystal.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) SizeDist%d_hp_min
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO

    WRITE(iu_stdout, '(/A)') &
      'Enter the number density at 100 microns and the first slope.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) SizeDist%a_hp_1, SizeDist%b_hp_1
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!
    WRITE(iu_stdout, '(/A)') &
      'Enter the number density at 1000 microns and the second slope.'
    DO
      READ(iu_stdin, *, IOSTAT=ios) SizeDist%a_hp_2, SizeDist%b_hp_2
      IF (ios /= 0) THEN
        WRITE(iu_err, '(A)') '+++ Erroneous input:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(A)') 'Please re-enter.'
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
!   Form the actual intercept parameters from the normalized
!   number densities and determine the transitional dimension.
    SizeDist%a_hp_1 = SizeDist%a_hp_1 &
                    / (1.0E-04_RealK**SizeDist%b_hp_1)
    SizeDist%a_hp_2 = SizeDist%a_hp_2 &
                    / (1.0E-03_RealK**SizeDist%b_hp_2)
    IF (ABS(SizeDist%b_hp_1-SizeDist%b_hp_2) > tol_zero) THEN
      SizeDist%d_hp_0 = EXP(LOG(SizeDist%a_hp_2/SizeDist%a_hp_1) &
                      / (SizeDist%b_hp_1-SizeDist%b_hp_2))
    ELSE
      SizeDist%d_hp_0 = 5.0E-04_RealK
    ENDIF
    IF ( (SizeDist%d_hp_0 < 1.0E-04_RealK) .OR. &
         (SizeDist%d_hp_0 > 1.0E-03_RealK) ) THEN
      WRITE(IU_ERR, '(/A)') &
        '*** Error: The transitional dimension is out of range.'
      ierr = i_err_fatal
      RETURN
    ENDIF
!   For want of a more suitable value of the typical dimension we take
!   100 microns.
    SizeDist%dimen_0=1.0E-04_RealK
!
  ELSE IF ( (SizeDist%i_distribution == ip_mitchell_97)      .OR. &
            (SizeDist%i_distribution == ip_mitchell_trop_00) .OR. &
            (SizeDist%i_distribution == ip_ivanova_mlat_00) ) THEN
!
!   These parameterizations apply to various shapes which must
!   be specified explicitly.
    DO 
      WRITE(iu_stdout, '(A)') 'Enter the shape of the crystals.'
      READ(iu_stdin, *) SizeDist%i_shape_particle
      IF ( (SizeDist%i_shape_particle == ip_shape_polycrystal) .OR. &
           (SizeDist%i_shape_particle == ip_shape_plate)       .OR. &
           (SizeDist%i_shape_particle == ip_shape_rosette)     .OR. &
           (SizeDist%i_shape_particle == ip_shape_column) ) THEN
        EXIT
      ELSE
        WRITE(iu_err, '(/A)') '*** Error: Illegal shape.'
        ierr = i_err_fatal
        RETURN
      ENDIF
    ENDDO
!
!   The distribution is entirely constrained by the mean maximum
!   dimension of the large mode, which serves to charactize it.
    WRITE(iu_stdout, '(A)') 'Enter the mean maximum dimension.'
    READ(iu_stdin, *) SizeDist%d_mean_max_l
!
!   These size distibutions are all bimodal gamma-distributions:
    IF (SizeDist%i_distribution == ip_mitchell_97) THEN
!
!     Small mode:
      SizeDist%nu_s=0.0E+00_RealK
      SizeDist%d_mean_max_l=0.0E+00_RealK
      IF (SizeDist%i_shape_particle == ip_shape_polycrystal) THEN
        SizeDist%d_mean_max_s = 40.0E-06_RealK
        kt                    = 2.5E+00_RealK
      ELSE IF ( (SizeDist%i_shape_particle == ip_shape_plate)   .OR. &
                (SizeDist%i_shape_particle == ip_shape_rosette) .OR. &
                (SizeDist%i_shape_particle == ip_shape_column) ) &
               THEN
        SizeDist%d_mean_max_s = 20.0E-06_RealK
        kt                    = 5.0E+00_RealK
      ENDIF
      SizeDist%lambda_s = (SizeDist%nu_s+1.0E+00_RealK) &
                         / SizeDist%d_mean_max_s
!
!     Large mode:
      SizeDist%nu_l = 2.0E+00_RealK * &
        (1.0E+00_RealK - EXP(-(SizeDist%d_mean_max_l &
                             /(kt*SizeDist%d_mean_max_s))**3))
      SizeDist%lambda_l = (SizeDist%nu_l + 1.0E+00_RealK) &
                         / SizeDist%d_mean_max_l
!
    ENDIF
!
  ENDIF
!
!
!
  RETURN
END

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to perform scattering calculations.
!
PROGRAM scatter_90
!
! Description:
!   This program calculates the monochromatic single scattering
!   properties of spherical particles averaged over a size distribution
!   at a range of specified wavelengths.
!
! Method:
!   A file of wavelengths for the calculation and a file of
!   refractive indices are read in. A distribution and a 
!   scattering algorithm are supplied. The program calculates
!   single scattering properties at each of the wavelengths given.
!
! IMPORTANT: This program is not yet complete in its treatment of
!            non-spherical particles.
!
!
! Modules used:
  USE realtype_rd
  USE dimensions_spec_ucf
  USE dimensions_pp_ucf
  USE parm_integ_acf
  USE def_std_io_icf
  USE scatter_algorithm_pcf
  USE shape_particle_pcf
  USE def_size_dist
  USE scatter_pp_pcf
  USE distribution_pcf
  USE rad_pcf
  USE rad_ccf, ONLY: pi
  USE file_type_pcf
  USE calc_gauss_weight_90_mod, ONLY: calc_gauss_weight_90
  USE spline_fit_mod, ONLY: spline_fit
!
!
  IMPLICIT NONE
!
!
!
!
  CHARACTER  (LEN=80) :: file_out*80
!   Name of output file
  CHARACTER  (LEN=80) :: file_phase*80
!   File containing phase function
  CHARACTER  (LEN=80) :: line*80
!   Line of input
  CHARACTER  (LEN=20) :: name_component*20
!   Name of scattering component
  CHARACTER  (LEN=1)  :: char_yn
!   Character response variable
  LOGICAL :: l_interactive
!   Flag for interactive use
!
!
  INTEGER :: iu_scatter
!   Unit number for output of scattering data
  INTEGER :: iu_phase
!   Unit number for output of the phase function
  INTEGER :: iu_file_in
!   Unit number for input file
  INTEGER :: ierr = i_normal
!   Error flag
  INTEGER :: ios
!       I/O error flag
  INTEGER :: n_wavelength
!       Number of wavelengths
  INTEGER :: n_refract
!       Number of pairs of refractive indices
  INTEGER :: n_refract_water
!       Number of pairs of refractive indices for water
  INTEGER :: n_humidities
!       Number of humidities
  INTEGER :: n_angle
!       Number of scattering angles
  INTEGER :: n_phf_term
!       Number of terms in the phase function
  INTEGER :: length_out
!       Length of name of output file
  INTEGER :: i_scatter_type
!       Index of scattering type
  INTEGER :: ichem_type
!       Chemical type of aerosol
  INTEGER :: i_scatter_algorithm
!       Indexing number of algorithm
  INTEGER :: i_component
!       Index of scattering component
  INTEGER :: i_output_type
!       Type of output file
  INTEGER :: i
!       Loop `index'
  INTEGER :: j
!       Loop `index'
  INTEGER :: k
!       Loop `index'
  LOGICAL :: l_data_region
!       Flag for input of data
  LOGICAL :: l_phase
!       Calculate phase function
  LOGICAL :: l_exist
!       File existence flag
  LOGICAL :: l_humid
!       Compute for varied humidities
  LOGICAL :: l_interp_spline
!       Use cubic splines to interpolate the refractive `index'
  LOGICAL :: l_stokes
!       Flag for the calculation of Stokes's first parameter
!
  TYPE (STR_size_dist) :: size_dist
!       Size distribution
!
  REAL  (RealK) :: size
!       Size parameter
  REAL  (RealK), Dimension(npd_wavelength_scat) :: wavelength
!       Wavelengths for scattering parameters
  REAL  (RealK), Dimension(npd_refract) :: wavelength_refract
!       Wavelengths of refractive indices
  REAL  (RealK), Dimension(npd_refract) :: wavelength_refract_water
!       Wavelengths of refractive indices of water
  REAL  (RealK), Dimension(npd_refract) :: re_refract
!       Real part of refractive `index'
  REAL  (RealK), Dimension(npd_refract) :: re_refract_water
!       Real part of water's refractive `index'
  REAL  (RealK), Dimension(npd_refract) :: im_refract
!       Imaginary part of refractive `index'
  REAL  (RealK), Dimension(npd_refract) :: im_refract_water
!       Imag. part of water's refractive `index'
  REAL  (RealK), Dimension(npd_refract) :: d2_re_refract
!       Second deriv. of re(refract)
  REAL  (RealK), Dimension(npd_refract) :: d2_re_refract_water
!       Second deriv. of re(refract) water
  REAL  (RealK), Dimension(npd_refract) :: d2_im_refract
!       Second deriv. of im(refract)
  REAL  (RealK), Dimension(npd_refract) :: d2_im_refract_water
!       Second deriv. of im(refract) water
  REAL  (RealK), Dimension(npd_humidities) :: humidity
!       Humidities
  REAL  (RealK) :: number_total
!       Number density of particles
  REAL  (RealK) :: radius_eff
!       Effective radius
  REAL  (RealK) :: radius_mean
!       Mean radius (no used here, but for consistency in argument list)
  REAL  (RealK) :: radius_typical
!       Typical radius
  REAL  (RealK) :: number_0
!       Number at peak radius
  REAL  (RealK) :: number_eff
!       Number at effective radius
  REAL  (RealK) :: volume_fraction
!       Volume fraction of spheres
  REAL  (RealK) :: proj_area_total
!       Overall projected area
  REAL  (RealK) :: number_total_dry
!       Number density of dry particles
  REAL  (RealK) :: radius_eff_dry
!       Dry effective radius
  REAL  (RealK) :: volume_fraction_dry
!       Volume fraction of dry spheres
  REAL  (RealK) :: growth_factor
!       Humidity growth factor.
  REAL  (RealK) :: dilution_factor
!       Humidity dilution factor.
!
! Parameters for analytic size distributions
  REAL  (RealK), Allocatable, Dimension(:) :: ln_r0_ln_dry
!       Dry value of ln_r0_ln
  REAL  (RealK), Allocatable, Dimension(:) :: ln_sigma_ln_dry
!       Dry value of ln_sigma_ln
  REAL  (RealK) :: radius_0
!       Typical radius of distribution
!
! Dimensions set by the size distribution
  INTEGER :: nd_size
!       Size allocated for points in the size distribution
  INTEGER :: nd_mode
!       Size allocated for modes of the size distribution
  INTEGER :: nd_integral
!       Size of array containing scattering integrals
  REAL  (RealK) :: q_scatter
!       Scattering efficiency
  REAL  (RealK) :: q_ext
!       Extinction efficiency
  REAL  (RealK), Dimension(npd_scatt_angle) :: angle
!       Scattering angles
  REAL  (RealK), Dimension(npd_scatt_angle) :: mu_angle
!       Cosines of scattering angles
  REAL  (RealK), Dimension(npd_scatt_angle) :: weight_angle
!       Gaussian weights applied to scattering angles
  REAL  (RealK) :: xs_scatter
!       Scattering cross-section
  REAL  (RealK) :: xs_ext
!       Extinction cross-section
  COMPLEX  (RealK), Dimension(npd_scatt_angle) :: s1
!       Matrix element
  COMPLEX  (RealK), Dimension(npd_scatt_angle) :: s2
!       Matrix element
  REAL  (RealK) :: absorption
!       Averaged absorption
  REAL  (RealK) :: extinction
!       Averaged extinction
  REAL  (RealK) :: scattering
!       Averaged scattering
  REAL  (RealK), Dimension(npd_scatt_angle) :: i_stokes
!       Stokes's first parameter at this angle
  REAL  (RealK), Dimension(npd_phase_term) :: phase_fnc_term
!       Averaged moments of the phase function
  COMPLEX  (RealK) :: n_relative
!       Relative `index' of refraction
  COMPLEX  (RealK) :: n_relative_water
!       Refractive `index' of water
!
! Subroutines called:
  EXTERNAL &
    get_free_unit, open_file_in, &
    particle_size_90, size_integral_90, grow_particles, &
    mie_scatter, scatter_integral_90, adt_integral, &
    decompose_phf_90
!
! Functions called:
  COMPLEX  (RealK) :: refractive_index
!   Function to find refractive `index'
  REAL  (RealK) :: number_particle_90
!   Function for number density
  LOGICAL :: set_interactive
!   Function to determine whether operation is interactive
  EXTERNAL &
      set_interactive, refractive_index, number_particle_90
!
!
!
! Set the flag for interactive operation
  l_interactive=set_interactive()
!
! Set a unit number for generic input
  CALL get_free_unit(ierr, iu_file_in)
  IF (ierr /= i_normal) STOP
!
!
! Obtain the wavelengths where scattering calculations are required.
  CALL get_wavelengths &
    (npd_wavelength_scat, n_wavelength, wavelength, ierr)
!
!
! Obtain the relevant refractive indices.
  CALL get_refract_index(l_interactive, npd_refract, n_refract, &
    wavelength_refract, re_refract, im_refract, &
    l_interp_spline, d2_re_refract, d2_im_refract, &
    ierr)
  IF (ierr /= i_normal) STOP
!
! Ask whether humidity effects are required and set the type of
! output file.
!
  WRITE(iu_stdout, '(/a)') 'Is growth with humidity required?'
  DO
    READ(iu_stdin, '(a)') char_yn
    IF ( (char_yn == 'Y') .OR. (char_yn.eq.'y') ) THEN
      l_humid =  .TRUE. 
      WRITE(iu_stdout, '(a)') &
        'A list of the humidities must also be supplied.'
      EXIT
    ELSE IF ( (char_yn == 'N') .OR. (char_yn.eq.'n') ) THEN
      l_humid =  .FALSE. 
      n_humidities = 1
      EXIT
    ELSE
      WRITE(iu_err, '(/a)') '+++ unrecognized response:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a/)') 'please re-enter.'
      ELSE
        STOP
      ENDIF
    ENDIF
  ENDDO
!
  WRITE(iu_stdout, '(/)')
!
! Specification of refractive indices for water.
! This is only required if humidity effects are being taken
! into account.
!
  IF (l_humid) THEN
    CALL open_file_in(ierr, iu_file_in, &
      'Specify the file containing refractive indices of water.')
    IF (ierr /= i_normal) stop
    l_data_region   = .FALSE. 
    n_refract_water = 0
    DO
      READ(iu_file_in, '(a)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (l_data_region) THEN
        IF (line(1:4) /= '*END') THEN
          backspace(iu_file_in)
          n_refract_water = n_refract_water + 1
          IF (n_refract_water > npd_refract) THEN
            WRITE(iu_err, '(/a)') '*** Error: The file of ' // &
              'refractive indices for water is too long.'
            ierr=i_err_fatal
            STOP
          ENDIF
          READ(iu_file_in, *) &
            wavelength_refract_water(n_refract_water), &
            re_refract_water(n_refract_water), &
            im_refract_water(n_refract_water)
        ELSE
          l_data_region= .FALSE. 
        ENDIF
      ELSE
        IF (line(1:11) == '*BEGIN_DATA') THEN
          l_data_region= .TRUE. 
        ENDIF
      ENDIF
    ENDDO
!
    CLOSE(iu_file_in)
!
!   Specification of humidities.
!
    CALL open_file_in(ierr, iu_file_in, &
      'Specify the file of relative humidities.')
    IF (ierr /= i_normal) STOP
    l_data_region = .FALSE. 
    n_humidities  = 0
    DO
      READ(iu_file_in, '(a)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (l_data_region) THEN
        IF (line(1:4) /= '*END') THEN
          backspace(iu_file_in)
          n_humidities = n_humidities + 1
          IF (n_humidities > npd_humidities) THEN
            WRITE(iu_err, '(/a)')  &
              '*** Error: There are too many humidities: ' // &
              'increase npd_humidities and recompile.'
            ierr=i_err_fatal
            STOP
          ENDIF
          READ(iu_file_in, *) humidity(n_humidities)
        ELSE
          l_data_region = .FALSE. 
        ENDIF
      ELSE
        IF (line(1:11) == '*BEGIN_DATA') THEN
          l_data_region= .TRUE. 
        ENDIF
      ENDIF
    ENDDO
!
    CLOSE(iu_file_in)
!
!
!   Find the second derivatives for the spline fits to the
!   refractive index of water.
!
    CALL spline_fit(n_refract_water, wavelength_refract_water, &
      re_refract_water, d2_re_refract_water)
    CALL spline_fit(n_refract_water, wavelength_refract_water, &
      im_refract_water, d2_im_refract_water)
!
!   End of block of code for humidity effects only.
!
  ENDIF
!
  IF (l_interp_spline) THEN
!
!   Find the second derivatives for the spline fits to the 
!   refractive `index'.
    CALL spline_fit(n_refract, wavelength_refract, &
      re_refract, d2_re_refract)
    CALL spline_fit(n_refract, wavelength_refract, &
      im_refract, d2_im_refract)
!
  ELSE
!
!   If cubic splines were not requested, zero the second
!   derivatives to obtain linear interpolation instead.
    d2_re_refract(:) = 0.0_RealK
    d2_im_refract(:) = 0.0_RealK
!
  ENDIF
!
! Set the size distribution.
! If humidity variation is being used, this is the size distribution
! of the dry aerosol, before humidity growth takes place.
!
  CALL particle_size_90(l_interactive, nd_size, nd_mode, &
    size_dist, ierr)
  IF (ierr /= i_normal) THEN
    ierr=i_err_fatal
    stop
  ENDIF
!
! Check that the type of distribution is consistent with the
! assumptions made in the program.
!
  IF (l_humid .AND. &
       (size_dist%i_distribution /= ip_log_normal) ) THEN
    WRITE(iu_stdout,'(/a)') 'Humidity effects are only available ' // &
      'with log-normal size distributions at present.'
    STOP
  ENDIF
  IF (l_humid .AND. &
       (size_dist%n_mode /= 1) ) THEN
    WRITE(iu_stdout,'(/a)') 'Humidity effects are only available ' // &
      'with a monomodal size distribution at present.'
    STOP
  ENDIF
!
! Select the scattering algorithm.
  WRITE(iu_stdout, '(/a)') &
    'Select the algorithm for the scattering calculation.'
  DO i=1, npd_algorithm_scatter
    WRITE(iu_stdout, '(3x, i3, a3, a40)') &
      i, '.) ', name_scatter_algorithm(i)
  ENDDO
  WRITE(iu_stdout, '(/)')
  DO
    READ(iu_stdin, *) i_scatter_algorithm
!
!   ADT is a permitted option only with a gamma distribution.
    IF (i_scatter_algorithm == ip_algorithm_adt) THEN
      IF (size_dist%i_distribution /= ip_modified_gamma) THEN
        WRITE(iu_err, '(/a)') &
          '+++ ADA is permitted only with a modified ' // &
          'gamma-distribution.'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(a/)') 'Please re-enter.'
        ELSE
          STOP
        ENDIF
      ELSE
        EXIT
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
! Determine the consistency of the scattering algorithm with the
! shape of the particle.
  IF ( ( (i_scatter_algorithm == ip_algorithm_mie) .AND.  &
         (size_dist%i_shape_particle /= ip_shape_sphere) ) .OR.  &
       ( (i_scatter_algorithm == ip_algorithm_adt) .AND.  &
         (size_dist%i_shape_particle /= ip_shape_sphere) ) ) THEN
    WRITE(iu_err, '(/a)') &
      '*** Eror: the scattering algorithm selected ' // &
      'is not applicable to this shape.'
    STOP
  ENDIF
!
!
!
! Determine whether the phase function is to be calculated. This is
! only permitted for an external distribution with one size and
! only for certain algorithms.
!
! Note also that this refers to an explicit respresentation of the
! phase function at specific angles, not to the specification of 
! moments.
!
  IF ( (size_dist%i_distribution == ip_external) .AND. &
       (size_dist%n_size == 1) .AND. &
       (i_scatter_algorithm == ip_algorithm_mie) ) THEN
    WRITE(iu_stdout, '(/a)') &
      'Is the phase function to be calculated ? (y/n)'
    READ(iu_stdin, '(a)') char_yn
    IF ( (char_yn == 'Y') .OR. (char_yn.eq.'y') ) THEN
      l_phase= .TRUE. 
      CALL open_file_in(ierr, iu_file_in, &
        'Specify the file containing the scattering angles.')
      IF (ierr /= i_normal) STOP
      l_data_region = .FALSE. 
      n_angle       = 0
      DO
        READ(iu_file_in, '(a)', IOSTAT=ios) line
        IF (ios /= 0) EXIT
        IF (l_data_region) THEN
          IF (line(1:4) /= '*END') THEN
            BACKSPACE(iu_file_in)
            n_angle = n_angle + 1
            IF (n_angle > npd_scatt_angle) THEN
              WRITE(iu_err, '(/a)') '*** Error: There are ' // &
                'too many angles in the input file.'
              ierr=i_err_fatal
              STOP
            ENDIF
            READ(iu_file_in, *) angle(n_angle)
          ELSE
            l_data_region= .FALSE. 
          ENDIF
        ELSE
          IF (line(1:11) == '*BEGIN_DATA') THEN
            l_data_region= .TRUE. 
          ENDIF
        ENDIF
      ENDDO
!
      CLOSE(iu_file_in)
!
      DO j=1, n_angle
        mu_angle(j)=cos(angle(j))
      ENDDO
    ELSE
      l_phase= .FALSE. 
    ENDIF
  ELSE
     l_phase= .FALSE. 
  ENDIF
!
!
! Determine the number of moments of the phase function required.
  WRITE(iu_stdout, '(/a)') &
    'Enter the number of moments of the phase function.'
  READ(iu_stdin, *, iostat=ios) n_phf_term
  IF (ios /= 0) THEN
    WRITE(iu_err, '(/a)') &
      '*** Error: This number could not be read.'
    ierr=i_err_io
    STOP
  ENDIF
!
! If only one moment of the phase function (the asymmetry) is
! required it is calculated directly, but if other moments are
! needed we calculate the phase function and use Gaussian 
! quadrature at order order above the term retained: given the
! properties of Gausian quadrature this will be sufficient.
! We need the array of angles for special purposes in this
! case, so the arrays S1 and S2 at arbitrary angles cannot be
! calculated: for now we forid the calculation of the phase
! function in full generality in this case. To calculate these
! moments we need the Stokes's first parameter.
  IF (n_phf_term > 1) THEN
    IF (l_phase) THEN
      WRITE(iu_err, '(/a, /a)') &
        '*** Error: The simultaneous calculation of a full', &
        'phase function and its moments is not yet allowed.'
      STOP
    ENDIF
    l_stokes= .TRUE. 
!
!   If too small a value of n_angle is used the quadrature becomes
!   unreliable. "9" is a heuristic lower limit based a single test.
!   I thought Gaussian quadrature was supposed to do better than this.
    n_angle=MAX(9, n_phf_term+1)
!
    IF (n_angle > npd_scatt_angle) THEN
      WRITE(iu_err, '(/a, /a,/a)') &
        '*** Error: The allowed size of arrays of scattering', &
        'angles is too small: recompile with a larger value of', &
        'npd_scatt_angle'
      stop
    ENDIF
    CALL calc_gauss_weight_90(ierr, n_angle, mu_angle, weight_angle)
  ELSE
    l_stokes= .FALSE. 
  ENDIF
!
! Set the size of the dummy array of integrands.
  IF (n_phf_term > 1) THEN
    nd_integral=n_angle+3
  ELSE
    nd_integral=n_phf_term+2
  ENDIF
!
!
! Read in the name of the file for output.
  WRITE(iu_stdout, '(/a/)') &
    'Give the name of the file to contain the scattering data.'
  READ(iu_stdin, '(a)', iostat=ios) file_out
  IF (ios /= 0) THEN
    WRITE(iu_err, '(/a)') &
      '*** error: this name could not be read.'
    ierr=i_err_io
    STOP
  ENDIF
! Check for the existence of the file.
  DO
    INQUIRE(FILE=file_out, EXIST=l_exist)
    IF (l_exist) THEN
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/a)') 'This file already exists.'
        WRITE(iu_stdout, '(/a/)') 'Do you wish to overwrite? (y/n).'
        READ(iu_stdin, '(a)') char_yn
        IF ( (char_yn /= 'Y') .AND. (char_yn /= 'y') ) THEN
          WRITE(iu_stdout, '(/a/)') &
            'Please specify another name or !quit.'
          READ(iu_stdin, '(a)') file_out
          IF ( (file_out(1:5) == '!quit') .OR. &
               (file_out(1:5) == '!Quit') .OR. &
               (file_out(1:5) == '!QUIT') ) THEN
            WRITE(iu_stdout, '(/a/)') '*** program terminated.'
            ierr=i_err_io
            STOP
          ENDIF
        ELSE
          EXIT
        ENDIF
      ELSE
        STOP
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
! The file may now safely be opened.

  CALL get_free_unit(ierr, iu_scatter)

  OPEN(UNIT=iu_scatter, FILE=file_out, IOSTAT=ios, STATUS='UNKNOWN')
  IF (ios /= 0) THEN
    WRITE(iu_err, '(/, 3a)') &
      '*** Error: File ', file_out, ' could not be opened.'
    ierr=i_err_io
    STOP
  ENDIF
!
! We now have enough information to determine the type of output
! file.
  IF ( (i_scatter_algorithm == ip_algorithm_mie) ) THEN
    IF (l_humid) THEN
      IF (n_phf_term > 1) THEN
        i_output_type=it_file_phf_mie_humid
      ELSE
        i_output_type=it_file_mie_humid
      ENDIF
    ELSE
      IF (n_phf_term > 1) THEN
        i_output_type=it_file_phf_mie_dry
      ELSE
        i_output_type=it_file_mie_dry
      ENDIF
    ENDIF
  ELSE
    IF (l_humid) THEN
      i_output_type=it_file_adt_humid
    ELSE
      i_output_type=it_file_adt_dry
    ENDIF
  ENDIF
!
! Open a file for the phase function.
  IF (l_phase) THEN
!   Find the length of the output filename.
    length_out=LEN(file_out)
    DO
      IF (file_out(length_out: length_out) /= ' ') EXIT
      length_out=length_out-1
    ENDDO
    file_phase(1: length_out+6)=file_out(1: length_out)//'.phase'
    WRITE(iu_stdout, '(/a)') &
      'The phase function will be written to the file:'
    WRITE(iu_stdout, '(6x, a)') file_phase
    OPEN(UNIT=iu_phase, FILE=file_phase, IOSTAT=ios, STATUS='unknown')
    IF (ios /= 0) THEN
      WRITE(iu_err, '(/, 3a)') &
        '*** Error: The file ', file_phase, ' could not be opened.'
      ierr=i_err_io
      STOP
    ENDIF
  ENDIF
!
! Set the component indices for the file.
  WRITE(iu_stdout, '(/a/)') &
    'Enter the type number for the scatterer.'
  DO i=0, npd_scatter_type
    WRITE(iu_stdout, '(3x, i3, a2, 3x, a20)') &
      i, '.)', name_scatter_type(i)
  ENDDO
  WRITE(iu_stdout, '(/)')
!
  DO
    READ(iu_stdin, *, iostat=ios) i_scatter_type
    IF (ios /= 0) THEN
      WRITE(iu_err, '(a)') '+++ illegal response:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'please re-type.'
      ELSE
        STOP
      ENDIF
    ELSE
      IF ( (i_scatter_type < 0) .OR. &
           (i_scatter_type > npd_scatter_type) ) THEN
        WRITE(iu_err, '(a)') '+++ The response is out of range: '
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'please retype.'
        ELSE
          STOP
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDIF
  ENDDO
!
  IF (i_scatter_type == ip_type_unassigned) THEN
!   Type 0 is unassigned to allow use of the program as a pure Mie
!   code for different substances when there is no intention of
!   proceeding to assemble a spectral file.
    i_component=0
    name_component='unassigned          '
  ELSE IF (i_scatter_type == ip_type_aerosol) THEN
    WRITE(iu_stdout, '(/a)') 'Specify the aerosol component:'
    DO i=1, npd_aerosol_component
      WRITE(iu_stdout, '(3x, i3, a2, 3x, a20)') &
        i, '.)', name_aerosol_component(i)
    ENDDO
    WRITE(iu_stdout, '(/)')
!
    DO
      READ(iu_stdin, *, iostat=ios) i_component
      IF ( (ios /= 0) .OR. &
           (i_component > npd_aerosol_component) ) THEN
        WRITE(iu_err, '(a)') '+++ illegal response:'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-type.'
        ELSE
          STOP
        ENDIF
      ELSE
        IF ( (i_component <= 0) .OR. &
             (i_component > npd_aerosol_component) ) THEN
          WRITE(iu_err, '(a)') '+++ The response is out of range: '
          IF (l_interactive) THEN
            WRITE(iu_stdout, '(a)') 'Please re-type.'
          ELSE
            STOP
          ENDIF
        ELSE
          name_component=name_aerosol_component(i_component)
          EXIT
        ENDIF
      ENDIF
    ENDDO
!
  ENDIF
!
!
! Calculate the effective radius and volume fraction.
! When the effects of humidity are included these refer to the
! dry distribution.
!
  IF ( (size_dist%i_distribution == ip_external) .AND. &
       (size_dist%n_size == 1) ) THEN
!
    IF (size_dist%i_shape_particle == ip_shape_sphere) THEN
!
      proj_area_total=pi * size_dist%number(size_dist%n_size) * &
        size_dist%dimen(size_dist%n_size)**2
      volume_fraction=4.0_RealK*pi * &
        size_dist%number(size_dist%n_size) * &
        size_dist%dimen(size_dist%n_size)**3/3.0_RealK
      number_total=size_dist%number(size_dist%n_size)
!
    ELSE IF (size_dist%i_shape_particle == ip_shape_hexcyl) THEN
!
      proj_area_total=(3.0_RealK/2.0_RealK) * &
        size_dist%number(size_dist%n_size) * &
        size_dist%dimen(size_dist%n_size)**2 * &
        (size_dist%aspect(size_dist%n_size) + &
        SQRT(3.0_RealK)/4.0_RealK) 
      volume_fraction=(SQRT(2.7e+01_RealK) / 4.0_RealK) * &
        size_dist%number(size_dist%n_size) * &
        size_dist%dimen(size_dist%n_size)**3 * &
        size_dist%aspect(size_dist%n_size)
      number_total=size_dist%number(size_dist%n_size)
!
    ENDIF
!
  ELSE
!
    CALL size_integral_90(nd_integral, npd_panel, npd_refinement, &
      npd_point, npd_panel_point, p_panel_ratio, &
      size_dist, &
      number_total, proj_area_total, volume_fraction, radius_mean, &
      ierr)
    IF (ierr /= i_normal) STOP
!
  ENDIF
!
  radius_eff= (3.0_realK/4.0_RealK) * volume_fraction / &
    proj_area_total
!
  IF (l_phase) THEN
    WRITE(iu_phase, '(27x, a///)') 'Scattering phase functions'
  ENDIF
!
! Store the parameters of the dry distribution here.
! At present it is assumed that a lognormal distribution is used.
!
  IF (l_humid) THEN
!
    IF (size_dist%i_distribution /= ip_log_normal) THEN
      WRITE(iu_err, '(/a)') &
        '*** Error: At present only a log-normal' // &
        'distribution may be used for moist aerosols.'
      STOP
    ENDIF
!
    allocate(ln_r0_ln_dry(nd_mode))
    allocate(ln_sigma_ln_dry(nd_mode))
    ln_r0_ln_dry = size_dist%ln_r0_ln
    ln_sigma_ln_dry = size_dist%ln_sigma_ln
!
    volume_fraction_dry=volume_fraction
    radius_eff_dry=radius_eff
    number_total_dry=number_total
!
  ENDIF
!
!
! Begin loop over humidity values.
!
  Loop_Humidity: DO k = 1, n_humidities
!
    IF (l_humid) THEN
!
!     Adjust size of particles to allow for the effect of 
!     humidity. Note that it is assumed that both the dry and 
!     moist size distributions are lognormal.
!
      ichem_type = i_component
      CALL grow_particles(humidity(k),size_dist%ln_r0_ln, &
        size_dist%ln_sigma_ln, &
        ichem_type,growth_factor,humidity(n_humidities-1),ierr)
      IF (ierr /= i_normal) STOP
!
      CALL size_integral_90(nd_integral, npd_panel, npd_refinement, &
        npd_point, npd_panel_point, p_panel_ratio, size_dist, &
        number_total, proj_area_total, volume_fraction, radius_mean, &
        ierr)
      IF (ierr /= i_normal) STOP
!
    ENDIF
!
    radius_eff= (3.0_realK/4.0_RealK) * volume_fraction / &
      proj_area_total
!
!   Write the headers to the file
!
    WRITE(iu_scatter, '(/, a13, i5)') '*FILE TYPE = ', i_output_type
    WRITE(iu_scatter, '(/, 16x, a, //)') &
      'Single scattering parameters.'
    WRITE(iu_scatter, '(a23, i5, 3x, a1, a40, a1)') &
      'Scattering algorithm = ', i_scatter_algorithm, &
      '(', name_scatter_algorithm(i_scatter_algorithm), ')'
    WRITE(iu_scatter, '(a23, i5, 3x, a1, a20, a1)') &
      'Scattering type      = ', i_scatter_type, &
      '(', name_scatter_type(i_scatter_type), ')'
    IF (i_scatter_type == ip_type_aerosol) THEN
      WRITE(iu_scatter, '(a23, i5, 3x, a1, a20, a1)') &
        'Component            = ', i_component, &
        '(', name_component, ')'
    ENDIF
    IF (l_humid) THEN
      WRITE(iu_scatter,'(/,a11,f7.5)') 'humidity = ',humidity(k)
    ENDIF
    WRITE(iu_scatter, '(/, a18, 1x, 1pe12.5, 1x, a3)') &
      'Number density   =', number_total, 'm-3'
    WRITE(iu_scatter, '(a18, 1x, 1pe12.5, 1x, a1)') &
      'Effective radius =', radius_eff, 'm'
    IF (l_humid) THEN
      WRITE(iu_scatter, '(a35, 1x, 1pe12.5, 1x, a1)') &
        'Effective radius of dry component =', &
        radius_eff_dry, 'm'
    ENDIF
    WRITE(iu_scatter, '(a18, 1x, 1pe12.5)') &
      'Volume fraction  =', volume_fraction
    IF (l_humid) THEN
      WRITE(iu_scatter, '(a35, 1x, 1pe12.5)') &
        'Volume fraction of dry component  =', &
        volume_fraction_dry
    ENDIF
    IF (n_phf_term > 1) THEN
      WRITE(iu_scatter, '(/, a40, 1x, i3, /)') &
        'Number of terms in the phase function =', n_phf_term
    ENDIF
    WRITE(iu_scatter, '(/, 4x, a14, 6x, a16, 4x, a16, 4x, a14)') &
      'Wavelength (m)', 'Absorption (m-1)', &
      'Scattering (m-1)', 'Moments of pf'
!
!
!   Loop over wavelengths.
!
    DO i=1, n_wavelength
!
!     Find the refractive `index' for the current wavelength.
!
      n_relative=refractive_index(ierr, npd_refract, wavelength(i), &
        n_refract, wavelength_refract, &
        re_refract, im_refract, &
        d2_re_refract, d2_im_refract )
!
!     For moist aerosol, the effect of dilution is included.
!     it is assumed that all particles grow by the same factor.
!
      IF (l_humid) THEN
        n_relative_water = refractive_index(ierr, &
          npd_refract, wavelength(i), &
          n_refract_water, wavelength_refract_water, &
          re_refract_water, im_refract_water, &
          d2_re_refract_water, d2_im_refract_water )
!
!       Dilution produces a volume-weighted average of the
!       complex refractive `index'.
!
        dilution_factor = 1.0_RealK/(growth_factor)**3
        n_relative = n_relative*dilution_factor + &
          n_relative_water*(1.0_RealK - dilution_factor)
      ENDIF
!
!
      IF ( (size_dist%i_distribution == ip_external) .AND. &
           (size_dist%n_size == 1) ) THEN
        size=2.0_RealK*pi*radius_eff/wavelength(i)
        IF (i_scatter_algorithm == ip_algorithm_mie) THEN
          CALL mie_scatter(ierr, size, n_relative, &
            q_scatter, q_ext, phase_fnc_term(1), &
            l_phase.OR.l_stokes, n_angle, mu_angle, s1, s2, &
            npd_scatt_angle, npd_log_deriv)
          DO j=1, n_angle
            i_stokes(j) = &
              size_dist%number(size_dist%n_size) * &
              (abs(s1(j))**2+abs(s2(j))**2)
          ENDDO
          IF (ierr /= i_normal) THEN
            ierr=i_err_fatal
            STOP
          ENDIF
          extinction=pi*size_dist%number(size_dist%n_size) * &
            q_ext*size_dist%dimen(size_dist%n_size)**2
          scattering=pi*size_dist%number(size_dist%n_size) * &
            q_scatter*size_dist%dimen(size_dist%n_size)**2
        ENDIF
      ELSE
!       Notice that the effective radius is now normally supplied
!       as the typical value: this should be more efficient. 
!       However, in the case of non-spherical particles the 
!       distribution may be 0 at the effective radius.
        IF (size_dist%i_shape_particle == ip_shape_sphere) THEN
          radius_typical=radius_eff
        ELSE
          number_0=number_particle_90(radius_0, size_dist, ierr)
          number_eff=number_particle_90(radius_eff, size_dist, ierr)
          IF (number_eff*radius_eff**2 > &
               1.0e-02_RealK*number_0*radius_0**2) THEN
            radius_typical=radius_eff
          ELSE
            radius_typical=radius_0
          ENDIF
        ENDIF
        IF (i_scatter_algorithm == ip_algorithm_mie) THEN
          CALL scatter_integral_90(size_dist, &
             wavelength(i), n_relative, &
             i_scatter_algorithm, &
             n_angle, mu_angle, &
             p_panel_ratio, &
             extinction, scattering, phase_fnc_term(1), &
             l_stokes, i_stokes, &
             npd_scatt_angle, npd_log_deriv, nd_integral, &
             npd_panel, npd_refinement, &
             npd_point, npd_panel_point, &
             ierr)
        ELSE IF (i_scatter_algorithm == ip_algorithm_adt) THEN
          extinction=0.0_RealK
          scattering=0.0_RealK
          DO j=1, size_dist%n_mode
            size=2.0_RealK*pi * &
              size_dist%rm_mg(j)**size_dist%beta_mg(j) / &
              wavelength(i)
            CALL adt_integral(size_dist%alpha_mg(j), size, &
              n_relative, q_ext, q_scatter)
            extinction=extinction+size_dist%weight_mode(j) + &
              q_ext*7.5e-01_RealK*volume_fraction/radius_eff
            scattering=scattering+size_dist%weight_mode(j) + &
              q_scatter*7.5e-01_RealK*volume_fraction/radius_eff
          ENDDO
!         ADA does not give an 
!         asymmetry factor: zero is printed out as a flag.
          phase_fnc_term(1)=0.0_RealK
        ENDIF
        IF (ierr /= i_normal) THEN
          ierr=i_err_fatal
          STOP
        ENDIF
      ENDIF
!
      IF (l_phase) THEN
        WRITE(iu_phase, '(/a13, 1pe12.5, 6x, a17, 1pe12.5)') &
          'Wavelength = ', wavelength(i), 'Size parameter = ', size
        WRITE(iu_phase, '(6x, a19, 2(1x, 1pe12.5), /)') &
          'Refractive index = ', n_relative
        WRITE(iu_phase, '(9x, a5, 8x, a6, 16x, a2, 23x, a2)') &
          'Angle', 'Cosine', 's1', 's2'
        DO j=1, n_angle
          WRITE(iu_phase, '(6x, 1pe10.3, 4x, 1pe10.3, 3x, ' // &
            '2(1x, 1pe10.3), 3x, 2(1x, 1pe10.3))') &
            angle(j), mu_angle(j), s1(j), s2(j)
        ENDDO
        WRITE(iu_phase, '(/)')
      ENDIF
!
      absorption=extinction-scattering
!
      IF (n_phf_term > 1) THEN
!       Decompose the phase function into its moments.
        CALL decompose_phf_90(n_angle, mu_angle, weight_angle, &
          n_phf_term, i_stokes, &
          phase_fnc_term, ierr)
      ENDIF
!
      WRITE (iu_scatter, '(4(4x, 1pe16.9), :, /, (20x, 3(4x, 1pe16.9)))') &
        wavelength(i), absorption, scattering, &
        (phase_fnc_term(j), j=1, n_phf_term)
!
    ENDDO
!
    WRITE(iu_scatter, '(1x)')
!
!   Remember to restore the dry values of the distribution
!   parameters before going around the loop again (or indeed
!   leaving the loop.)
!
    IF (l_humid) THEN
      size_dist%ln_r0_ln = ln_r0_ln_dry
      size_dist%ln_sigma_ln = ln_sigma_ln_dry
    ENDIF
!
!   End of humidity loop.
!
  ENDDO Loop_Humidity
!
!
  CLOSE(iu_scatter)
!
!
!
  STOP
END PROGRAM scatter_90

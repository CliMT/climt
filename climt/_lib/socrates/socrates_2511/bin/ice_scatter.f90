! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to calculate monochromatic scattering for ice crystals
!
PROGRAM Ice_Scatter
!
!
! Description:
!   This program calculates the monochromatic single scattering 
!   properties of ice crystals averaged over a size distribution
!   at a range of specified wavelengths.
!
! Method:
!   Various algorithms are provided.
!
!
! Modules used:
!
  USE realtype_rd
  USE error_pcf
  USE def_std_io_icf
  USE scatter_pp_pcf
  USE scatter_algorithm_pcf
  USE shape_particle_pcf
  USE measure_particle_pcf
  USE distribution_pcf
  USE def_size_dist
  USE rad_ccf, ONLY: pi
  USE file_type_pcf
  USE dimensions_pp_ucf
  USE dimensions_spec_ucf
  USE parm_integ_acf
  USE def_sct_db
  USE def_db_crystal_geometry
  USE def_db_ss_mono
  USE calc_gauss_weight_90_mod, ONLY: calc_gauss_weight_90
!
!
  IMPLICIT NONE
!

  INTEGER :: ierr
!              Error flag
  INTEGER :: ios
!              I/O error flag
  INTEGER :: i_scatter_algorithm
!              Algorithm used to calculate scattering properties
  INTEGER :: n_phf_term
!              Number of terms in the phase function
  INTEGER :: nd_integral
!              Size allocated for integrals
  INTEGER :: iu_db_input
!              Unit number for reading the scattering database
  INTEGER :: iu_scatter
!              Unit number for output of scattering data
  INTEGER :: i_output_type
!              Type number of output file
  INTEGER :: i_scatter_type
!              Type of scattering particle
  INTEGER :: n_wavelength
!              Number of wavelengths
  INTEGER :: n_refract
!              Number of wavlengths for refractive indices
  INTEGER :: nd_size
!       Size allocated for points in the size distribution
  INTEGER :: nd_mode
!       Size allocated for number of modes
  INTEGER :: i
!              Loop variable
  INTEGER :: j
!              Loop variable
  LOGICAL :: l_interactive
!              Flag for interactive use
  LOGICAL :: l_phase
!              Flag to calculate the phase function
  LOGICAL :: l_exist
!              Flag for existence of file
  LOGICAL :: l_spline_interp
!              Flag for interpolation using cubic splines
  LOGICAL :: l_stokes
!              Flag to calculate Stokes parameters
  LOGICAL :: l_matched
!              Flag for checking that the scattering algforithm
!              is consistent with the shape of the particle
  CHARACTER (LEN=80) :: file_out
!                         Name of output file
  CHARACTER (LEN=1)  :: char_yn
!                         Character response variable
  TYPE (STR_size_dist) :: SizeDist
!           Size distribution
!
  CHARACTER  (LEN=80) :: ice_bin_database
!   Filename of the binary database
!
  TYPE (STR_db_cryst_geom) :: DBGeom
!           Geometry for specifying crystals in a database
  REAL (RealK) :: dm_average
!           Mean maximum dimension
!
  TYPE (STR_db_ss_mono) :: ice_db_mono_info
!           Total monochromatic single scattering information
!           at a single wavelength
!
  REAL (RealK) :: volume_total
!                   Total volume of the size distribution
  REAL (RealK) :: number_total
!                   Total number density of the size distribution
  REAL (RealK) :: proj_area_total
!                   Total projected area of the size distribution
!
  REAL (RealK) :: weight_angle(npd_scatt_angle)
!           Gaussian weights applied to scattering angles
  COMPLEX (RealK) :: n_relative
!                      Relative refractive index
  REAL (RealK) :: size
!                   Size parameter of the particle
  INTEGER :: n_angle
!              Number of angles where the phase function is 
!              calculated
!
! General variables for the database
  INTEGER, Dimension(npd_wavelength_scat) :: n_rec_block
!   Number of blocks of data for each wavelength
  INTEGER, Dimension(npd_wavelength_scat, npd_size_scat) :: db_record
!   Pointers to the records for each wavelength and each block at
!   that wavelength
  INTEGER :: db_n_angle
!   Common number of angles in the database

  REAL (RealK) :: q_ext
!                   Extinction efficiency
  REAL (RealK) :: q_scat
!                   Scattering efficiency
  REAL (RealK) :: extinction
!                   Actual extinction
  REAL (RealK) :: scattering
!                   Actual scattering
  REAL (RealK) :: absorption
!                   Actual absorption
  REAL (RealK) :: asymmetry
!                   Actual asymmetry
  REAL (RealK) :: i_stokes(npd_scatt_angle)
!                   Stokes''s first parameter at the specified angles
  REAL (RealK) :: phase_fnc_term(npd_phase_term)
!           Averaged moments of the phase function
  REAL (RealK) :: mu_angle(npd_scatt_angle)
!           Cosines of scattering angles
  COMPLEX (RealK) :: s1(npd_scatt_angle)
!                      Element of the scattering matrix
  COMPLEX (RealK) :: s2(npd_scatt_angle)
!                      Element of the scattering matrix
!
!
  REAL (RealK) :: wavelength(npd_wavelength_scat)
!                   Wavelengths for refractive indices
  REAL (RealK) :: wavelength_refract(npd_refract)
!                   Wavelengths for refractive indices
  REAL (RealK) :: re_refract(npd_refract)
!                   Real parts of the refractive indices
  REAL (RealK) :: im_refract(npd_refract)
!                   Imaginary parts of the refractive
!                   indices
  REAL (RealK) :: d2_re_refract(npd_refract)
!                   Second derivatives of the real parts
!                   of the refractive indices
  REAL (RealK) :: d2_im_refract(npd_refract)
!                   Second derivatives of the imaginary
!                   parts of the refractive indices
!
  LOGICAL :: set_interactive
!              Function to determine whether operation is interactive
!
  EXTERNAL set_interactive
!
  COMPLEX (RealK) :: refractive_index
!                   Function to calculate refractive indices
!
  EXTERNAL refractive_index
!
  REAL (RealK) :: proj_area_particle
!                   Function to calculate the projected area 
!                   of a particle
!
  EXTERNAL proj_area_particle
!
  REAL (RealK) :: volume_particle
!                   Function to calculate the volume of a particle
!
  EXTERNAL volume_particle
!
  INTERFACE

     SUBROUTINE size_integral_90                  &
         (nd_integral, npd_panel, npd_refinement, &
          npd_point, npd_panel_point,             &
          p_panel_ratio, SizeDist,                &
          number_total, proj_area_total,          &
          volume_total, dm_average,               &
          ierr, DBGeom)

          USE realtype_rd
          USE def_std_io_icf
          USE shape_particle_pcf
          USE def_size_dist
          USE prec_integral_tcf
          USE error_pcf
          USE rad_ccf, ONLY: pi
          USE def_db_crystal_geometry

          INTEGER, Intent(IN) :: nd_integral
          INTEGER, Intent(IN) :: npd_panel
          INTEGER, Intent(IN) :: npd_refinement
          INTEGER, Intent(IN) :: npd_point
          INTEGER, Intent(IN) :: npd_panel_point
          REAL  (RealK), Intent(IN) :: p_panel_ratio
          TYPE (STR_size_dist), Intent(IN) :: SizeDist
          REAL  (RealK), INTENT(OUT) :: number_total
          REAL  (RealK), INTENT(OUT) :: proj_area_total
          REAL  (RealK), INTENT(OUT) :: volume_total
          REAL  (RealK), INTENT(OUT) :: dm_average
          INTEGER, Intent(InOut) :: ierr      
          TYPE  (STR_db_cryst_geom), Intent(IN), Optional :: DBGeom

     END SUBROUTINE SIZE_INTEGRAL_90

  END INTERFACE


!
!
! Set the flag for interactive operation
  l_interactive=set_interactive()
!
! Calculate the sizes of the particles.
  CALL particle_size_90(l_interactive, nd_size, nd_mode, SizeDist, ierr)
  IF (ierr /= i_normal) STOP
!
! Select the scattering algorithm
  WRITE(iu_stdout, '(/A)') &
 &  'Enter the scattering algorithm.'
  DO
    READ(iu_stdin, *, IOSTAT=ios) i_scatter_algorithm
    IF ( (ios /= 0)                                   .OR. &
         ( (i_scatter_algorithm /= ip_algorithm_mie)  .AND. &
           (i_scatter_algorithm /= ip_adtmitchell96)  .AND. &
           (i_scatter_algorithm /= ip_scat_database) ) ) THEN
      WRITE(iu_err, '(A)') '+++ Erroneous input:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(A)') 'Please re-enter.'
      ELSE
        ierr = i_err_fatal
        STOP
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
!
! Define the appropriate number of moments of the phase function.
  WRITE(iu_stdout, '(/A)') &
 &  'Enter the number of moments of the phase function required.'
  DO
    READ(iu_stdin, *, IOSTAT=ios) n_phf_term
    IF (ios /= 0) THEN
      WRITE(iu_err, '(A)') '+++ Erroneous input:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(A)') 'Please re-enter.'
      ELSE
        ierr = i_err_fatal
        STOP
      ENDIF
    ELSE
      IF (n_phf_term > npd_phase_term) THEN
        WRITE(iu_err, '(A)') &
          '*** Error: Recompile with a larger ' // &
          'setting for npd_phase_term.'
        STOP
      ENDIF
      EXIT
    ENDIF
  ENDDO
!
! Open the scattering database.
  IF (i_scatter_algorithm == ip_scat_database) THEN
!
    CALL get_free_unit(ierr, iu_db_input)
    IF (ierr /= i_normal) STOP
    WRITE(iu_stdout, '(A)') &
      'Enter the filename of the binary scattering database.'
    READ(iu_stdin, '(A)') ice_bin_database
    OPEN(UNIT=iu_db_input, FILE=ice_bin_database, FORM='UNFORMATTED', &
      ACCESS='DIRECT', RECL=sct_db_recl)
    IF ( ierr /= i_normal ) STOP
!
  ENDIF
!

! Check the consistency of the algorithm and the shape of the particle.
  l_matched = .FALSE.
  IF (i_scatter_algorithm == ip_algorithm_mie) THEN
    l_matched = (SizeDist%i_shape_particle == ip_shape_sphere) 
  ELSE IF (i_scatter_algorithm == ip_adtmitchell96) THEN
    l_matched = (SizeDist%i_shape_particle == ip_shape_sphere)    .OR. &
                (SizeDist%i_shape_particle == ip_shape_column)    .OR. &
                (SizeDist%i_shape_particle == ip_shape_plate)     .OR. &
                (SizeDist%i_shape_particle == ip_shape_rosette)   .OR. &
                (SizeDist%i_shape_particle == ip_shape_polycrystal)
  ELSE IF (i_scatter_algorithm == ip_scat_database) THEN
    l_matched = (SizeDist%i_shape_particle == ip_shape_db_undefined)
  ENDIF
  IF (.NOT.l_matched) THEN
    WRITE(iu_err, '(/A, /A)') &
      '+++ Error: The scattering algorithm does not match ', &
      'the shape of the particle.'
    ierr = i_err_fatal
    STOP
  ENDIF
!
!
! The logic here is a little fuzzy: where should we have l_phase and
! l_stokes, and what are sensible conditions?
  l_phase  = (n_phf_term > 1)
  l_stokes = l_phase
!
! Set the file for output. 
  WRITE(iu_stdout, '(/A)') &
    'Give the name of the file to contain the scattering data.'
  DO
    READ(iu_stdin, '(A)', IOSTAT=ios) file_out
    IF (ios /= 0) THEN
      WRITE(iu_err, '(/A)') '*** Error: This name could not be read.'
      ierr = i_err_fatal
      STOP
    ENDIF
!   Check whether the file exists.
    INQUIRE(FILE=file_out, EXIST=l_exist)
    IF (l_exist) THEN
      WRITE(iu_err, '(/A)') '*** Error: This file already exists.'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(/A/)') 'Do you wish to overwrite? (Y/N).'
        READ(iu_stdin, '(A)') char_yn
        IF ( (char_yn /= 'Y') .AND. (char_yn /= 'y') ) THEN
          WRITE(iu_stdout, '(/A/)') &
            'Please specify another name or type !QUIT.'
          READ(iu_stdin, '(A)') file_out
          IF ( (file_out(1:5).EQ.'!QUIT').OR. &
               (file_out(1:5).EQ.'!quit').OR. &
               (file_out(1:5).EQ.'!Quit') ) THEN
            WRITE(iu_stdout, '(/A/)') '*** Program terminated.'
            ierr = i_err_io
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
!
! Open the output file now.
  CALL get_free_unit(ierr, iu_scatter)
  OPEN(UNIT=iu_scatter, FILE=file_out, IOSTAT=ios, STATUS='UNKNOWN')
  IF (ios /= 0) THEN
    WRITE(IU_ERR, '(/, 3A)') '*** Error: File ', FILE_OUT &
      , ' could not be opened.'
    ierr = i_err_io
    STOP
  ENDIF
!
! We now have enough information to determine the type of output file.
  IF (i_scatter_algorithm == ip_algorithm_mie) THEN
    i_output_type = it_file_mie_dry
  ELSE IF (i_scatter_algorithm == ip_adtmitchell96) THEN
    i_output_type = it_file_adt_dry
  ELSE IF (i_scatter_algorithm == ip_scat_database) THEN
    i_output_type = it_file_scat_database
  ENDIF
  i_scatter_type = ip_type_ice
!

  IF ( (SizeDist%i_distribution   == ip_external) .AND. &
       (SizeDist%i_shape_particle == ip_shape_db_undefined) ) THEN
    CALL ice_db_read_geometry(ierr, DBGeom)
  ENDIF


  IF ( (SizeDist%i_distribution == ip_external) .AND. &
       (SizeDist%n_size == 1) ) Then
    proj_area_total = proj_area_particle(SizeDist%dimen(1), &
                                         SizeDist%i_measure, &
                                         SizeDist%i_shape_particle, &
                                         ierr, DBGeom)
    volume_total  = volume_particle   (SizeDist%dimen(1), &
                                         SizeDist%i_measure, &
                                         SizeDist%i_shape_particle, &
                                         ierr, DBGeom)
    number_total    = SizeDist%number(1)
  ELSE
!   To do the geometrical calculations we require 4 integrals.
    nd_integral = 4
    CALL size_integral_90(nd_integral, &
      npd_panel, npd_refinement, npd_point, npd_panel_point, &
      p_panel_ratio, SizeDist, &
      number_total, proj_area_total, volume_total, dm_average, &
      ierr, DBGeom &
      )
  ENDIF
!
!
!
! Write the headers to the file
!
  WRITE(iu_scatter, '(/, A13, I5)') &
    '*FILE TYPE = ', I_OUTPUT_TYPE
  WRITE(iu_scatter, '(/, 16X, A, //)') &
    'SINGLE SCATTERING PARAMETERS.'
  WRITE(iu_scatter, '(A23, I5, 3X, A1, A40, A1)') &
    'SCATTERING ALGORITHM = ', i_scatter_algorithm &
    , '(', name_scatter_algorithm(i_scatter_algorithm), ')'
  WRITE(iu_scatter, '(A23, I5, 3X, A1, A20, A1)') &
    'SCATTERING TYPE      = ', i_scatter_type, &
    '(', name_scatter_type(i_scatter_type), ')'
  WRITE(iu_scatter, '(/, A18, 1X, 1PE12.5, 1X, A3)') &
    'NUMBER DENSITY   =', number_total, 'm-3'
  WRITE(iu_scatter, '(A18, 1X, 1PE12.5, 1X, A6)') &
    'PROJECTED AREA   =', proj_area_total, 'm2.m-3'
  WRITE(iu_scatter, '(A18, 1X, 1PE12.5, 1X, A6)') &
    'VOLUME FRACTION  =', volume_total, 'm3.m-3'
  WRITE(iu_scatter, '(A18, 1X, 1PE12.5, 1X, A1)') &
    'DE               =', 1.5_RealK*volume_total/proj_area_total, 'm'
  WRITE(iu_scatter, '(A18, 1X, 1PE12.5, 1X, A1)') &
    'DM (MEANMAXDIM)  =', dm_average, 'm'
  WRITE(iu_scatter, '(/, A41, I3, /)') &
    'NUMBER OF TERMS IN THE PHASE FUNCTION =  ', N_PHF_TERM
  WRITE(iu_scatter, '(/, 4X, A14, 6X, A16, 4X, A16, 4X, A14)') &
    'Wavelength (m)', 'Absorption (m-1)' &
    , 'Scattering (m-1)', 'Moments of PF'
!
!
  IF (i_scatter_algorithm == ip_scat_database) THEN
!
!   Read the wavelengths from the database
    CALL get_db_wavelengths(iu_db_input, &
      npd_wavelength_scat, npd_size_scat, &
      n_wavelength, wavelength, &
      n_rec_block, db_record, db_n_angle, ierr) 
!
!   Setup Gaussian roots and weights
!   (The code seems not to work if n_angle is set to too high a value:
!   this requires investigation.)
!   An odd order of quadrature is required by the integration.
    n_angle = n_phf_term + 1 + MOD(n_phf_term, 2)
    IF (n_angle > npd_scatt_angle) THEN
      WRITE(iu_err, '(/A)') &
        '*** Error: Too many angles are required to evaluate the ' // &
        'phase function.'
      STOP
    ENDIF
    CALL calc_gauss_weight_90(ierr, n_angle, mu_angle, weight_angle)
!
!
  ELSE
!
!   Read in the file of wavelengths 
!   where scattering is to be calculated.
    CALL get_wavelengths(npd_wavelength_scat, n_wavelength, &
      wavelength, ierr)
!
!   Obtain the relevant refractive indices.
    CALL get_refract_index(l_interactive, npd_refract, n_refract, &
      wavelength_refract, re_refract, im_refract, &
      l_spline_interp, d2_re_refract, d2_im_refract, &
      ierr)
!
  ENDIF
!
!
!
! Integration is done over the range of scattering angles required.
  nd_integral=n_angle + 3
!
! Loop over wavelengths:
!
  DO i = 1, n_wavelength
!
    IF (i_scatter_algorithm == ip_scat_database) THEN
!
!     Read all information valid at this wavelength from the database.
!     Set arrays for interpolation.
      CALL db_read_single_wavelength(iu_db_input, &
        npd_wavelength_scat, npd_scatt_angle, npd_size_scat, &
        i, n_wavelength, wavelength, &
        n_rec_block, db_record, &
        l_phase, n_angle, mu_angle, &
        ice_db_mono_info, &
        ierr)
!
      CALL db_scatter_integral(npd_wavelength_scat, npd_size_scat, &
                               SizeDist,i,n_wavelength,wavelength, &
                               DBGeom, ice_db_mono_info, &
                               n_angle,mu_angle, &
                               p_panel_ratio, &
                               extinction,scattering,asymmetry, &
                               l_stokes,i_stokes, &
                               npd_scatt_angle, &
                               nd_integral,npd_panel,npd_refinement, &
                               npd_point,npd_panel_point, &
                               ierr)
!
!     Deallocate the information on ice crystals.
      DEALLOCATE(ice_db_mono_info%dm)
      DEALLOCATE(ice_db_mono_info%ss)
      DEALLOCATE(ice_db_mono_info%d2_ss)
!
    ELSE 
!
      n_relative = refractive_index(ierr, npd_refract, wavelength(i), &
        n_refract, wavelength_refract, re_refract, im_refract, &
        d2_re_refract, d2_im_refract)
!
      IF ( (SizeDist%i_distribution == ip_external) .AND. &
           (SizeDist%n_size == 1) ) THEN
!
        IF (SizeDist%i_measure == ip_measure_radius) THEN
          size = 2.0_RealK * pi * SizeDist%dimen_0 / wavelength(i)
        ELSE IF (SizeDist%i_measure == ip_measure_max_dimen) THEN
          size = pi * SizeDist%dimen_0 / wavelength(i)
        ENDIF
!
        IF (i_scatter_algorithm == ip_algorithm_mie ) THEN
          CALL mie_scatter(ierr, size, n_relative, &
            q_scat, q_ext, asymmetry, &
            l_phase, n_angle, mu_angle, s1, s2, &
            npd_scatt_angle, npd_log_deriv)
          IF (ierr /= i_normal) THEN
            ierr = i_err_fatal
            STOP
          ENDIF
          DO j = 1, n_angle
            i_stokes(j) = SizeDist%number(1) * &
              (ABS(s1(j))**2 + ABS(s2(j))**2)
          ENDDO
        ELSE IF (i_scatter_algorithm == ip_adtmitchell96 ) THEN
          CALL mie_scatter(ierr, size, n_relative, &
            q_scat, q_ext, asymmetry, &
            l_phase, n_angle, mu_angle, s1, s2, &
            npd_scatt_angle, npd_log_deriv)
          IF (ierr /= i_normal) THEN
            ierr = i_err_fatal
            STOP
          ENDIF
          DO j = 1, n_angle
            i_stokes(j) = SizeDist%number(1) * &
              (ABS(s1(j))**2 + ABS(s2(j))**2)
          ENDDO
        ENDIF
!
        extinction = q_ext  * proj_area_total
        scattering = q_scat * proj_area_total
!
      ELSE
!
        IF ( (i_scatter_algorithm == ip_algorithm_mie ) .OR. &
             (i_scatter_algorithm == ip_adtmitchell96 ) ) THEN
!
          CALL scatter_integral_90(SizeDist, wavelength(i), n_relative, &
            i_scatter_algorithm, &
            n_angle, mu_angle, &
            p_panel_ratio, &
            extinction, scattering, asymmetry, &
            l_stokes, i_stokes, &
            npd_scatt_angle, &
            npd_log_deriv, nd_integral, npd_panel, npd_refinement, &
            npd_point, npd_panel_point, &
            ierr &
            )
!
        ENDIF
!
      ENDIF
!
    ENDIF
!
!
    absorption = extinction - scattering
!
    phase_fnc_term(1) = asymmetry
!

!   Decompose the phase function into its moments.
    IF (l_phase) THEN
      CALL decompose_phf_90(n_angle, mu_angle, weight_angle, &
                         n_phf_term, i_stokes, phase_fnc_term, &
                         ierr &
                         )
!     Renormalize the moments to the asymmetry. This is done mainly
!     for consistency with the external database. Direct integration
!     of the phase function tends to give low values (because of the
!     highly peaked nature of the phase function), while the actual
!     valuse are available by direct integration. This renormalization
!     is entirely heuristic.
      DO j=1, n_phf_term
        phase_fnc_term(j) = phase_fnc_term(j) * &
          asymmetry / phase_fnc_term(1)
      ENDDO
!
    ENDIF
!
!   Add the results to the output file.
    WRITE(iu_scatter, '(4(4X, 1PE16.9), (20X, 3(4X, 1PE16.9)))') &
      wavelength(i), absorption, scattering, &
      (phase_fnc_term(j), j=1, n_phf_term)
!
  ENDDO
!
  WRITE(iu_scatter, '(1x)')
!
  CLOSE(iu_scatter)
!
!
END

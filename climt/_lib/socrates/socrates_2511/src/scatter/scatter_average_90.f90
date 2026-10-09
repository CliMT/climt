! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to average scattering properties across bands.
!
PROGRAM scatter_average_90
!
! Description:
!    This program reads a file of monochromatic single-scattering
! data and averages them across the bands specified in a spectral
! file. Optionally, the data may be fitted.
!
! Method:
!    A file containing blocks of monochromatic single 
! scattering properties is read in. These monochromatic
! values are averaged across the bands given in a 
! spectral file. The averaged values may be written to
! a file or fitted using a recognized parametrization.
!
!
!
! Modules used
  USE realtype_rd
  USE def_spectrum
  USE def_inst_flt
  USE dimensions_pp_ucf
  USE dimensions_spec_ucf
  USE scatter_pp_pcf
  USE def_solarspec
  USE def_s_scat_prop
  USE method_weight_pcf
  USE def_std_io_icf
  USE file_type_pcf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
! Declaration of variables.
  CHARACTER (LEN=256) :: file_spectral
!   Name of spectral file
  TYPE  (StrSpecData) :: Spectrum
!   Spectral data used for weighting
!
! I/O varaibles
  INTEGER :: iu_scatter
!   Unit number for reading scattering data
  INTEGER :: iu_average
!   Unit number for writing out averaged data
  INTEGER :: ios
!   I/O error flag
  LOGICAL :: l_interactive
!   Flag for interactive use
!
  CHARACTER (LEN=1) :: char_yn
!   Character response
  INTEGER :: ierr = i_normal
!       Error flag
  INTEGER :: n_int_weight
!       Number of weighting points
  INTEGER :: i_weight
!       Type of weighting function
  INTEGER :: i_method_weight
!       Method of weighting
  INTEGER :: n_block
!       Number of data block
  INTEGER :: n_block_p1
!       Number of data block plus 1
  INTEGER :: i_begin
!       Initial point in band
  INTEGER :: i_band
!       Index of band
  INTEGER :: i
!       Loop variable
  INTEGER :: j
!       Loop variable
  INTEGER :: k
!       Loop variable
  INTEGER :: ipb
!       Offset loop variable
  INTEGER :: i_input_type
!       Type of input file
  INTEGER :: i_output_type
!       Type of output file
  LOGICAL :: l_end = .FALSE.
!       End of file flag
  LOGICAL :: l_write
!       Writing flag
  LOGICAL :: l_wave_number =.FALSE.
!       Flag for weighting by wavenumber
  REAL  (RealK) :: t_weight
!       Temperature for weighting
  REAL  (RealK) :: weight_coeff(0: npd_wavelength_scat)
!       Coefficients for weighting
  REAL  (RealK) :: work(0: npd_wavelength_scat)
!       Work array
!
  TYPE (StrSingleScatSzD) :: ss_data
!   Single scattering data averaged over a size distribution
!
  INTEGER :: i_scatter_type(npd_mie_block)
!       Type of scatterer in each block
  REAL  (RealK) :: volume_fraction(npd_mie_block)
!       Volume fraction in each block of spectral data
  REAL  (RealK) :: dim_char(npd_mie_block)
!       Characteristic dimensions in each block of spectral data
  REAL  (RealK) :: particle_density(npd_mie_block)
!       Particle density in each block of spectral data
  CHARACTER  (LEN=5) :: fit_species
!   Phase of condensate to be fitted
!
  TYPE  (StrSolarSpec) :: SolarSpec
!       Spectral solar irradiance
!
  LOGICAL :: include_instrument_response
!   Flag to include the instrumental response function
  TYPE  (StrFiltResp) :: filter
!   Instrumental response function
!
  INTEGER :: n_phf_term_fit
!       Number of terms in the phase function provided
  REAL  (RealK) :: absorption_band(npd_band, npd_mie_block)
!       Band-averaged absorption
  REAL  (RealK) :: scattering_band(npd_band, npd_mie_block)
!       Band-averaged absorption
  REAL  (RealK) :: phase_fnc_band(npd_band, npd_phase_term, npd_mie_block)
!       Band-averaged absorption
!
  REAL  (RealK) :: exclude_low(npd_exclude)
!       Lower limits of excluded intervals (unit: m-1)
  REAL  (RealK) :: exclude_high(npd_exclude)
!       Upper limits of excluded intervals (unit: m-1)
  REAL  (RealK) :: nu_min, nu_max
!       Limits of wavenumbers required by spectral file
!
!
! External functions:
  LOGICAL, EXTERNAL :: set_interactive
!   Function to set the flag for interactive operation


  INTERFACE
!
    SUBROUTINE cloud_fit_90 &
!
      (l_interactive, fit_species, n_band, n_data, vol_frac, d, &
       particle_density, &
       absorption_ave, scattering_ave, n_phf_term, phf_fnc_ave, &
       ierr)
!
      USE realtype_rd
!
!
!
      LOGICAL, Intent(IN) :: l_interactive
      CHARACTER  (LEN=5), Intent(IN) :: fit_species
      INTEGER, Intent(IN) :: n_band
      INTEGER, Intent(IN) :: n_data
!
      INTEGER, Intent(INOUT) :: ierr
      REAL  (RealK), Intent(IN), Dimension(:) :: vol_frac
      REAL  (RealK), Intent(IN), Dimension(:) :: d
      REAL  (RealK), Intent(IN), Dimension(:) :: particle_density
!
      REAL  (RealK), Intent(IN), Dimension(:, :) :: absorption_ave
      REAL  (RealK), Intent(IN), Dimension(:, :) :: scattering_ave
      INTEGER, Intent(IN) :: n_phf_term
      REAL  (RealK), Intent(IN), Dimension(:, :, :) :: phf_fnc_ave
!
!
!
    END SUBROUTINE cloud_fit_90
!
!
!
  END INTERFACE
!
!
!
!
! Set the flag for interactive operation
  l_interactive=set_interactive()
!
! Read in the spectral file.
  WRITE(*, "(a)") "Enter the name of the spectral file."
  READ(*, "(a)") file_spectral
  CALL read_spectrum(file_spectral, Spectrum)

! Find limits of wavenumbers required by spectral file
  nu_min=MINVAL(1.0_RealK/Spectrum%Basic%wavelength_long( &
                        1:Spectrum%Basic%n_band )) - EPSILON(nu_min)
  nu_max=MAXVAL(1.0_RealK/Spectrum%Basic%wavelength_short( &
                        1:Spectrum%Basic%n_band )) + EPSILON(nu_max)

! Open the scattering files.
  CALL get_free_unit(ierr, iu_scatter)
  IF (ierr /= i_normal) STOP
  CALL open_file_in(ierr, iu_scatter, &
    'Give the name of the file of scattering data to be averaged.')
  IF (ierr /= i_normal) STOP
!
! Select the weighting function.
  CALL select_weight_scatter_90(i_weight, SolarSpec, t_weight, &
    i_method_weight, l_interactive, ierr)
  IF (ierr /= i_normal) stop
!
! Read in instrument response if required
  CALL get_inst_response_int
!
! Determine whether the averages across each band should be written
! to a file.
  CALL get_free_unit(ierr, iu_average)
  IF (ierr /= i_normal) STOP
  CALL open_average_90(l_write, iu_average, l_interactive, ierr)
  IF (ierr /= i_normal) stop
!
! How many moments of the phase function for the fit?
  WRITE(iu_stdout, '(a)') &
    'How many moments of the phase function are to be processed?'
  DO
    READ(iu_stdin, *, IOSTAT=ios) n_phf_term_fit
    IF (ios /= 0) THEN
      WRITE(iu_err, "(a)") "+++ Erroneous response."
      IF (l_interactive) THEN
        WRITE(*, "(a)") "Please re-specify"
      ELSE
        STOP
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
! Each block of data in the scattering file is read and averaged
! separately. We allocate a large amount of space and then read,
! rather then trying to do the allocation dynamically for now.

  ALLOCATE(ss_data%wavenum(npd_wavelength_scat))
  ALLOCATE(ss_data%abs(npd_wavelength_scat))
  ALLOCATE(ss_data%scat(npd_wavelength_scat))
  ALLOCATE(ss_data%phf(npd_phase_term, npd_wavelength_scat))

  n_block=0
  DO
!   Check for the end of the input file.
    READ(iu_scatter, '(A)', IOSTAT=ios) char_yn
    IF (ios /= 0) EXIT
    BACKSPACE(iu_scatter)
!
    n_block = n_block + 1
    CALL read_scatter_block_90(iu_scatter, i_input_type, &
      ss_data, nu_min, nu_max, &
      npd_wavelength_scat, npd_phase_term, ierr)
    IF (ierr /= i_normal) STOP
!
!
!   Average across each band.
    DO i_band=1, Spectrum%Basic%n_band
!
      CALL set_scatter_weight_int
!
!     Calculate the mean scattering properties. For later fitting we keep
!     the volume fraction and characteristic dimension of each block.
      CALL calculate_means_int
!
    ENDDO
!
    i_scatter_type(n_block)=ss_data%i_scatter_type
    volume_fraction(n_block)=ss_data%vol_frac
    dim_char(n_block)=ss_data%dim_char
    IF (i_input_type == it_file_scat_mass) THEN
      particle_density(n_block)=ss_data%particle_density
    ELSE
      particle_density(n_block)=0.0_RealK
    END IF
!
    IF (l_write) THEN
!     Set the output type for the file of results.
      IF (i_input_type == it_file_mie_dry) THEN
        i_output_type=it_file_ave_mie_dry
      ELSE IF (i_input_type == it_file_phf_mie_dry) THEN
        i_output_type=it_file_ave_phf_mie_dry
      ELSE IF (i_input_type == it_file_mie_humid) THEN
        i_output_type=it_file_ave_mie_humid
      ELSE IF (i_input_type == it_file_phf_mie_humid) THEN
        i_output_type=it_file_ave_phf_mie_humid
      ELSE IF (i_input_type == it_file_scat_database) THEN
        i_output_type=it_file_ave_phf_mie_dry
      ELSE IF (i_input_type == it_file_scat_mass) THEN
        i_output_type=it_file_ave_scat_mass
      ENDIF
      CALL write_average_90(iu_average, Spectrum%Basic%n_band, &
        n_block, ss_data, i_output_type, &
        absorption_band(1, n_block), &
        scattering_band(1, n_block), &
        n_phf_term_fit, phase_fnc_band(1, 1, n_block), &
        npd_band, npd_phase_term &
        )
    ENDIF
  ENDDO
!
  DEALLOCATE(ss_data%wavenum)
  DEALLOCATE(ss_data%abs)
  DEALLOCATE(ss_data%scat)
  DEALLOCATE(ss_data%phf)
!
!
  CLOSE(iu_average)
!
!
! Fit the data if required.
  WRITE(iu_stdout, '(/a)') &
    'Do you want to parametrize these data ? (y/n)'
  DO
!
    READ(iu_stdin, '(a)') char_yn
!
    IF ( (char_yn == 'y') .OR. (char_yn == 'Y') ) THEN
!
!     The data are checked to ensure that they are all for 
!     the same type of scatterer.
      DO i=2, n_block
        IF (i_scatter_type(i) /= i_scatter_type(1)) THEN
          WRITE(iu_err, '(/a, /a)') &
            '*** Error: these blocks do not all represent ' // &
            'the same type of scatterer.', &
            'Fitting will not be performed.'
          ierr=i_err_fatal
          STOP
        ENDIF
      ENDDO
!
      IF ( (i_scatter_type(n_block) == IP_type_droplet) .OR. &
           (i_scatter_type(n_block) == IP_type_ice) ) THEN
!
        IF (i_scatter_type(n_block) == IP_type_ice) THEN
          fit_species = "Ice  "
        ELSE IF (i_scatter_type(n_block) == IP_type_droplet) THEN
          fit_species = "Water"
        ENDIF
!
        CALL cloud_fit_90(l_interactive, fit_species, &
          Spectrum%Basic%n_band, n_block, &
          volume_fraction, dim_char, particle_density, &
          absorption_band, scattering_band, &
          n_phf_term_fit, phase_fnc_band, &
          ierr)
      ELSE
        WRITE(iu_err, "(a)") &
          "*** Error: Fits for this type of scatterer are not provided."
      ENDIF
!
      EXIT
!
    ELSE IF ( (char_yn == 'n') .OR. (char_yn == 'N') ) THEN
      STOP
    ELSE 
!
      WRITE(iu_err, '(a)') '+++ Illegal response.'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        STOP
      ENDIF
!
    ENDIF
!
  ENDDO


CONTAINS
!
!
!
  SUBROUTINE get_inst_response_int
!
!
!
    CHARACTER  (LEN=1) :: char_if
!     I/O response variable
!
!
!
    WRITE(iu_stdout, "(/A)") &
      "Is an instrumental response required? (Y/N)"
    Inst: DO
!
      READ(iu_stdin, '(A)') char_if
!
      IF ( (char_if == 'Y') .OR. (char_if == 'y') ) THEN
!
        include_instrument_response=.TRUE.
        CALL read_instrument_response_90(filter, ierr)
        IF (ierr /= i_normal) STOP
        EXIT
!
      ELSE IF ( (char_if == 'N') .OR. (char_if == 'n') ) THEN
!
        include_instrument_response=.FALSE.
        EXIT
!
      ELSE
!
        WRITE(iu_err, '(/A)') '*** Error: Unrecognized input'
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(/A)') 'Please re-enter.'
        ELSE
          STOP
        ENDIF
!
      ENDIF
!
    ENDDO Inst  
!
!
!
  END SUBROUTINE get_inst_response_int
!
!
!
  SUBROUTINE set_scatter_weight_int
!
!
!
!   Calculate the weightings.
    DO k=1, Spectrum%Basic%n_band_exclude(i_band)
      exclude_high(k) = 1.0 / &
        Spectrum%Basic%wavelength_short( &
          Spectrum%Basic%index_exclude(k, i_band))
      exclude_low(k) = 1.0 / &
        Spectrum%Basic%wavelength_long( &
          Spectrum%Basic%index_exclude(k, i_band))
    ENDDO
    CALL weightings_90(ierr, &
      i_weight, &
      SolarSpec, &
      t_weight, &
      include_instrument_response, filter, &
      1.0_RealK/Spectrum%Basic%wavelength_long(i_band), &
      1.0_RealK/Spectrum%Basic%wavelength_short(i_band), &
      Spectrum%Basic%n_band_exclude(i_band), &
      exclude_low, exclude_high, &
      ss_data%n_wavenumber, ss_data%wavenum, &
      i_begin, n_int_weight, weight_coeff &
      )
    IF (ierr /= i_normal) STOP
!
!
!
  END SUBROUTINE set_scatter_weight_int
!
!
!
  SUBROUTINE calculate_means_int
!

!   Variables for thick weighting:
    REAL  (RealK) :: phi
!     Square root of (1-omega)
    REAL  (RealK) :: psi
!     Square root of (1-g*omega)
    REAL  (RealK) :: reflect_inf
!     Monochromatic reflection coefficient for an infinite cloud
    REAL  (RealK) :: extinction_tot
!     Total monochromatic extinction coefficient
    REAL  (RealK) :: mean_reflect_inf
!     Mean infinite reflection coefficient
    REAL  (RealK) :: mean_extinction
!     Mean total extinction
    REAL  (RealK) :: mean_scattering
!     Mean scattering
    REAL  (RealK) :: mean_scattering_phase_fnc(npd_phase_term)
!     Mean of product of scattering and phase function
    REAL  (RealK) :: omega_eff
!     Effective albedo of single scattering
!
!
!    
!   Calculate the mean scattering.
    mean_scattering = 0.0_RealK
    DO k = 0, n_int_weight
      mean_scattering = &
        mean_scattering + weight_coeff(k) * ss_data%scat(i_begin + k)
    ENDDO
    DO j = 1, n_phf_term_fit
      mean_scattering_phase_fnc(j) = 0.0_RealK
      DO k = 0, n_int_weight
        mean_scattering_phase_fnc(j) = &
          mean_scattering_phase_fnc(j) + &
          weight_coeff(k) * ss_data%scat(i_begin + k) * &
          ss_data%phf(j, i_begin + k)
      ENDDO
      IF (mean_scattering > 0.0_RealK) THEN
        phase_fnc_band(i_band, j, n_block) = &
          mean_scattering_phase_fnc(j) / mean_scattering
      ELSE
        phase_fnc_band(i_band, j, n_block) = 0.0_RealK
      END IF
    ENDDO
!
!   The absorption may be weighted in different ways for
!   different purposes.
!
    IF (i_method_weight == ip_weight_thin) THEN
!
      scattering_band(i_band, n_block) = &
        mean_scattering
      absorption_band(i_band, n_block) = 0.0_RealK
      DO k = 0, n_int_weight
        absorption_band(i_band, n_block) = &
          absorption_band(i_band, n_block) + &
          weight_coeff(k) * ss_data%abs(i_begin + k)
      ENDDO
!
    ELSE IF (i_method_weight == ip_weight_thick) THEN
!
!     The reflection coefficient of an infinite cloud is
!     calculated and weighted, and the total extinction is
!     weighted as opposed to the individual coefficients.
      mean_extinction = 0.0_RealK
      mean_reflect_inf = 0.0_RealK
      DO k = 0, n_int_weight
        ipb = k + i_begin
        extinction_tot = ss_data%abs(ipb) + ss_data%scat(ipb)
        phi = SQRT(ss_data%abs(ipb) / extinction_tot)
        psi = SQRT((ss_data%abs(ipb) + ss_data%scat(ipb) * &
                   (1.0_RealK - ss_data%phf(1, ipb))) / &
                   extinction_tot)
        reflect_inf = (psi - phi)/(psi + phi)
        mean_extinction = mean_extinction + &
          weight_coeff(k) * extinction_tot
        mean_reflect_inf = mean_reflect_inf + &
          weight_coeff(k) * reflect_inf
      ENDDO
!     Determine the mean value of omega.
      omega_eff = 4.0_RealK * mean_reflect_inf / &
        ((1.0_RealK + mean_reflect_inf)**2 - &
        phase_fnc_band(i_band, 1, n_block) * &
        (1.0_RealK - mean_reflect_inf)**2)
      absorption_band(i_band, n_block) = (1.0_RealK - omega_eff) * &
        mean_extinction
      scattering_band(i_band, n_block) = omega_eff * &
        mean_extinction
!
    ELSE IF (i_method_weight == ip_weight_thin_thick) THEN
!
      WRITE(iu_err, '(/a)') &
        '*** Error: the implementation of this scheme ' // &
        'is not yet complete.'
      ierr = i_err_fatal
!
    ENDIF
!
!
!
  END SUBROUTINE calculate_means_int
!
!
!
END PROGRAM scatter_average_90

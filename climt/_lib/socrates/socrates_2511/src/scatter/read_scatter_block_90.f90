! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read one block of scattering data.
!
SUBROUTINE read_scatter_block_90(iu_scatter, i_input_type, ss_data, &
  nu_min, nu_max, nd_wavelength_scat, nd_phase_term, ierr)
!
!- ---------------------------------------------------------------------
!
!
!
! Modules to set types of variables:
  USE realtype_rd
  USE def_std_io_icf
  USE def_s_scat_prop
  USE file_type_pcf
  USE scatter_pp_pcf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  INTEGER, Intent(IN) :: iu_scatter
!   Uint number for reading single scattering data
!
  TYPE (StrSingleScatSzD), Intent(OUT) :: ss_data
!   Single scattering data averaged over a size distribution

  REAL (RealK), INTENT(IN) :: nu_min, nu_max
!   Limits of wavenumbers required by spectral file

  INTEGER, Intent(OUT) :: i_input_type
!   Type of input file
!
  INTEGER, Intent(IN) :: nd_wavelength_scat
!   Size allocated for wavelengths
  INTEGER, Intent(IN) :: nd_phase_term
!   Size allocated for terms in the phase function
!        
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
!
! Local variables.
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i_algorithm
!   Type of scattering algorithm
  INTEGER :: j
!   Loop variable
  INTEGER :: k
!   Loop variable
  CHARACTER  (LEN=80) :: line
!   Input line of data
  LOGICAL :: non_blank
!   Test for blank line
  INTEGER , Allocatable :: map(:)
!   Mapping array for sorting by wavenumber
  REAL  (RealK), Allocatable :: dummy(:)
!   Temporary array used to reorder single scattering properties
  INTEGER :: loc_min, loc_max
!   Location of min/max wavenumber values in data
!
!
! Functions called:
!
  INTERFACE
!
    SUBROUTINE map_heap_func(a, map)
!
      USE realtype_rd
!
      REAL  (RealK), Intent(IN), Dimension(:) :: a
!
      INTEGER, Intent(OUT), Dimension(:) :: map
!
    END SUBROUTINE map_heap_func
!
  END INTERFACE
!
!
!
!
! Read the input type. This is used to control the
! subsequent reading in.
!
  READ(iu_scatter, '(/, 13x, i5)', IOSTAT=ios) i_input_type
  IF (ios /= 0) THEN
    WRITE(iu_err, '(/A)') &
      '*** Error: Unable to determine type of input file.'
  ENDIF
!
! In this developmental phase we set the code to the input type:
! this may not be the best practice in the long run.
  ss_data%code = i_input_type
!
  SELECT CASE(i_input_type)
!
    CASE (IT_file_mie_dry) 
!     Read header without humidity variable.
      READ(iu_scatter, '(////, 23x, i5, /, 23x, i5, 4x, a20)') &
        i_algorithm, ss_data%i_scatter_type, ss_data%name_type
      IF (ss_data%i_scatter_type == IP_type_aerosol) THEN
        READ(iu_scatter, '(23x, i5, 4x, a20)') &
          ss_data%i_component, ss_data%name_component
      ENDIF
      READ(iu_scatter, '(/, 3(19x, 1pe12.5, /), /)') &
        ss_data%number_density, ss_data%dim_char, ss_data%vol_frac
      ss_data%n_phf_term=1
!
    CASE(IT_file_phf_mie_dry)
!     Read header without humidity variable.
      READ(iu_scatter, '(////, 23x, i5, /, 23x, i5, 4x, a20)') &
        i_algorithm, ss_data%i_scatter_type, ss_data%name_type
      IF (ss_data%i_scatter_type == IP_type_aerosol) THEN
        READ(iu_scatter, '(23x, i5, 4x, a20)') &
          ss_data%i_component, ss_data%name_component
      ENDIF
      READ(iu_scatter, '(/, 3(19x, 1pe12.5, /))') &
        ss_data%number_density, ss_data%dim_char, ss_data%vol_frac
      READ(iu_scatter, '(41x, i3, ///)') ss_data%n_phf_term
!
    CASE (IT_file_mie_humid)
!     Read header with humidity variable.
      READ(iu_scatter, '(////, 23x, i5, /, 23x, i5, 4x, a20)') &
        i_algorithm, ss_data%i_scatter_type, ss_data%name_type
      IF (ss_data%i_scatter_type == IP_type_aerosol) THEN
        READ(iu_scatter, '(23x, i5, 4x, a20)') &
          ss_data%i_component, ss_data%name_component
      ENDIF
      READ(iu_scatter, '(/, 11x, f7.5, //, 19x, 1pe12.5, /, ' // &
        ' 2(19x, 1pe12.5, /, 36x, 1pe12.5, /), /)') &
          ss_data%humidity, ss_data%number_density, &
          ss_data%dim_char, ss_data%re_dry, &
          ss_data%vol_frac, ss_data%vol_frac_dry
      ss_data%n_phf_term=1
!
    CASE(IT_file_phf_mie_humid)
!     Read header with humidity variable.
      READ(iu_scatter, '(////, 23x, i5, /, 23x, i5, 4x, a20)') &
        i_algorithm, ss_data%i_scatter_type, ss_data%name_type
      IF (ss_data%i_scatter_type == IP_type_aerosol) THEN
        READ(iu_scatter, '(23x, i5, 4x, a20)') &
          ss_data%i_component, ss_data%name_component
      ENDIF
      READ(iu_scatter, '(/, 11x, f7.5, //, 19x, 1pe12.5, /,  ' // &
        ' 2(19x, 1pe12.5, /, 36x, 1pe12.5, /))') &
          ss_data%humidity, ss_data%number_density, &
          ss_data%dim_char, ss_data%re_dry, &
          ss_data%vol_frac, ss_data%vol_frac_dry
      READ(iu_scatter, '(41x, i3, ///)') ss_data%n_phf_term
!
    CASE(IT_file_scat_database)
!     Read headers from a file constructed from a database.
      READ(iu_scatter, '(////, 23x, i5, /, 23x, i5, 4x, a20)') &
        i_algorithm, ss_data%i_scatter_type, ss_data%name_type
      IF (ss_data%i_scatter_type == IP_type_aerosol) THEN
        READ(iu_scatter, '(23x, i5, 4x, a20)') &
          ss_data%i_component, ss_data%name_component
      ENDIF
      READ(iu_scatter, '(/, 19x, 1pe12.5, a3)') ss_data%number_density
      READ(iu_scatter, '(/, 19x, 1pe12.5, a3)') ss_data%vol_frac
      READ(iu_scatter, '(19x, 1pe12.5, a3)') ss_data%dim_char
      READ(iu_scatter, '(//, 41x, i3, ///)') ss_data%n_phf_term

    CASE(IT_file_scat_mass)
      READ(iu_scatter, '(////, 23x, i5, /, 23x, i5, 4x, a20)') &
        i_algorithm, ss_data%i_scatter_type, ss_data%name_type
      READ(iu_scatter, '(/, 30x, 1pe12.5)') ss_data%mass_mixing_ratio
      READ(iu_scatter, '(/, 30x, 1pe12.5)') ss_data%air_density
      READ(iu_scatter, '(30x, 1pe12.5)') ss_data%particle_density
      ss_data%vol_frac = ss_data%mass_mixing_ratio * ss_data%air_density &
                       / ss_data%particle_density
      READ(iu_scatter, '(30x, 1pe12.5)') ss_data%dim_char
      READ(iu_scatter, '(40x, i3, //)') ss_data%n_phf_term

    CASE DEFAULT
      WRITE(iu_err, '(/a)')  &
        '*** Error: the input file is of an invalid type.'
      ierr=i_err_fatal
      RETURN
!
  END SELECT
!
! Read and process each line in turn.
  ss_data%n_wavenumber=0
  DO
    READ(iu_scatter, '(a)') line
    IF (non_blank(line)) THEN
      IF (ss_data%n_wavenumber < nd_wavelength_scat) THEN
        ss_data%n_wavenumber=ss_data%n_wavenumber+1
      ELSE
        WRITE(iu_err, '(/a)') &
          '*** Error: there are too many wavelengths in the ' // &
          'input file.'
        WRITE(iu_err, '(a)') &
          'Increase npd_wavelength_scat and recompile.'
        ierr=i_err_fatal
        RETURN
      ENDIF
      BACKSPACE(iu_scatter)
      READ(iu_scatter, '(4(4x, 1pe16.9), :, /, (20x, 3(4x, 1pe16.9)))') &
        ss_data%wavenum(ss_data%n_wavenumber), &
        ss_data%abs(ss_data%n_wavenumber), &
        ss_data%scat(ss_data%n_wavenumber), &
        (ss_data%phf(j, ss_data%n_wavenumber), j=1, ss_data%n_phf_term)
!     Convert to wavenumber.
      ss_data%wavenum(ss_data%n_wavenumber) = 1.0_RealK / &
        ss_data%wavenum(ss_data%n_wavenumber)
    ELSE
      EXIT
    ENDIF
  ENDDO

! Extend boundary values to the limits required by the spectral file
  loc_min=MINLOC(ss_data%wavenum(1:ss_data%n_wavenumber),1)
  IF (ss_data%wavenum(loc_min) > nu_min) THEN
    IF (ss_data%n_wavenumber < nd_wavelength_scat) THEN
      ss_data%n_wavenumber=ss_data%n_wavenumber+1
      ss_data%wavenum(ss_data%n_wavenumber) = nu_min
      ss_data%abs(ss_data%n_wavenumber)=ss_data%abs(loc_min)
      ss_data%scat(ss_data%n_wavenumber)=ss_data%scat(loc_min)
      ss_data%phf(:,ss_data%n_wavenumber)=ss_data%phf(:,loc_min)
    ELSE
      WRITE(iu_err, '(/a)') &
        '*** Error: there are too many wavelengths in the ' // &
        'input file.'
      WRITE(iu_err, '(a)') &
        'Increase npd_wavelength_scat and recompile.'
      ierr=i_err_fatal
      RETURN
    END IF
  END IF
  loc_max=MAXLOC(ss_data%wavenum(1:ss_data%n_wavenumber),1)
  IF (ss_data%wavenum(loc_max) < nu_max) THEN
    IF (ss_data%n_wavenumber < nd_wavelength_scat) THEN
      ss_data%n_wavenumber=ss_data%n_wavenumber+1
      ss_data%wavenum(ss_data%n_wavenumber) = nu_max
      ss_data%abs(ss_data%n_wavenumber)=ss_data%abs(loc_max)
      ss_data%scat(ss_data%n_wavenumber)=ss_data%scat(loc_max)
      ss_data%phf(:,ss_data%n_wavenumber)=ss_data%phf(:,loc_max)
    ELSE
      WRITE(iu_err, '(/a)') &
        '*** Error: there are too many wavelengths in the ' // &
        'input file.'
      WRITE(iu_err, '(a)') &
        'Increase npd_wavelength_scat and recompile.'
      ierr=i_err_fatal
      RETURN
    END IF
  END IF

! Order the data by increasing wavenumber.
  ALLOCATE(dummy(ss_data%n_wavenumber))
  ALLOCATE(map(ss_data%n_wavenumber))
  CALL map_heap_func(ss_data%wavenum(1:ss_data%n_wavenumber), map)
  DO k=1, ss_data%n_wavenumber
    dummy(k)=ss_data%wavenum(map(k))
  ENDDO
  ss_data%wavenum(1:ss_data%n_wavenumber)=dummy(1:ss_data%n_wavenumber)
  DO k=1, ss_data%n_wavenumber
    dummy(k)=ss_data%abs(map(k))
  ENDDO
  ss_data%abs(1:ss_data%n_wavenumber)=dummy(1:ss_data%n_wavenumber)
  DO k=1, ss_data%n_wavenumber
    dummy(k)=ss_data%scat(map(k))
  ENDDO
  ss_data%scat(1:ss_data%n_wavenumber)=dummy(1:ss_data%n_wavenumber)
  DO j=1, ss_data%n_phf_term
    DO k=1, ss_data%n_wavenumber
      dummy(k)=ss_data%phf(j, map(k))
    ENDDO
    ss_data%phf(j, 1:ss_data%n_wavenumber)=dummy(1:ss_data%n_wavenumber)
  ENDDO
  DEALLOCATE(dummy)
  DEALLOCATE(map)
!
!
!
END SUBROUTINE read_scatter_block_90

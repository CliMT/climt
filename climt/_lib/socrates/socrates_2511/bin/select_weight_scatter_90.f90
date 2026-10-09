! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to select weighting for scattering.
!
SUBROUTINE select_weight_scatter_90(i_weight, SolarSpec, t_weight, &
      i_method_weight, l_interactive, ierr)
!
! Description:
!   The weighting function and method of weighting for single
!   scattering propeties are selected.
!
! Method:
!   The weighting function for the monochromatic single
!   scattering properties is selected by the user from a
!   menu of possible functions. The user also sets the method 
!   of weighting, that is, thin or thick averaging.
!      
!
!
!
!
! Modules to set types of variables:
  USE realtype_rd
  USE dimensions_pp_ucf
  USE def_std_io_icf
  USE def_solarspec
  USE method_weight_pcf
  USE weighting_pcf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive use
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(OUT) :: i_weight
!   Type of weighting function
  TYPE (StrSolarSpec), Intent(OUT) :: SolarSpec
!   Solar spectral irradiance data
  REAL  (RealK), Intent(OUT) :: t_weight
!   Weighting temperature
  INTEGER, Intent(OUT) :: i_method_weight
!   Method of weighting
!
!
!
! Local variables
!
  INTEGER :: ios
!   I/O flag
  LOGICAL :: l_solar_spectrum
!       Logical for reading of spectrum
!
!
!
! Display menu of weighting functions.
  WRITE(iu_stdout, '(/a)') &
    'Select the type of weighting functions.'
  WRITE(iu_stdout, '(3x, i2, a)' ) IP_weight_planck, &
    '. Planckian weighting at a given temperature.'
  WRITE(iu_stdout, '(3x, i2, a)' ) IP_weight_d_planck, &
    '. Differential planckian weighting at a given temperature.'
  WRITE(iu_stdout, '(3x, i2, a)' ) IP_weight_solar, &
    '. TOA solar spectral weighting.'
  WRITE(iu_stdout, '(3x, i2, a)' ) IP_weight_uniform, &
    '. Uniform weighting.'
!
  DO
    READ(iu_stdin, *, iostat=ios) i_weight
    IF ( (ios /= 0)                         .OR.   &
         ( (i_weight /= IP_weight_planck)   .AND.  &
           (i_weight /= IP_weight_d_planck) .AND.  &
           (i_weight /= IP_weight_solar)    .AND.  &
           (i_weight /= IP_weight_uniform) ) ) then
      WRITE(iu_err, '(a)') '+++ Erroneous response:'
      IF (l_interactive) then
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
  IF ( (i_weight == IP_weight_planck) .OR.  &
       (i_weight == IP_weight_d_planck) ) THEN
    WRITE(iu_stdout, '(/a)') &
      'Specify the temperature for weighting.'
    DO
      READ(iu_stdin, *, iostat=ios) t_weight
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Unrecognized response: '
        IF (l_interactive) THEN
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
  ENDIF
!
! If solar weighting is used read the data in.
  IF (i_weight == ip_weight_solar) THEN
    CALL read_solar_spectrum(SolarSpec, ierr)
    IF (ierr /= i_normal) return
  ENDIF
!
! Display menu of methods of weighting.
  WRITE(iu_stdout, '(/a)') &
    'Select the method of weighting the optical properties.'
  WRITE(iu_stdout, '(3x, i2, a)' ) IP_weight_thin, &
    '. Simple averaging of single scattering properties.'
  WRITE(iu_stdout, '(3x, i2, a)' ) IP_weight_thick, &
    '. Thick averaging of the reflection.'
!
  WRITE(iu_stdout, '(a/)') 'Enter the required number.'
  DO
    READ(iu_stdin, *, IOSTAT=ios) i_method_weight
    IF (ios /= 0) THEN
      WRITE(iu_err, '(a)') '+++ Erroneous response:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE IF ( (i_method_weight /= IP_weight_thin) .AND. &
              (i_method_weight /= IP_weight_thick) ) THEN
      WRITE(iu_err, '(a)') '+++ response out of range:'
      IF (l_interactive) THEN
        WRITE(iu_stdout, '(a)') 'Please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
!
!
  RETURN
END SUBROUTINE select_weight_scatter_90

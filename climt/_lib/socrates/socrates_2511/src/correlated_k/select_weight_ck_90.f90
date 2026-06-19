! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the weighting function for correlated-k.
!
SUBROUTINE select_weight_ck_90 &
!
(i_weight, SolarSpec, l_interactive, ierr)
!
! Description:
!   A list of possible weighting functions is displayed and
!   the user selects one. If solar weighting is used the
!   solar spectrum is read in.
!
!
!
! Modules used:
  USE realtype_rd
  USE def_solarspec
  USE weighting_pcf
  USE error_pcf
  USE def_std_io_icf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(OUT) :: i_weight
!   Method of weighting
  TYPE (StrSolarSpec), Intent(OUT) :: SolarSpec
!   Solar spectral irradiance data
!
! Local variables
!
!
  INTEGER :: ios
!   I/O error flag
!
!
! Display menu of weightings.
  WRITE(iu_stdout, '(/a)') &
    'Select the method of weighting the transmittances.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_planck, &
    '. Planckian weighting at transmission temperature.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_d_planck, &
    '. Differential planckian weighting at transmission temperature.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_solar, &
    '. TOA solar spectral weighting.'
  WRITE(iu_stdout, '(3x, i2, a)' ) ip_weight_uniform, &
    '. Uniform weighting.'
!
  WRITE(iu_stdout, '(a/)') 'Enter required number.'
  DO
    READ(iu_stdin, *, iostat=ios) i_weight
    IF ( (ios /= 0)                         .OR.   &
         ( (i_weight /= ip_weight_planck)   .AND.  &
           (i_weight /= ip_weight_d_planck) .AND.  &
           (i_weight /= ip_weight_solar)    .AND.  &
           (i_weight /= ip_weight_uniform) ) ) then
      WRITE(iu_err, '(a)') '+++ Erroneous response:'
      IF (l_interactive) then
        WRITE(iu_stdout, '(a)') 'please re-enter.'
      ELSE
        ierr=i_err_fatal
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
  ENDDO
!
! If solar weighting is used read the data in.
  if (i_weight == ip_weight_solar) then
    CALL read_solar_spectrum(SolarSpec, ierr)
    if (ierr /= i_normal) RETURN
  endif
!
!
!
  RETURN
END SUBROUTINE select_weight_ck_90

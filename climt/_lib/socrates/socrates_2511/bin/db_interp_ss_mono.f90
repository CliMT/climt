! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to interpolate single scattering properties.
!
SUBROUTINE db_interp_ss_mono &
!
(nd_wavelength, nd_scatt_angle, nd_size_scat, &
 dimen, &
 l_phase, n_angle, &
 ice_db_mono_info, &
 scattering_point, extinction_point, asymmetry_point, i_stokes, &
 ierr)

  USE realtype_rd
  USE def_std_io_icf
  USE error_pcf
  USE def_db_ss_mono
  USE spline_evaluate_mod, ONLY: spline_evaluate
!
!
  IMPLICIT NONE
!
!
!
!
! Dummy arguments
!
! Sizes of arrays
  INTEGER, Intent(IN) :: nd_wavelength
!   Size allocated for wavelengths of scattering calculations
  INTEGER, Intent(IN) :: nd_scatt_angle
!   Size allocated for angles where scattering is carried out
  INTEGER, Intent(IN) :: nd_size_scat
!   Size allocated for the number of blocks of data in the database
!   at each frequency
!
  TYPE  (STR_db_ss_mono), Intent(IN) :: ice_db_mono_info
!   Total monochromatic information in the database
!
  REAL  (RealK), Intent(IN) :: dimen
!   Size (maximum dimension)
  INTEGER, Intent(IN) :: n_angle
!   Number of scattering angles
  LOGICAL, Intent(IN) :: l_phase
!   Calculate phase function
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
  REAL  (RealK), Intent(OUT) :: scattering_point 
!   Scattering (at a single wavelength)
  REAL  (RealK), Intent(OUT) :: extinction_point
!   Extinction (at a single wavelength)
  REAL  (RealK), Intent(OUT) :: asymmetry_point
!   Asymmetry factor (at a single wavelength)
  REAL  (RealK), Intent(OUT), Dimension(nd_scatt_angle) :: i_stokes
!   First Stokes parameter
!
!
!
! Local arguments.
  INTEGER  ::  i
!   Loop variable
!
!
!
! Check if input for size is within the range of the database.
  IF ( dimen < MINVAL(ice_db_mono_info%dm, &
               MASK=ice_db_mono_info%dm>0.0) .OR. &
       dimen > MAXVAL(ice_db_mono_info%dm) )  THEN
!
    scattering_point = 0.0_RealK
    extinction_point = 0.0_RealK
    asymmetry_point  = 0.0_RealK
    IF (l_phase) i_stokes(1:n_angle) = 0.0_RealK
!
    RETURN
!
 ENDIF
!
!
!
! Apply a spline fit to the quantities contained in the database
! as a function of size.
!
  CALL spline_evaluate(ierr, SIZE(ice_db_mono_info%dm), &
                          ice_db_mono_info%dm, &
                          ice_db_mono_info%ss(:, 1), &
                          ice_db_mono_info%d2_ss(:, 1), &
                          dimen, scattering_point)
  IF (ierr == i_err_range) THEN
!   Recover harmlessly form errors of interpolation.
    scattering_point = 0.0_RealK
    ierr =  i_normal
  ENDIF
  CALL spline_evaluate(ierr, SIZE(ice_db_mono_info%dm), &
                          ice_db_mono_info%dm, &
                          ice_db_mono_info%ss(:, 2), &
                          ice_db_mono_info%d2_ss(:, 2), &
                          dimen, extinction_point)
  IF (ierr == i_err_range) THEN
!   Recover harmlessly form errors of interpolation.
    extinction_point = 0.0_RealK
    ierr =  i_normal
  ENDIF
  CALL spline_evaluate(ierr, SIZE(ice_db_mono_info%dm), &
                          ice_db_mono_info%dm, &
                          ice_db_mono_info%ss(:, 3), &
                          ice_db_mono_info%d2_ss(:, 3), &
                          dimen, asymmetry_point)
  IF (ierr == i_err_range) THEN
!   Recover harmlessly form errors of interpolation.
    asymmetry_point = 0.0_RealK
    ierr =  i_normal
  ENDIF
!
  IF (l_phase) THEN
    DO i=1, n_angle
      CALL spline_evaluate(ierr, SIZE(ice_db_mono_info%dm), &
                              ice_db_mono_info%dm, &
                              ice_db_mono_info%ss(:, i+4), &
                              ice_db_mono_info%d2_ss(:, i+4), &
                              dimen, i_stokes(i))
      IF (ierr == i_err_range) THEN
!       Recover harmlessly form errors of interpolation.
        i_stokes(i) = 0.0_RealK
        ierr =  i_normal
      ENDIF
    ENDDO
  ENDIF
!
!
!
  RETURN
END SUBROUTINE db_interp_ss_mono

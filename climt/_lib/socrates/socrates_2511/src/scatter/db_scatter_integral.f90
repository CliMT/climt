! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE db_scatter_integral &
!
(nd_wavelength, nd_size_scat, &
 SizeDist,wavelength_index,n_wavelength,wavelength,DBGeom, &
 ice_db_mono_info, &
 n_angle,mu_angle, & 
 panel_ratio, &
 extinction, scattering, asymmetry, l_stokes, i_stokes, &
 nd_scatt_angle, &
 nd_integral, nd_panel, nd_refinement, &
 nd_point, nd_panel_point, &
 ierr &
)
!
! Method:
!   The extinction, the scattering and the asymmetry are 
!   initialized to 0. Initial estimates of these quantities are
!   made and the range of integration is extended until further
!   extension does not sensibly alter these estimates. This is
!   done by dividing the range of integration into panels and
!   adding new panels as required. the resolution within a panel 
!   is then increased until the integrals converge.
!
!
!
! Modules used
  USE realtype_rd
  USE def_size_dist  
  USE def_std_io_icf  
  USE prec_integral_tcf  
  USE scatter_algorithm_pcf  
  USE shape_particle_pcf  
  USE error_pcf  
  USE rad_ccf, ONLY: pi  
  USE def_db_crystal_geometry  
  USE def_db_ss_mono  
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments
!
  INTEGER, Intent(IN) :: n_wavelength 
!           Number of wavelengths
!
! Sizes of dummy arrays:
  INTEGER, Intent(IN) :: nd_wavelength
!   Size allocated for wavelengths of scattering calculations
  INTEGER, Intent(IN) :: nd_scatt_angle
!   Size allocated for angles where scattering is carried out
  INTEGER, Intent(IN) :: nd_size_scat
!   Size allocated for the number of blocks of data in the database
!   at each frequency
  INTEGER, Intent(IN) :: nd_integral
!           Size allocated for integrals to be evaluated concurrently
  INTEGER, Intent(IN) :: nd_panel
!           Size allocated for panels in initial integration
  INTEGER, Intent(IN) :: nd_refinement
!           Size allocated for number of refinements
  INTEGER, Intent(IN) :: nd_point
!           Size allocated for points in the interval of integration
  INTEGER, Intent(IN) :: nd_panel_point
!           Size allocated for points in each panel
!
  TYPE (STR_size_dist), Intent(IN) :: SizeDist
!       Size distribution
!
  TYPE (STR_db_cryst_geom), Intent(IN) :: DBGeom
!       Ice crystal geometry for use with databases
  TYPE (STR_db_ss_mono), Intent(IN) :: ice_db_mono_info
!       All information in the database at this wavelength
!
!
!
  LOGICAL, INTENT(IN) :: l_stokes
!           Flag to calculate Stokes''s first parameter
  INTEGER, INTENT(IN) :: n_angle
!           Number of scattering angles
  INTEGER, INTENT(IN) :: wavelength_index
!           Wavelength of light (database index)
  REAL (RealK), INTENT(IN) :: wavelength(nd_wavelength)
!
  REAL  (RealK), INTENT(IN) :: mu_angle(nd_scatt_angle)
!           Cosines of scattering angles
  REAL  (RealK), INTENT(IN) :: panel_ratio
!           Ratio of sizes for expansion or contraction of the panels
!
  REAL  (RealK), INTENT(OUT) :: extinction
!           Mean extinction
  REAL  (RealK), Intent(OUT) :: scattering
!           Mean scattering
  REAL  (RealK), Intent(OUT) :: asymmetry
!           Mean asymmetry
  REAL  (RealK), Intent(OUT) :: i_stokes(nd_scatt_angle)
!           Mean value of Stokes''s first parameter
!
  INTEGER, Intent(INOUT) :: ierr
!           Error flag
!
!
!     Local variables.
  INTEGER :: n_integral
!           Number of integrals to evaluate
  INTEGER :: n_panel
!           Number of panels
  INTEGER :: i_panel
!           Index of panel
  INTEGER :: i
!           Loop variable
  INTEGER :: j
!           Loop variable
  INTEGER :: K
!           Loop variable
  INTEGER :: n_point
!           Number of points in interval
  INTEGER :: i_refinement
!           Index of refinement
  LOGICAL :: l_add_upper
!           Logical to add an upper panel
  LOGICAL :: l_add_lower
!           Logical to add a lower panel
  LOGICAL :: l_refine
!           Logical to perform refinement
  LOGICAL :: l_phase
!           Flag to calculate the phase function
  REAL  (RealK) :: dimen_low
!           Lower limit of integrals
  REAL  (RealK) :: dimen_high
!           Upper limit of integrals
  REAL  (RealK) :: dimen_panel_low
!           Lower limit of panel
  REAL  (RealK) :: dimen_panel_high
!           Upper limit of panel
  REAL  (RealK) :: integral_estimate(nd_integral)
!           Estimates of integrals
  REAL  (RealK) :: integral_panel_temp(nd_integral)
!           Estimates of panel integrals
  REAL  (RealK) :: integral_panel(nd_integral, nd_panel, nd_refinement)
!           Integrals over panels
  REAL  (RealK) :: dimen(nd_point)
!           Radii of points of integration
  REAL  (RealK) :: y(nd_point, nd_integral)
!           Point values of integrand
  REAL  (RealK) :: dimen_coarse(nd_panel, nd_panel_point)
!           Starting values for romberg step
  REAL  (RealK) :: y_coarse(nd_integral, nd_panel, nd_panel_point)
!           Starting values for romberg step
  REAL  (RealK) :: refinement(nd_integral)
!           Value of refinement
!
  REAL  (RealK) :: number_point
!           Point value of number density
  REAL  (RealK) :: proj_area_point
!           Point value of projected area
  REAL  (RealK) :: asymmetry_point
!           Point value of asymmetry
  REAL  (RealK) :: volume_point
!           Point value of volume
  REAL  (RealK) :: extinction_point
!           Point value of extinction
  REAL  (RealK) :: scattering_point
!           Point value of scattering
  REAL  (RealK) :: s(nd_scatt_angle)
!           Element of scattering function
!
! Functions called:
  REAL  (RealK) :: number_particle_90
!       Function for number density
  REAL  (RealK) :: trapezoid
!       Integration routine

  EXTERNAL number_particle_90, &
           trapezoid

  INTERFACE

     FUNCTION proj_area_particle &
              (dimen, i_measure, i_shape, ierr, DBGeom) RESULT(proj_area)

        USE realtype_rd
        USE def_db_crystal_geometry

        INTEGER, Intent(IN) :: i_shape
        INTEGER, Intent(IN) :: i_measure
        REAL  (RealK), Intent(IN) :: dimen
        TYPE  (STR_db_cryst_geom), Intent(IN), Optional :: DBGeom
        REAL  (RealK) :: proj_area
        INTEGER, Intent(InOut) :: ierr

     END FUNCTION proj_area_particle

     FUNCTION volume_particle &
              (dimen, i_measure, i_shape, ierr, DBGeom) RESULT(volume)

        USE realtype_rd
        USE def_db_crystal_geometry

        INTEGER, Intent(IN) :: i_shape
        INTEGER, Intent(IN) :: i_measure
        REAL  (RealK), Intent(IN) :: dimen
        TYPE  (STR_db_cryst_geom), Intent(IN), Optional :: DBGeom
        REAL  (RealK) :: volume
        INTEGER, Intent(InOut) :: ierr

     END FUNCTION volume_particle
  END INTERFACE

!
!
!
  extinction = 0.0_RealK
  scattering = 0.0_RealK
  asymmetry = 0.0_RealK
  IF (l_stokes) i_stokes(1:n_angle) = 0.0_RealK
!
! Use a variable for the number of integrals to facilitate
! subsequent modification.
  n_integral = 3
  IF (l_stokes) THEN
    n_integral = n_integral + n_angle
  ENDIF
  IF (n_integral > nd_integral) THEN
    WRITE(iu_err, '(/A)') &
      '*** Error: The number of integrals to be evaluated ' &
      //'exceeds the allocated size.'
    ierr = i_err_fatal
    RETURN
  ENDIF
!
! The phase function is needed if Stokes''s parameter is.
  l_phase = l_stokes
!
!
!
! Preliminary integration
  integral_estimate(1:n_integral) = 0.0_RealK
  dimen_low   = SizeDist%dimen_0
  dimen_high  = SizeDist%dimen_0
  n_panel     = 0
  l_add_upper = .TRUE.
  l_add_lower = .TRUE.
!
! Include extra panels in the range of integration at either limit
! until the inclusion of further panels makes little change in the
! estimate of the integral.
!
  Set_range: DO ; IF ( .NOT.(l_add_upper .OR. l_add_lower) ) EXIT
!
    IF (l_add_upper) THEN
      dimen_panel_low = dimen_high
      dimen_panel_high = panel_ratio * dimen_high
      DO i = 1, nd_panel_point
!
        dimen(i) = dimen_panel_low * &
          EXP((REAL(i-1, RealK) / REAL(nd_panel_point-1, RealK)) * &
          LOG(dimen_panel_high / dimen_panel_low))
!
        number_point    = number_particle_90(dimen(i), SizeDist, ierr)
        proj_area_point = proj_area_particle(dimen(i), &
                                             SizeDist%i_measure, &
                                             SizeDist%i_shape_particle, &
                                             ierr, DBGeom)
        volume_point    = volume_particle   (dimen(i), &
                                             SizeDist%i_measure, &
                                             SizeDist%i_shape_particle, &
                                             ierr, DBGeom)
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF
!
!
        IF (number_point > 0.0_RealK) THEN
!         The scattering code is not called if there are
!         no particles at this size. 
          CALL db_interp_ss_mono(nd_wavelength, &
                                 nd_scatt_angle, nd_size_scat, &
                                 dimen(i), &
                                 l_phase, n_angle, &
                                 ice_db_mono_info, &
                                 scattering_point, &
                                 extinction_point, &
                                 asymmetry_point, &
                                 s, &
                                 ierr)
          IF (ierr /= i_normal) THEN
            WRITE(iu_err, '(/A, 3(/A23, 1PE12.5))') &
              '*** Error: Failure in database scattering calculations: ', &
              '      Radius           =  ', dimen(i),      &
              '      Wavelength_index =  ', real(wavelength_index),    &
              '      Number density   =  ', number_point
            RETURN
          ENDIF
        ELSE
          scattering_point = 0.0_RealK
          extinction_point = 0.0_RealK
          asymmetry_point = 0.0_RealK
          IF (l_phase) s(1:n_angle) = 0.0_RealK
        ENDIF
        IF (l_stokes) &
          y(i, 4:3+n_angle) = number_point * s(1:n_angle)
!
!
        y(i, 1) = number_point * extinction_point
        y(i, 2) = number_point * scattering_point
        y(i, 3) = asymmetry_point * y(i, 2)
!
      ENDDO
!

!     The panel is included if any one of the integrals is of
!     sufficient size.
      l_add_upper = .FALSE.
      DO j = 1, n_integral
        integral_panel_temp(j) = trapezoid(nd_panel_point &
          , dimen, y(1, j))
        l_add_upper = l_add_upper .OR. &
          ( ABS(integral_panel_temp(j)) > tol_panel * &
                                          ABS(integral_estimate(j)) )
      ENDDO
!
      IF (l_add_upper) THEN
        IF (n_panel >= nd_panel) THEN
          WRITE(iu_err, '(/A, /A)') &
            '*** Error: Too many panels are required to evaluate' &
            , 'the integrals reliably: increase npd_panel.'
          ierr = i_err_fatal
          RETURN
        ENDIF
        n_panel = n_panel + 1
        dimen_high = dimen_high * panel_ratio
        DO i = 1, nd_panel_point
          dimen_coarse(n_panel, i) = dimen(i)
        ENDDO
        DO j = 1, n_integral
          integral_estimate(j) = integral_estimate(j) + &
            integral_panel_temp(j)
          integral_panel(j, n_panel, 1) = integral_panel_temp(j)
          DO i = 1, nd_panel_point
            y_coarse(j, n_panel, i) = y(i, j)
          ENDDO
        ENDDO
      ENDIF
!
    ENDIF
!
!
    IF (l_add_lower) THEN
      dimen_panel_high = dimen_low
      dimen_panel_low  = dimen_panel_high / panel_ratio
      DO i = 1, nd_panel_point

        dimen(i) = dimen_panel_low &
           * EXP((REAL(i-1, RealK)/REAL(nd_panel_point-1, RealK)) &
           * LOG(dimen_panel_high/dimen_panel_low))

        number_point = number_particle_90(dimen(i), SizeDist, ierr)
        proj_area_point = proj_area_particle(dimen(i), &
                                             SizeDist%i_measure, &
                                             SizeDist%i_shape_particle, &
                                             ierr, DBGeom)
        volume_point    = volume_particle   (dimen(i), &
                                             SizeDist%i_measure, &
                                             SizeDist%i_shape_particle, &
                                             ierr, DBGeom)
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF

        IF (number_point > 0.0_RealK) THEN
!         The scattering code is not called if there are
!         no particles at this size. 
          CALL db_interp_ss_mono(nd_wavelength, &
                                 nd_scatt_angle, nd_size_scat, &
                                 dimen(i), &
                                 l_phase, n_angle, &
                                 ice_db_mono_info, &
                                 scattering_point, &
                                 extinction_point, &
                                 asymmetry_point, &
                                 s, &
                                 ierr)
          IF (ierr /= i_normal) THEN
            WRITE(iu_err, '(/A, 3(/A23, 1PE12.5))') &
              '*** Error: Failure in database scattering calculations: ', &
              '      Radius           =  ', dimen(i),      &
              '      Wavelength_index =  ', real(wavelength_index),    &
              '      Number density   =  ', number_point
            RETURN
          ENDIF
        ELSE
          scattering_point = 0.0_RealK
          extinction_point = 0.0_RealK
          asymmetry_point  = 0.0_RealK
          IF (l_phase) s(1:n_angle) = 0.0_RealK
        ENDIF
        IF (l_stokes) &
          y(i, 4:3+n_angle) = number_point * s(1:n_angle)
!
!
        y(i, 1) = number_point * extinction_point 
        y(i, 2) = number_point * scattering_point
        y(i, 3) = asymmetry_point * y(i, 2)
!
      ENDDO
      DO j = 1, n_integral
        integral_panel_temp(j) = trapezoid(nd_panel_point, &
          dimen, y(1, j))
      ENDDO
!     The panel is included if any one of the integrals is of
!     sufficient size.
      l_add_lower = .FALSE.
      DO j = 1, n_integral
        integral_panel_temp(j) = trapezoid(nd_panel_point, &
          dimen, y(1, j))
        l_add_lower = l_add_lower .OR. &
          (ABS(integral_panel_temp(j)) > tol_panel * &
                                         ABS(integral_estimate(j)))
      ENDDO

      IF (l_add_lower) THEN
        IF (n_panel >= nd_panel) THEN
          WRITE(iu_err, '(/A, /A)') &
            '*** Error: Too many panels are required to evaluate' &
            , 'the integrals reliably: increase npd_panel.'
          ierr = i_err_fatal
          RETURN
        ENDIF
        n_panel = n_panel + 1
        dimen_low = dimen_low/panel_ratio
        DO i = 1, nd_panel_point
          dimen_coarse(n_panel, i) = dimen(i)
        ENDDO
        DO j = 1, n_integral
          integral_estimate(j) = integral_estimate(j) + &
            integral_panel_temp(j)
          integral_panel(j, n_panel, 1) = integral_panel_temp(j)
          DO i = 1, nd_panel_point
            y_coarse(j, n_panel, i) = y(i, j)
          ENDDO
        ENDDO
      ENDIF
!
    ENDIF
!
  ENDDO Set_range
!
!
! Now use Romberg extrapolation to determine the more exact values
! of the integrals.
  DO i_panel = 1, n_panel
!
    i_refinement = 1
    n_point = nd_panel_point
    DO i = 1, n_point
      dimen(i) = dimen_coarse(i_panel, i)
      DO j = 1, n_integral
        y(i, j) = y_coarse(j, i_panel, i)
      ENDDO
    ENDDO
!
    DO
      i_refinement = i_refinement + 1
!     Add new points between old ones.
      DO i = n_point, 1, -1
        dimen(2 * i-1) = dimen(i)
        DO j = 1, n_integral
          y(2 * i-1, j) = y(i, j)
        ENDDO
      ENDDO
      n_point = 2 * n_point - 1
      DO i = 2, n_point-1, 2

        dimen(i) = SQRT(dimen(i + 1) * dimen(i - 1))
        number_point = number_particle_90(dimen(i), SizeDist, ierr)
        proj_area_point = proj_area_particle(dimen(i), &
                                             SizeDist%i_measure, &
                                             SizeDist%i_shape_particle, &
                                             ierr, DBGeom)
        volume_point    = volume_particle   (dimen(i), &
                                             SizeDist%i_measure, &
                                             SizeDist%i_shape_particle, &
                                             ierr, DBGeom)
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF

        IF (number_point > 0.0_RealK) THEN
!         The scattering code is not called if there are
!         no particles at this size.
          CALL db_interp_ss_mono(nd_wavelength, &
                                 nd_scatt_angle, nd_size_scat, &
                                 dimen(i), &
                                 l_phase, n_angle, &
                                 ice_db_mono_info, &
                                 scattering_point, &
                                 extinction_point, &
                                 asymmetry_point, &
                                 s, &
                                 ierr)
          IF (ierr /= i_normal) THEN
            WRITE(iu_err, '(/A, 3(/A23, 1PE12.5))') &
              '*** Error: Failure in database scattering calculations: ' &
              , '      Radius          =  ', dimen(i) &
              , '      Wavelength      =  ', real(wavelength_index) &
              , '      Number density  =  ', number_point
            RETURN
          ENDIF
        ELSE
          scattering_point = 0.0_RealK
          extinction_point = 0.0_RealK
          asymmetry_point = 0.0
          IF (l_phase) s(1:n_angle) = 0.0_RealK
        ENDIF
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF
        IF (l_stokes) &
          y(i, 4:3+n_angle) = number_point * s(1:n_angle)
!
        y(i, 1) = number_point * extinction_point
        y(i, 2) = number_point * scattering_point
        y(i, 3) = asymmetry_point * y(i, 2)
!
      ENDDO
!
      l_refine = .FALSE.
      DO j = 1, n_integral
!       In relative terms there is little to be gained by using
!       the earlier value of the integral for the sum
!       over alternate points.
        refinement(j) = trapezoid(n_point, dimen, y(1, j)) - &
          integral_panel(j, i_panel, i_refinement-1)
        integral_panel(j, i_panel, i_refinement) = &
          integral_panel(j, i_panel, i_refinement-1) + &
          refinement(j) / 0.75_RealK
        l_refine = l_refine .OR. &
          (ABS(refinement(j)) > tol_refinement * &
                                ABS(integral_estimate(j)))
      ENDDO
      IF (l_refine) THEN
        IF ( (2 * n_point >= nd_point) .OR. &
             (i_refinement >= nd_refinement) ) THEN
          WRITE(iu_err, '(/A, 1X, 1PE10.3, 1X, A)') &
            '+++ Warning: Scattering integrals at wavelength index' , &
            wavelength_index, 'have not converged.'
          WRITE(iu_err, '(A)') 'The last values will be used.'
        EXIT
        ENDIF
      ELSE
          EXIT
      ENDIF
!
    ENDDO
!
    extinction = extinction + &
      integral_panel(1, i_panel, i_refinement)
    scattering = scattering + &
      integral_panel(2, i_panel, i_refinement)
    asymmetry = asymmetry + &
      integral_panel(3, i_panel, i_refinement)
    IF (l_stokes) &
      i_stokes(1:n_angle) = i_stokes(1:n_angle) + &
        integral_panel(4:n_angle+3, i_panel, i_refinement)
!
  ENDDO
!
  asymmetry = asymmetry / scattering
!
!
!
  RETURN
END SUBROUTINE db_scatter_integral

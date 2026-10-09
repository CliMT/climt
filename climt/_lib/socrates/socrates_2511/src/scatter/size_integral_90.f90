! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate properties of a size distribution.
!
Subroutine size_integral_90 &
!
(nd_integral, nd_panel, nd_refinement, nd_point, &
 nd_panel_point, panel_ratio, SizeDist, &
 number_total, proj_area_total, volume_total, dimen_average, &
 ierr, DBGeom &
)
!
! Description:
!   This subroutine performs integration of geometrical information
!   over the size distribution.
!
! Method:
!   The total number density, the projected area and the volume 
!   fraction are initialized to 0. Initial estimates of these 
!   quantities are made and the range of integration is extended 
!   until further extension does not sensibly alter these estimates. 
!   This is done by dividing the range of integration into panels 
!   and adding new panels as required. The resolution within a panel
!   is then increased until the integrals converge.
!
!
!
! Modules used.
  USE realtype_rd
  USE def_std_io_icf
  USE shape_particle_pcf
  USE def_size_dist
  USE prec_integral_tcf
  USE error_pcf
  USE rad_ccf, ONLY: pi
  USE def_db_crystal_geometry
!
!
  IMPLICIT NONE
!
!
! Dummy variables:
  INTEGER, Intent(IN) :: nd_integral
!       Size allocated for integrals to be evaluated concurrently
  INTEGER, Intent(IN) :: nd_panel
!       Size allocated for panels in initial integration
  INTEGER, Intent(IN) :: nd_refinement
!       Size allocated for number of refinements
  INTEGER, Intent(IN) :: nd_point
!       Size allocated for points in the interval of integration
  INTEGER, Intent(IN) :: nd_panel_point
!       Size allocated for points in each panel
!
  TYPE (STR_size_dist), Intent(IN) :: SizeDist
!       Size distribution
!
!
!
  REAL  (RealK), Intent(IN) :: panel_ratio
!         Ratio of sizes for expansion of the panels
!
  REAL  (RealK), INTENT(OUT) :: number_total
!         Total number of particles per unit volume
  REAL  (RealK), INTENT(OUT) :: proj_area_total
!         Total projected area per unit volume
  REAL  (RealK), INTENT(OUT) :: volume_total
!         Total volume of particles per unit volume
  REAL  (RealK), INTENT(OUT) :: dimen_average
!
  INTEGER, Intent(InOut) :: ierr
!           Error flag
!
  TYPE  (STR_db_cryst_geom), Intent(IN), Optional :: DBGeom
!   Ice crystal geometry for use with databases
!
!
!
! Local variables.
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
  INTEGER :: n_point
!           Number of points in interval
  INTEGER :: i_refinement
!           Index of refinement
  LOGICAL :: l_dimen_average
!           Flag for evaluation of mean dimension
  LOGICAL :: l_add_upper
!           Logical to add an upper panel
  LOGICAL :: l_add_lower
!           Logical to add a lower panel
  LOGICAL :: l_refine
!           Logical to perform refinement
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
  REAL  (RealK) :: number_point
!           Point value of number density
!
!     Functions called:
  REAL  (RealK) :: number_particle_90
!           Function for number density
  REAL  (RealK) :: trapezoid
!           Integration routine
!
  EXTERNAL number_particle_90, &
           
           trapezoid

! External Functions and Subroutines
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
        INTEGER, Intent(InOut) :: ierr
        TYPE  (STR_db_cryst_geom), Intent(IN), Optional :: DBGeom
        REAL  (RealK) :: volume

     END FUNCTION volume_particle
  END INTERFACE

!  
! - End Of Header
!
!
!
! Use a variable number of integrals to allow for expansion.
  n_integral = nd_integral
! Determine whether the average dimension is to be set from the
! number of integrals: this would be more elegantly done using
! optional arguments, but that requires interface blocks.
  l_dimen_average = (n_integral==4)
!
  number_total       = 0.0_RealK
  volume_total       = 0.0_RealK
  proj_area_total    = 0.0_RealK
  IF (l_dimen_average) dimen_average = 0.0_RealK
!
! Preliminary integration
  DO j = 1, n_integral
    integral_estimate(j) = 0.0_RealK
  ENDDO
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
  DO ; IF ( .NOT.(l_add_upper .OR. l_add_lower) ) EXIT
!
    IF (l_add_upper) THEN
      dimen_panel_low  = dimen_high
      dimen_panel_high = panel_ratio*dimen_high
      DO i = 1, nd_panel_point
        dimen(i) = dimen_panel_low * &
          EXP( (REAL(i-1, RealK) / REAL(nd_panel_point-1, RealK) ) * &
          LOG(dimen_panel_high / dimen_panel_low) )
        number_point = number_particle_90(dimen(i), SizeDist, ierr)
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF
        y(i, 1) = number_point * proj_area_particle(dimen(i), &
          SizeDist%i_measure, &
          SizeDist%i_shape_particle, &
          ierr, DBGeom)
        y(i, 2) = number_point * volume_particle(dimen(i), &
          SizeDist%i_measure, &
          SizeDist%i_shape_particle, &
          ierr, DBGeom)
        y(i, 3) = number_point 
        IF (l_dimen_average) y(i, 4) = number_point * dimen(i)
      ENDDO
!     The panel is included if any one of the integrals is of
!     sufficient size.
      l_add_upper = .FALSE.
      DO j = 1, n_integral
        integral_panel_temp(j) = trapezoid(nd_panel_point, &
          dimen, y(1, j))
        l_add_upper = l_add_upper .OR. &
          ( ABS(integral_panel_temp(j)) > tol_panel * &
                                          ABS(integral_estimate(j)) )
      ENDDO
!     No more panels if limit is reached
      l_add_upper = l_add_upper .AND. (n_panel < nd_panel-1)

      IF (l_add_upper) THEN
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
      dimen_panel_low = dimen_panel_high/panel_ratio
      DO i = 1, nd_panel_point
        dimen(i) = dimen_panel_low * &
          EXP( (REAL(i-1, RealK) / REAL(nd_panel_point-1, RealK)) * &
          LOG( dimen_panel_high / dimen_panel_low) )
        number_point = number_particle_90(dimen(i), SizeDist, ierr)
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF
        y(i, 1) = number_point * proj_area_particle(dimen(i), &
          SizeDist%i_measure, &
          SizeDist%i_shape_particle, &
          ierr, DBGeom)
        y(i, 2) = number_point * volume_particle(dimen(i), &
          SizeDist%i_measure, &
          SizeDist%i_shape_particle, &
          ierr, DBGeom)
        y(i, 3) = number_point
        IF (l_dimen_average) y(i, 4) = number_point * dimen(i)
      ENDDO
      DO j = 1, n_integral
        integral_panel_temp(j) = &
          trapezoid(nd_panel_point, dimen, y(1, j))
      ENDDO
!     The panel is included if any one of the integrals is 
!     of sufficient size.
      l_add_lower = .FALSE.
      DO j = 1, n_integral
        integral_panel_temp(J) = trapezoid(nd_panel_point, dimen, y(1, j))
        l_add_lower = l_add_lower .OR. &
          (ABS(integral_panel_temp(j)) > tol_panel * &
                                         ABS(integral_estimate(j)))
      ENDDO
!     No more panels if limit is reached
      l_add_lower = l_add_lower .AND. (n_panel < nd_panel-1)

      IF (l_add_lower) THEN
        n_panel = n_panel + 1
        dimen_low = dimen_low / panel_ratio
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
  ENDDO
!
!
! Now use Romberg extrapolation to determine the exact value
! of the integral.
  DO i_panel = 1, n_panel
!
    i_refinement = 1
    n_point = nd_panel_point
    DO i = 1, n_point
      dimen(i) = dimen_coarse(i_panel, i)
      DO j = 1, n_integral
        y(i, j) = y_coarse(J, i_panel, i)
      ENDDO
    ENDDO
!
    DO
      i_refinement = i_refinement + 1
!     Add new points between old ones.
      DO i  =  n_point, 1, -1
        dimen(2*i-1) = dimen(i)
        DO j = 1, n_integral
          y(2*i-1, j) = y(i, j)
        ENDDO
      ENDDO
      n_point = 2*n_point - 1
      DO i = 2, n_point-1, 2
        dimen(i) = SQRT(dimen(i+1) * dimen(i-1))
        number_point = number_particle_90(dimen(i), SizeDist, ierr)
        IF (ierr /= i_normal) THEN
          ierr = i_err_fatal
          RETURN
        ENDIF
        y(i, 1) = number_point * proj_area_particle(dimen(i), &
          SizeDist%i_measure, &
          SizeDist%i_shape_particle, &
          ierr, DBGeom)
        y(i, 2) = number_point * volume_particle(dimen(i), &
          SizeDist%i_measure, &
          SizeDist%i_shape_particle, &
          ierr, DBGeom)
        y(i, 3) = number_point
        IF (l_dimen_average) y(i, 4) = number_point * dimen(i)
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
        IF ( (2*n_point >= nd_point) .OR. &
             (i_refinement >= nd_refinement) ) THEN
          WRITE(iu_err, '(/A)') &
            '+++ Warning: The integral has not converged.'
          WRITE(iu_err, '(A)') 'The last value will be used.'
          EXIT
        ENDIF
      ELSE
        EXIT
      ENDIF
!
    ENDDO
!

    proj_area_total = proj_area_total + &
      integral_panel(1, i_panel, i_refinement)
    volume_total = volume_total + &
      integral_panel(2, i_panel, i_refinement)
    number_total = number_total + &
      integral_panel(3, i_panel, i_refinement)
    IF (l_dimen_average) THEN
      dimen_average = dimen_average + &
        integral_panel(4, i_panel, i_refinement)
    ENDIF

   ENDDO
!
   IF (l_dimen_average) THEN
     IF (number_total > epsilon(number_total)) THEN
       dimen_average = dimen_average/number_total
     ELSE
       dimen_average = 0.0
     ENDIF
   ENDIF
!
!
  RETURN
END

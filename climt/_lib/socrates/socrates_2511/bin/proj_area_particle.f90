! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate projected area of particles.
!
FUNCTION proj_area_particle &
!
(dimen, i_measure, i_shape, ierr, DBGeom)
!
! Description:
!   This routine returns the projected area of a particle 
!   as evaluated from the shape.
!
! Modules used.
  USE realtype_rd
  USE shape_particle_pcf
  USE measure_particle_pcf
  USE def_std_io_icf
  USE def_db_crystal_geometry
  USE error_pcf
  USE rad_ccf, ONLY: pi
  USE spline_evaluate_mod, ONLY: spline_evaluate
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  INTEGER, Intent(IN) :: i_shape
!           Shape of particle
  INTEGER, Intent(IN) :: i_measure
!           Mode of measurement of particle
  REAL  (RealK), Intent(IN) :: dimen
!           Dimension for evaluation
  TYPE  (STR_db_cryst_geom), Intent(IN), Optional :: DBGeom
!           Geometrical information for using a database
!
  REAL  (RealK) :: proj_area_particle
!           Projected area of particle
  INTEGER, Intent(InOut) :: ierr
!           Error flag
!
! Local variables.
  INTEGER :: i_size_range
!           Region, in range of sizes, where the crystal lies
!
!
! Coeffieicents in power laws for projected area
  REAL (RealK), Dimension(0: npd_shape, 2) :: Sigma
!   Coefficient in power law for projected area
  REAL (RealK), Dimension(0: npd_shape, 2) :: Delta
!   Exponents in power laws for projected area
  REAL (RealK), Dimension(0: npd_shape)    :: Dtrans
!   Transitional dimensions
  REAL (RealK), Dimension(0: npd_shape)    :: Dscale
!   Scaling dimensions
!
!- End of Header
!
  Dtrans(ip_shape_null)        = 1.0E-04_RealK
  Dtrans(ip_shape_sphere)      = 1.0E-04_RealK
  Dtrans(ip_shape_hexcyl)      = 1.0E-04_RealK
  Dtrans(ip_shape_polycrystal) = 1.0E-04_RealK
  Dtrans(ip_shape_plate)       = 1.0E-04_RealK
  Dtrans(ip_shape_rosette)     = 1.0E-04_RealK
  Dtrans(ip_shape_column)      = 1.0E-04_RealK
!
  Dscale(ip_shape_null)        = 1.0E-06_RealK
  Dscale(ip_shape_sphere)      = 1.0E-06_RealK
  Dscale(ip_shape_hexcyl)      = 1.0E-06_RealK
  Dscale(ip_shape_polycrystal) = 1.0E-06_RealK
  Dscale(ip_shape_plate)       = 1.0E-06_RealK
  Dscale(ip_shape_rosette)     = 1.0E-06_RealK
  Dscale(ip_shape_column)      = 1.0E-06_RealK
!
  Sigma(ip_shape_null, :)        = (/ 0.0_RealK, 0.0_RealK /)
  Sigma(ip_shape_sphere, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Sigma(ip_shape_hexcyl, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Sigma(ip_shape_polycrystal, :) = (/ 2.285E-05_RealK, 2.285E-05_RealK /)
  Sigma(ip_shape_plate, :)       = (/ 2.395E-05_RealK, 2.395E-05_RealK /)
  Sigma(ip_shape_rosette, :)     = (/ 6.837E-05_RealK, 8.687E-06_RealK /)
  Sigma(ip_shape_column, :)      = (/ 6.837E-05_RealK, 4.590E-06_RealK /)
!
  Delta(ip_shape_null, :)        = (/ 0.0_RealK, 0.0_RealK /)
  Delta(ip_shape_sphere, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Delta(ip_shape_hexcyl, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Delta(ip_shape_polycrystal, :) = (/ 1.880_RealK, 1.880_RealK /)
  Delta(ip_shape_plate, :)       = (/ 1.855_RealK, 1.855_RealK /)
  Delta(ip_shape_rosette, :)     = (/ 2.000_RealK, 1.568_RealK /)
  Delta(ip_shape_column, :)      = (/ 2.000_RealK, 1.415_RealK /)
!
!
!
  IF (i_shape == ip_shape_sphere) THEN
!
    IF (i_measure == ip_measure_radius) THEN
      proj_area_particle = pi * dimen**2
    ELSE IF (i_measure == ip_measure_max_dimen) THEN
      proj_area_particle = 0.25_RealK * pi * dimen**2
    ELSE IF (i_measure == ip_measure_proj_radius) THEN
      proj_area_particle = pi * dimen**2
    ENDIF
!
  ELSE IF ( (i_shape == ip_shape_column)  .OR. &
            (i_shape == ip_shape_plate)   .OR. &
            (i_shape == ip_shape_rosette) .OR. &
            (i_shape == ip_shape_polycrystal) ) THEN
!
    IF (i_measure == ip_measure_max_dimen) THEN
      IF (dimen < Dtrans(i_shape) ) THEN
        i_size_range = 1
      ELSE
        i_size_range = 2
      ENDIF
      proj_area_particle = Sigma(i_size_range, i_shape) * &
                             (dimen/Dscale(i_shape))  ** &
                             Delta(i_size_range, i_shape)
    ELSE IF (i_measure == ip_measure_proj_radius) THEN
      WRITE(iu_err, '(/A)') '*** Error: Incomplete code.'
      ierr = i_err_fatal
      RETURN
    ELSE
      WRITE(iu_err, '(/A)') '*** Error: Illegal measure of size.'
      ierr = i_err_fatal
      RETURN
    ENDIF
!
  ELSE IF ( (i_shape == ip_shape_db_undefined) ) THEN

    IF (i_measure == ip_measure_max_dimen) THEN
!     Interpolate in the square of the maximum dimension.
      CALL spline_evaluate(ierr, DBGeom%n_geom, &
             DBGeom%dm2, DBGeom%proj_area, DBGeom%d2_proj_area, &
             dimen**2, proj_area_particle)
      IF (ierr == i_err_range) THEN
!       Recover harmlessly form errors of interpolation.
        proj_area_particle = 0.0_RealK
        ierr =  i_normal
      ENDIF
    ELSE IF (i_measure == ip_measure_radius) THEN
      proj_area_particle = pi * dimen**2
    ELSE IF (i_measure == ip_measure_proj_radius) THEN
      proj_area_particle = pi * dimen**2
    ELSE
      WRITE(iu_err, '(/A)') '*** Error: Illegal measure of size.'
      ierr = i_err_fatal
      RETURN
    ENDIF

  ENDIF
!
!
!
  RETURN
END

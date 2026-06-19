! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate volume of particles.
FUNCTION volume_particle &
!
(dimen, i_measure, i_shape, ierr, DBGeom)
!
! Description:
!   This routine returns the volume of a particle 
!   as evaluated from the shape.

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
  REAL  (RealK) :: volume_particle
!           Projected area of particle
  INTEGER, Intent(InOut) :: ierr
!           Error flag
!
!
! Local variables.
  INTEGER :: i_size_range
!           Region, in range of sizes, where the crystal lies
!
!
! Coeffieicents in power laws for volume
  REAL (RealK), Dimension(0: npd_shape, 2) :: Alpha
!   Coefficient in power law for volume
  REAL (RealK), Dimension(0: npd_shape, 2) :: Beta
!   Exponents in power laws for volume
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
  Alpha(ip_shape_null, :)        = (/ 0.0_RealK, 0.0_RealK /)
  Alpha(ip_shape_sphere, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Alpha(ip_shape_hexcyl, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Alpha(ip_shape_polycrystal, :) = (/ 5.801E-04_RealK, 7.389E-06_RealK /)
  Alpha(ip_shape_plate, :)       = (/ 4.953E-05_RealK, 7.389E-06_RealK /)
  Alpha(ip_shape_rosette, :)     = (/ 1.000E-04_RealK, 3.080E-06_RealK /)
  Alpha(ip_shape_column, :)      = (/ 2.515E-04_RealK, 1.658E-06_RealK /)
!
  Beta(ip_shape_null, :)        = (/ 0.0_RealK, 0.0_RealK /)
  Beta(ip_shape_sphere, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Beta(ip_shape_hexcyl, :)      = (/ 0.0_RealK, 0.0_RealK /)
  Beta(ip_shape_polycrystal, :) = (/ 2.897_RealK, 2.449_RealK /)
  Beta(ip_shape_plate, :)       = (/ 2.852_RealK, 2.449_RealK /)
  Beta(ip_shape_rosette, :)     = (/ 2.997_RealK, 2.260_RealK /)
  Beta(ip_shape_column, :)      = (/ 3.000_RealK, 1.910_RealK /)
!
!
!
  IF (i_shape == ip_shape_sphere) THEN
!
    IF (i_measure == ip_measure_radius) THEN
      volume_particle = (4.0_RealK / 3.0_RealK) * pi * dimen**3
    ELSE IF (i_measure == ip_measure_max_dimen) THEN
      volume_particle = (1.0_RealK / 6.0_RealK) * pi * dimen**3
    ELSE IF (i_measure == ip_measure_proj_radius) THEN
      WRITE(*, *) 'Incomplete Code'
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
      volume_particle = Alpha(i_size_range, i_shape) * &
                          (dimen/Dscale(i_shape))  ** &
                          Beta(i_size_range, i_shape)
    ELSE IF (i_measure == ip_measure_proj_radius) THEN
      WRITE(*, *) 'Incomplete Code'
    ELSE
      WRITE(iu_err, '(/A)') '*** Error: Illegal measure of size.'
      ierr = i_err_fatal
      RETURN
    ENDIF
!
  ELSE IF ( (i_shape == ip_shape_db_undefined) )  THEN
!
    IF (i_measure == ip_measure_radius) THEN
      volume_particle = (4.0_RealK / 3.0_RealK) * pi * dimen**3
    ELSE IF (i_measure == ip_measure_max_dimen) THEN
!     Interpolate in the cube of the maximum dimension.
      CALL spline_evaluate(ierr, DBGeom%n_geom, &
             DBGeom%dm3, DBGeom%volume, DBGeom%d2_volume, &
             dimen**3, volume_particle)
      IF (ierr == i_err_range) THEN
!       Recover harmlessly form errors of interpolation.
        volume_particle = 0.0_RealK
        ierr =  i_normal
      ENDIF
    ELSE IF (i_measure == ip_measure_proj_radius) THEN
      volume_particle = (4.0_RealK / 3.0_RealK) * pi * dimen**3
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

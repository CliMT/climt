! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculation of optical properties after Mitchell.
!
SUBROUTINE AdtMitchell96 &
!
 (WaveLength, ProjArea, Volume, IShape, Refract, &
  Qext, Qabs, Qscat, Asymmetry, ierr)
!
! Description:
!
! This subroutine returns the optical propeties of a particle
! based on Mitchell's scheme of 1996. The inputs are the volume
! and projected area to allow a wider variety of schemes to be
! used.
!
! Modules:
  USE realtype_rd
  USE error_pcf
  USE def_std_io_icf
  USE rad_ccf, ONLY: pi
  USE shape_particle_pcf
!
!
  IMPLICIT NONE
!
!
! Input variables:
!
  INTEGER, INTENT(In) :: IShape
!                            Shape of the crystal
  REAL (RealK), INTENT(In) :: WaveLength
!                               Wavelength of the radiation
  REAL (RealK), INTENT(In) :: ProjArea
!                               Projected area of the particle
  REAL (RealK), INTENT(In) :: Volume
!                               Volume of the particle
!
  COMPLEX (RealK), INTENT(In) :: Refract
!                               Index of refraction of ice
!
! Output variables:
  REAL (RealK), INTENT(Out) :: Qext
!                                Extinction efficiency
  REAL (RealK), INTENT(Out) :: Qabs
!                                Absorption efficiency
  REAL (RealK), INTENT(Out) :: Qscat
!                                Scattering efficiency
  REAL (RealK), INTENT(Out) :: Asymmetry
!                                Asymmetry Factor
  INTEGER, INTENT(InOut) :: ierr
!                           Error flag
!
!
! Local variables:
!
! Variables related to the treatment of ill-conditioning
  REAL (RealK) :: Eps_R
!                    The smallest real number such that 1.0-EPS_R 
!                    is not 1 to the computer's precision
  REAL (RealK) :: Fr_Eps_R
!                    The fourth root of the above
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
  REAL (RealK) :: De
!                   Effective diameter
  REAL (RealK) :: Dmm
!                   Inferred maximum dimension
  REAL (RealK) :: xe
!                   Effective size parameter
  REAL (RealK) :: Qabs_ADT
!                   Absorption efficiency for ADT
  REAL (RealK) :: a1
!                   Correction factor in ADT
  REAL (RealK) :: C1
!                   Correction factor for internal reflection and
!                   refraction
!
  COMPLEX (RealK) :: w
!                      Argument of ADT function
!
!
!- End of Header
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
! Set numerical tolerances:
  Eps_R=EPSILON(Dmm)
  Fr_Eps_R=SQRT(SQRT(Eps_R))
!
!
! Determine the equivalent photon path and the size parameter.
  De=Volume/ProjArea
  xe=Pi*De/WaveLength
!
!
! Extinction:
!
! The standard formula from van de Hulst (p.175) is used.
  w = (0.0_RealK, 2.0_RealK) * (Refract-1.0_RealK) * CMPLX(xe)
!
! The standard formula for Qext becomes ill-conditioned when w is
! small. The error in the standard formula will be roughly
! precision/w^2: we switch when this is comparable to the neglected
! second term in the power series expansion.
  IF (ABS(w) .GT. Fr_Eps_R) THEN 
    Qext = 4.0_RealK * (0.5_RealK + REAL( EXP(-w)/w &
            + (EXP(-w) - (1.0_RealK, 0.0_RealK))/(w * w) ) )
  ELSE
    Qext = 2.0_RealK * REAL(w) / 3.0_RealK
  ENDIF
!
!
! Absorption:
!
! Begin with the absorption efficiency in pure ADT and then add
! corrections.
  Qabs_ADT=1.0_RealK-EXP(-4.0_RealK * AIMAG(Refract) * xe)
!
! Internal Reflection and Refraction:
  IF (IShape == IP_shape_plate) THEN
    a1 = 0.8_RealK
  ELSE IF ( (IShape == IP_shape_polycrystal) .OR. &
            (IShape == IP_shape_rosette) .OR. &
            (IShape == IP_shape_column)  ) THEN
    a1 = 1.0_RealK
  ELSE
    WRITE(IU_ERR, '(/A)') &
      "*** Error: Mitchell's scheme has been selected with an " &
      //"invalid shape."
    ierr=I_err_fatal
    RETURN
  ENDIF
! Modify the correction for smaller size parameters.
  If (xe .LT. 22.0_RealK) &
    a1 = a1 * (1.0_RealK - EXP( 2.0E-09_RealK * xe**7 ) )
  C1=a1 * (1.0_RealK - Qabs_ADT)
!
!
! Include corrections:
  Qabs  = Qabs_ADT * (1.0_RealK + C1)
  Qscat = Qext - Qabs
!
!
! Asymmetry:
!
! The calculation depends on the shape, but really provides values
! only for the ranges 0.2 - 0.7, 0.7 - 1.3, 1.3 - 1.9, 1.9 - 2.7
! and 2.7 - 3.5 microns. We use uniform relationships across each
! range, and extend the outer intervals to cover other frequencies.
! In his original paper, Mitchell uses Da, the median area dimension;
! we therefore employ his power laws using a linear dimenion inferred
! from the projected area with his power laws for that quantity.
  IF ( (IShape == ip_shape_column)  .OR. &
       (IShape == ip_shape_plate)   .OR. &
       (IShape == ip_shape_rosette) .OR. &
       (IShape == ip_shape_polycrystal) ) THEN
!
!   Work out what Dmm is for the small range and change to the other
!   if needed.
    Dmm = ( ProjArea / Sigma(1, IShape) )  ** &
                            ( 1.0_RealK / Delta(1, IShape) )
    IF (Dmm > Dtrans(IShape) ) THEN
      Dmm = ( ProjArea / Sigma(2, IShape) )  ** &
                              ( 1.0_RealK / Delta(2, IShape) )
    ENDIF
!
  ELSE
    WRITE(iu_err, '(/A)') '*** Error: Illegal shape.'
    ierr = i_err_fatal
    RETURN
  ENDIF
!
!
  SELECT CASE (Ishape)
!
    CASE (ip_shape_column)
!
      IF (WaveLength .LT. 7.0E-07_RealK) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9712_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.8986_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 7.0E-07_RealK) .AND. &
                (WaveLength .LT. 1.3E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9780_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.9049_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.3E-06_RealK) .AND. &
                (WaveLength .LT. 1.9E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.007_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.9317_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.9E-06_RealK) .AND. &
                (WaveLength .LT. 2.7E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.046_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.9678_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF (WaveLength .GE. 2.7E-06_RealK) THEN 
        Asymmetry = 0.9595 + 1.076_RealK * (Dmm/Dscale(IShape))
!
      ENDIF
!
    CASE (ip_shape_rosette)
!
      IF (WaveLength .LT. 7.0E-07_RealK) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9712_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.8986_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 7.0E-07_RealK) .AND. &
                (WaveLength .LT. 1.3E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9780_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.9049_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.3E-06_RealK) .AND. &
                (WaveLength .LT. 1.9E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.007_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.9317_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.9E-06_RealK) .AND. &
                (WaveLength .LT. 2.7E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.046_RealK * (Dmm/Dscale(IShape))**0.04473
        ELSE
          Asymmetry =  0.9678_RealK * (Dmm/Dscale(IShape))**0.01942
        ENDIF
!
      ELSE IF (WaveLength .GE. 2.7E-06_RealK) THEN 
        Asymmetry = 0.9595 + 1.076_RealK * (Dmm/Dscale(IShape))
!
      ENDIF
!
!
!
  END SELECT

  
!
!
!
Return
End

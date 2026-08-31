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
 (WaveLength, Dmm, IShape, RhoIce, Refract, &
  Qext, Qabs, ierr)
!
! Description:
!
! This subroutine returns the optical propeties of a particle
! based on Mitchell's scheme of 1996.
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
  REAL (RealK), INTENT(In) :: Dmm
!                               Mean maximum dimension of the crystal
  REAL (RealK), INTENT(In) :: RhoIce
!                               Density of ice in crystals
!
  COMPLEX (RealK), INTENT(In) :: Refract
!                               Index of refraction of ice
!
! Output variables:
  REAL (RealK), INTENT(Out) :: Qext
!                                Extinction efficiency
  REAL (RealK), INTENT(Out) :: Qabs
!                                Absorption efficiency
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
  REAL (RealK) :: Dtrans
!                   Transitional dimension for change in power laws
  REAL (RealK) :: Dscale
!                   Dimension used to scale Dmm
  REAL (RealK) :: Alpha
!                   Coefficient in power law for mass
  REAL (RealK) :: Beta
!                   Exponent in power law for mass
  REAL (RealK) :: Sigma
!                   Coefficient in power law for projected area
  REAL (RealK) :: Delta
!                   Exponent in power law for projected area
!
  REAL (RealK) :: Mass
!                   Mass of the crystal
  REAL (RealK) :: ProjArea
!                   Projected area of the crystal
  REAL (RealK) :: De
!                   Effective diameter
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
!- --------------------------------------------------------------------
!
!
! Set numerical tolerances:
  Eps_R=EPSILON(Dmm)
  Fr_Eps_R=SQRT(SQRT(Eps_R))
!
!
! The mass and projected area are related to the mean maximum dimension
! by power laws; these are copied from Mitchell (1996) and converted to 
! give answers in SI units, but the use of a scaling dimension keeps
! the actual digits the same and prevents the arguments to power laws
! getting out of range. The character of the relationship changes at a 
! transitional size, typically of 100 microns.
!
  IF (IShape .EQ. IP_shape_column) THEN
    Dtrans=1.0E-04_RealK
    Dscale=1.0E-06_RealK
    IF (Dmm .LT. Dtrans) THEN
      Alpha = 2.515E-04_RealK
      Beta  = 3.0_RealK
      Sigma = 6.837E-05_RealK
      Delta = 2.0_RealK
    ELSE
      Alpha = 1.658E-06_RealK
      Beta  = 1.91_RealK
      Sigma = 4.59E-06_RealK
      Delta = 1.415_RealK
    ENDIF
  ELSE IF (IShape .EQ. IP_shape_rosette) THEN
    Dtrans=1.0E-04_RealK
    Dscale=1.0E-06_RealK
    IF (Dmm .LT. Dtrans) THEN
      Alpha = 1.0E-04_RealK
      Beta  = 2.997_RealK
      Sigma = 6.837E-05_RealK
      Delta = 2.0_RealK
    ELSE
      Alpha = 3.08E-06_RealK
      Beta  = 2.26_RealK
      Sigma = 8.687E-06_RealK
      Delta = 1.568_RealK
    ENDIF
  ELSE IF (IShape .EQ. IP_shape_plate) THEN
    Dtrans=1.0E-04_RealK
    Dscale=1.0E-06_RealK
    IF (Dmm .LT. Dtrans) THEN
      Alpha = 4.953E-05_RealK
      Beta  = 2.852_RealK
      Sigma = 2.395E-05_RealK
      Delta = 1.855_RealK
    ELSE
      Alpha = 7.389E-06_RealK
      Beta  = 2.449_RealK
      Sigma = 2.395E-05_RealK
      Delta = 1.855_RealK
    ENDIF
  ELSE IF (IShape .EQ. IP_shape_polycrystal) THEN
    Dtrans=1.0E-04_RealK
    Dscale=1.0E-06_RealK
    IF (Dmm .LT. Dtrans) THEN
      Alpha = 5.801E-04_RealK
      Beta  = 2.879_RealK
      Sigma = 2.285E-05_RealK
      Delta = 1.88_RealK
    ELSE
      Alpha = 7.389E-06_RealK
      Beta  = 2.449_RealK
      Sigma = 2.285E-05_RealK
      Delta = 1.88_RealK
    ENDIF
  ELSE
    WRITE(IU_ERR, '(/A)') &
      "*** Error: Mitchell's scheme has been selected with an " &
      //"invalid shape."
    ierr=I_err_fatal
    RETURN
  ENDIF
  Mass=Alpha*EXP(Beta*LOG(Dmm/Dscale))
  ProjArea=Sigma*EXP(Delta*LOG(Dmm/Dscale))
! Determine the equivalent photon path and the size parameter.
  De=Mass/(RhoIce*ProjArea)
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
  IF (IShape .EQ. IP_shape_plate) THEN
    a1 = 0.8_RealK
  ELSE IF ( (IShape .EQ. IP_shape_polycrystal) .OR. &
            (IShape .EQ. IP_shape_rosette) .OR. &
            (IShape .EQ. IP_shape_column)  ) THEN
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
  Qabs = Qabs_ADT * (1.0_RealK + C1)
!
!
! Asymmetry:
!
! The calculation depends on the shape, but really provides values
! only for the ranges 0.2 - 0.7, 0.7 - 1.3, 1.3 - 1.9, 1.9 - 2.7
! and 2.7 - 3.5 microns. We use uniform relationships across each
! range, and extend the outer intervals to cover other frequencies.
!
  SELECT CASE (Ishape)
!
    CASE (ip_shape_column)
!
      IF (WaveLength .LT. 7.0E-07_RealK) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9712_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.8986_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 7.0E-07_RealK) .AND. &
                (WaveLength .LT. 1.3E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9780_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.9049_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.3E-06_RealK) .AND. &
                (WaveLength .LT. 1.9E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.007_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.9317_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.9E-06_RealK) .AND. &
                (WaveLength .LT. 2.7E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.046_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.9678_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF (WaveLength .GE. 2.7E-06_RealK) THEN 
        Asymmetry = 0.9595 + 1.076_RealK * (Dmm/Dscale)
!
      ENDIF
!
    CASE (ip_shape_rosette)
!
      IF (WaveLength .LT. 7.0E-07_RealK) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9712_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.8986_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 7.0E-07_RealK) .AND. &
                (WaveLength .LT. 1.3E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  0.9780_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.9049_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.3E-06_RealK) .AND. &
                (WaveLength .LT. 1.9E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.007_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.9317_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF ( (WaveLength .GE. 1.9E-06_RealK) .AND. &
                (WaveLength .LT. 2.7E-06_RealK) ) THEN
        IF (Dmm .LT. 5.0E-04_RealK) THEN
          Asymmetry =  1.046_RealK * (Dmm/Dscale)**0.04473
        ELSE
          Asymmetry =  0.9678_RealK * (Dmm/Dscale)**0.01942
        ENDIF
!
      ELSE IF (WaveLength .GE. 2.7E-06_RealK) THEN 
        Asymmetry = 0.9595 + 1.076_RealK * (Dmm/Dscale)
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

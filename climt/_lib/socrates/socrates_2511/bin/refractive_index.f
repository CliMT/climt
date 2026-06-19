! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to return the complex refractive index.
!
! Method:
!   The real and imaginary parts of the refractive index are
!   represented as cubic splines.
!
!- ---------------------------------------------------------------------
      FUNCTION refractive_index(ierr, nd_refract
     &  , wavelength_point, n_refract, wavelength
     &  , re_refract, im_refract, d2_re_refract, d2_im_refract)
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
      USE spline_evaluate_mod, ONLY: spline_evaluate
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
!
!     Sizes of arrays
      INTEGER, Intent(IN) ::
     &    nd_refract
!           Size allocated for wavelengths
!
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n_refract
!           Number of wavelengths
      complex  (RealK) ::
     &    refractive_index
!           Value of refractive index
      REAL  (RealK), Intent(IN) ::
     &    wavelength_point
!           Evaluation point
     &  , wavelength(nd_refract)
!           Array of wavelengths
     &  , re_refract(nd_refract)
!           Real parts of index
     &  , im_refract(nd_refract)
!           Imaginary parts of index
     &  , d2_re_refract(nd_refract)
!           Real second derivative
     &  , d2_im_refract(nd_refract)
!           Imaginary second derivative
!
!     Local variables
      REAL  (RealK) ::
     &    re_refract_index
!           Interpolated real part
     &  , im_refract_index
!           Interpolated imaginary part
!
!
      CALL spline_evaluate(ierr, n_refract, wavelength, re_refract
     &  , d2_re_refract, wavelength_point, re_refract_index)
      IF (ierr /= i_normal) THEN
        IF (ierr == i_err_range) THEN
          WRITE(iu_err, '(/a, 1pe12.5, a, /a)')
     &      '+++ warning: wavelength', wavelength_point , ' m is '
     &      //'outside the range of refractive indices.'
     &      , 'a value of 1.0 will be assumed.'
          re_refract_index=1.0_RealK
          ierr=i_normal
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
!
      CALL spline_evaluate(ierr, n_refract, wavelength, im_refract
     &  , d2_im_refract, wavelength_point, im_refract_index)
      IF (ierr /= i_normal) THEN
        IF (ierr == i_err_range) THEN
          im_refract_index=0.0_RealK
          ierr=i_normal
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
      refractive_index=CMPLX(re_refract_index, im_refract_index)
!
!
!
      RETURN
      END

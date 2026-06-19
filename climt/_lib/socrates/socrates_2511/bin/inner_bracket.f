! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find points just with a band.
!
! Method:
!   A spectrum of wavelengths is passed to the routine. The indices
!   of points of the spectrum just within a given band are found
!   by taking the end-points of the spectrum and moving in.
!
!- ---------------------------------------------------------------------
      SUBROUTINE inner_bracket(ierr
     &  , wave_length_shorter, wave_length_longer
     &  , n_points, spectrum
     &  , i_short, i_long)

      USE realtype_rd
      USE rad_pcf
      USE def_std_io_icf

      IMPLICIT NONE


      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , i_short
!           Short point
     &  , i_long
!           Long point
      INTEGER, Intent(IN) ::
     &    n_points
!           Number of points in spectrum
      REAL  (RealK), Intent(IN) ::
     &    wave_length_shorter
!           Shorter wavelength
     &  , wave_length_longer
!           Longer wavelength
     &  , spectrum(n_points)
!           Spectrum of wavelengths
!
!
!     Check the order of the inputs.
      IF (wave_length_shorter >= wave_length_longer) THEN
        WRITE(iu_err, '(/a)')
     &    '*** error: wavelengths are incorrectly ordered.'
        ierr=i_abort_calculation
        RETURN
      ENDIF
!
!     Initialize I_SHORT and I_LONG and then search.
      i_short=1
      i_long=n_points
1     if (wave_length_shorter.ge.spectrum(i_short)) then
        i_short=i_short+1
        IF (i_short <= n_points) goto 1
      ENDIF
2     if (wave_length_longer.le.spectrum(i_long)) then
        i_long=i_long-1
        IF (i_long >= 1) goto 2
      ENDIF
!
!
!
      RETURN
      END

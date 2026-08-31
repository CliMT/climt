! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 6.
!
! Method:
!   Calculate the thermal coefficients for the LW spectral files.
!   This routine basically just determines whether a filter function
!   is required or not.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_6(Spectrum, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE dimensions_spec_ucf

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local arguments.
  INTEGER :: ios
!   IO status
  CHARACTER (LEN=1) :: l_filter
!   Character flag for filter function


  IF (ALLOCATED(Spectrum%Planck%thermal_coeff)) &
      DEALLOCATE(Spectrum%Planck%thermal_coeff)
  ALLOCATE(Spectrum%Planck%thermal_coeff(0:npd_thermal_coeff-1, &
                                         Spectrum%Dim%nd_band))
  IF (ALLOCATED(Spectrum%Planck%theta_planck_tbl)) &
      DEALLOCATE(Spectrum%Planck%theta_planck_tbl)
  ALLOCATE(Spectrum%Planck%theta_planck_tbl(0:npd_thermal_coeff-1))

  WRITE(*, '(/A)') 'Is a filter function required (Y/N)?'
  DO
    READ(*, *, IOSTAT=ios) l_filter
    IF (ios /= 0) THEN
      WRITE(*, '(A)') '***error: unrecognized response'
      WRITE(*, '(A)') 'Please re-enter.'
    ELSE
      EXIT
    END IF
  END DO

  IF (l_filter=='Y' .OR. l_filter=='y') THEN
    CALL make_block_6_2(ierr, &
      Spectrum%Basic%n_band, &
      Spectrum%Basic%wavelength_short, &
      Spectrum%Basic%wavelength_long, &
      Spectrum%Basic%l_present(14), &
      Spectrum%Basic%n_band_exclude, &
      Spectrum%Basic%index_exclude, &
      Spectrum%Planck%n_deg_fit, &
      Spectrum%Planck%t_ref_planck, &
      Spectrum%Planck%thermal_coeff, &
      Spectrum%Planck%theta_planck_tbl, &
      Spectrum%Basic%l_present(6), &
      Spectrum%Planck%l_planck_tbl )
  ELSE
    CALL make_block_6_1(ierr, &
      Spectrum%Basic%n_band, &
      Spectrum%Basic%wavelength_short, &
      Spectrum%Basic%wavelength_long, &
      Spectrum%Basic%l_present(14), &
      Spectrum%Basic%n_band_exclude, &
      Spectrum%Basic%index_exclude, &
      Spectrum%Planck%n_deg_fit, &
      Spectrum%Planck%t_ref_planck, &
      Spectrum%Planck%thermal_coeff, &
      Spectrum%Planck%theta_planck_tbl, &
      Spectrum%Basic%l_present(6), &
      Spectrum%Planck%l_planck_tbl )
  END IF

END SUBROUTINE make_block_6

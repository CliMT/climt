! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read data for a single wavelenth from the database. 
!
SUBROUTINE db_read_single_wavelength &
!
(iu_db_input, &
 nd_wavelength, nd_scatt_angle, nd_size_scat, &
 wavelength_index, n_wavelength, wavelength, &
 n_rec_block, db_record, &
 l_phase, n_angle, mu_angle, &
 ice_db_mono_info, &
 ierr)

  USE realtype_rd
  USE dimensions_pp_ucf
  USE def_std_io_icf
  USE error_pcf
  USE def_sct_db
  USE def_db_ss_mono
  USE rad_ccf, ONLY: pi
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
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
  INTEGER, Intent(IN) :: iu_db_input
!   Unit number for input from the database
!
  INTEGER, Intent(IN) :: wavelength_index
!   Index of current wavelength in the list of those in the database
  INTEGER, Intent(IN) :: n_wavelength
!   Total number of wavelengths
!
  REAL (RealK), Intent(IN), Dimension(nd_wavelength) :: wavelength
!   Full list of wavelengths
!
  INTEGER, Intent(IN) :: n_angle
!   Number of scattering angles
  LOGICAL, Intent(IN) :: l_phase
!   Calculate phase function
!
! General variables for the database
  INTEGER, Intent(IN), Dimension(nd_wavelength) :: n_rec_block
!   Number of blocks of data for each wavelength
  INTEGER, Intent(IN), Dimension(nd_wavelength, nd_size_scat) :: &
    db_record
!   Pointers to the records for each wavelength and each block at
!   that wavelength
!
  REAL (RealK), Intent(IN), Dimension(nd_scatt_angle) :: mu_angle
!   Points for Gaussian quadrature
!
  TYPE  (STR_db_ss_mono), Intent(OUT) :: ice_db_mono_info
!   Total monochromatic information in the database
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
!
!
! Local variables.
  INTEGER  ::  i
!   Loop variable
  INTEGER  ::  j
!   Loop variable
  INTEGER  ::  k
!   Loop variable
  INTEGER  ::  m
!   Loop variable
  INTEGER  ::  n_r4_per_rec
!   Number of 4-byte reals in a record
  INTEGER  ::  i_data_rec
!   Record of the data file being read
  INTEGER  ::  n_full_rec
!   Number of full records to hold the phase function
  LOGICAL  ::  l_partial_rec
!   Logical to indicate the need for a partial record
!
  TYPE (str_sctdb_entry)  ::  ice_sct
!   Structure holding the ice scattering data
     
  REAL (RealK) :: angle_mu
!   Angle in radians converted from mu
  REAL (RealK) :: gauss_left_off
!   Left-hand weighting for interpolation
  REAL (RealK) :: gauss_right_off
!   Right-hand weighting for interpolation
!
  INTEGER  :: phf_index
!   Index of point of phase function in array
!
  REAL (RealK), Allocatable, Dimension(:) :: phft
!   Windowed phase function

  INTEGER, Parameter :: nf=5
!   Width of window
  REAL (RealK), Dimension(-nf:+nf) :: wt
!   Windowing weights

  CHARACTER  (LEN=8) :: date
  CHARACTER  (LEN=10):: time
!
!
!
! Check if wavelength index is valid
  IF (wavelength_index < 1 .OR. wavelength_index > n_wavelength) THEN
    WRITE(iu_err, '(/A)') '*** Error: There are ' // &
      'not as many wavelengths in the database file as expected.'
    ierr = i_err_fatal
    RETURN
  ENDIF
!
!
! Allcoate space for all information at this frequency.
  ALLOCATE(ice_db_mono_info%dm(n_rec_block(wavelength_index)))
  ALLOCATE(ice_db_mono_info%ss(n_rec_block(wavelength_index), &
    n_angle+4))
  ALLOCATE(ice_db_mono_info%d2_ss(n_rec_block(wavelength_index), &
    n_angle+4))

! Set generic sizes.
  n_r4_per_rec=sct_db_recl/8
!
! Check the size of the record
  IF (n_r4_per_rec < 10) THEN
    WRITE(iu_err, '(/A)') &
      '*** Error: The predefined record length is too small.'
    ierr = i_err_fatal
    RETURN
  ENDIF
!
!
! Read data from those records of the database,
! which refer to the wavelength chosen
  DO j=1, n_rec_block(wavelength_index)

!   Read the header
    READ(iu_db_input, REC=db_record(wavelength_index,j)) &
      ice_sct%dm, ice_sct%wavelength, ice_sct%csca, ice_sct%cext, &
      ice_sct%omega, ice_sct%asymm, ice_sct%n_angle, &
      ice_sct%l_uniform, ice_sct%phf_angle(1), ice_sct%d_angle, &
      i_data_rec

!   Check if wavelength in database is as expected
    IF (ABS(wavelength(wavelength_index)-ice_sct%wavelength) > & 
          1.0E-5*MIN(wavelength(wavelength_index), &
          REAL(ice_sct%wavelength, RealK)) ) THEN
      WRITE(iu_err, '(/A)') '*** Error: The wavelength requested' // &
        ' was not found in the database file.'
      ierr = i_err_fatal
      RETURN
    ENDIF
!
    IF (l_phase) THEN
!
!     Now read the phase function. This may be given on a uniform
!     or non-uniform set of angles.
      n_full_rec=ice_sct%n_angle/n_r4_per_rec
      l_partial_rec=((ice_sct%n_angle - n_r4_per_rec*n_full_rec) > 0)
!
      IF (.NOT.ice_sct%l_uniform) THEN
!       We need to read the angles as well.
        DO i=1, n_full_rec
           READ(iu_db_input, REC=i_data_rec) &
             (ice_sct%phf_angle(k), k=1+(i-1)*n_r4_per_rec, i*n_r4_per_rec)
             i_data_rec=i_data_rec+1
        ENDDO
        IF (l_partial_rec) THEN
           READ(iu_db_input, REC=i_data_rec) &
             (ice_sct%phf_angle(k), &
              k=1+n_full_rec*n_r4_per_rec, ice_sct%n_angle)
             i_data_rec=i_data_rec+1
        ENDIF
      ELSE
!       Set the angels for use shortly.
        DO i=1, ice_sct%n_angle
          ice_sct%phf_angle(i) = ice_sct%d_angle * REAL((i-1), RealK)
        ENDDO
      ENDIF
!
!     Read the phase function.
      DO i=1, n_full_rec
        READ(iu_db_input, REC=i_data_rec) &
          (ice_sct%phf(k), k=1+(i-1)*n_r4_per_rec, i*n_r4_per_rec)
          i_data_rec=i_data_rec+1
      ENDDO
      IF (l_partial_rec) THEN
      READ(iu_db_input, REC=i_data_rec) &
        (ice_sct%phf(k), &
        k=1+n_full_rec*n_r4_per_rec, ice_sct%n_angle)
        i_data_rec=i_data_rec+1
      ENDIF

    ENDIF
!
!
!
!   Discontinuities in the information in the database cause
!   problems in calculating moments. Applying a smoothing window
!   is an interim solution and a triangular window is selected.
!
    DO m=-nf,+nf
      wt(m)=REAL(nf-ABS(m), RealK)/REAL(nf**2, RealK)
    ENDDO
!
    ALLOCATE(phft(ice_sct%n_angle))
    phft=0.0_RealK
!
    DO i=1, ice_sct%n_angle
      DO m=-nf,+nf
        IF ((i+m) < 1) THEN
          phft(i)=phft(i)+ice_sct%phf(1)*wt(m) 
        ELSE IF ((i+m) > ice_sct%n_angle) THEN
          phft(i)=phft(i)+ice_sct%phf(ice_sct%n_angle)*wt(m)
        ELSE
          phft(i)=phft(i)+ice_sct%phf(i+m)*wt(m)
        ENDIF
      ENDDO
    ENDDO
!
!
!
!   Copy across the information in the current database.
!
    ice_db_mono_info%dm(j)    = ice_sct%dm
    ice_db_mono_info%ss(j, 1) = ice_sct%csca
    ice_db_mono_info%ss(j, 2) = ice_sct%cext
    ice_db_mono_info%ss(j, 3) = ice_sct%asymm
    ice_db_mono_info%ss(j, 4) = ice_sct%omega
!
!
!   Interpolate to the required points. We assume that angle increases
!   monotonically.
    IF (l_phase) THEN
      phf_index = ice_sct%n_angle
      DO i=1, n_angle
!       The database is in raw angles, increasing with entry, while
!       the list of Guassian roots increases in cosine, so the 
!       corresponding angles decrease.
        angle_mu = ACOS(mu_angle(i))
        DO
          IF (ice_sct%phf_angle(phf_index) <= angle_mu)  EXIT
          phf_index = phf_index - 1
          IF (phf_index == 1) EXIT
        ENDDO
!       Ensure we have a valid interval.
        phf_index = MAX(phf_index, 1)
        phf_index = MIN(phf_index, ice_sct%n_angle-1)
!       phf_index is now to the left of the point we want.
        gauss_left_off = &
          (ice_sct%phf_angle(phf_index+1) - angle_mu) / &
          (ice_sct%phf_angle(phf_index+1) - &
           ice_sct%phf_angle(phf_index))
        gauss_right_off=1.0_RealK-gauss_left_off
        ice_db_mono_info%ss(j, 4+i) = &
          phft(phf_index)*gauss_left_off + &
          phft(phf_index+1)*gauss_right_off
      ENDDO 
    ENDIF
!
    DEALLOCATE(phft)
!
  ENDDO
!
!
! Set the second derivatives for spline fitting. At present, because
! of discontinuities in the database we zero the second derivatives,
! so the fits will effectively be linear.
  ice_db_mono_info%d2_ss = 0.0_RealK
!
!
!
  RETURN
END SUBROUTINE db_read_single_wavelength

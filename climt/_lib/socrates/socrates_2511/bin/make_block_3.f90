! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 3.
!
! Method:
!   A solar spectrum is read if necessary. The monochromatic
!   Rayleigh scattering coefficients are calculated and 
!   weighted with the solar spectrum (in make_block_3_1).
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_3(Sp, Sol, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE gas_list_pcf, ONLY: npd_gases, name_absorb
  USE dimensions_pp_ucf
  USE def_std_io_icf
  USE def_solarspec
  USE def_refract, ONLY: StrRefract, allocate_refract, deallocate_refract

  IMPLICIT NONE


  TYPE (StrSpecData), INTENT(INOUT) :: Sp
!   Spectral file to be assigned
  TYPE (StrSolarSpec), INTENT(INOUT) :: Sol
!   Solar spectrum
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables
  INTEGER :: i
!   Loop variable
  TYPE (StrRefract) :: Refract_H2
!   Refractive index of H2
  INTEGER :: iu_n_h2
!   File unit for H2 refractive index data
  INTEGER :: ios
!   I/O error flag
  CHARACTER(LEN=1) :: char_atm
!   Response to atmospheric composition question
  CHARACTER(LEN=1) :: char_yn
!   Character with yes/no response
 LOGICAL :: l_h2he_atm = .FALSE.
!   Calculate rayleigh scattering for H2/He atmosphere
  CHARACTER(LEN=80) :: line
!   Line with response from user
  INTEGER :: type_gas
!   Actual type of gas
  LOGICAL :: l_retain_gases
!   Flag for using the gases for which Rayleigh scattering coefficients are
!   already tabulated in the spectral file
  INTEGER :: type_index(npd_gases) = (/ (-1, i = 1, npd_gases) /)
!   The spectral indices of each gas in the file


! If a solar spectrum is already present that is used, otherwise
! one is read in.
  IF (Sol%n_points > 0) THEN
    WRITE(*, '(/a)') &
      'Rayleigh scattering coefficients will be averaged using '
    WRITE(*, '(a/)') &
      'the solar spectrum read in earlier.'
  ELSE
    CALL read_solar_spectrum(Sol, ierr)
    IF (ierr /= i_normal) RETURN
  ENDIF

! Is the atmosphere composed of air, H2-He gas or a custom composition?
  WRITE(*, '(/A)') 'Is the atmosphere composed of air, H2-He gas or ' // &
      'a custom composition? (A/H/C)'
  DO
    READ(*, *, IOSTAT=ios) char_atm
    IF (ios /= 0) THEN
      WRITE(*, '(A)') '***error: unrecognized response'
      WRITE(*, '(A)') 'Please re-enter.'
    ELSE
      EXIT
    ENDIF
  ENDDO

  IF (char_atm=='H' .OR. char_atm=='h') THEN
!   Rayleigh scattering by H2/He, read in the refractive index for H2

    Sp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_total
    IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
        DEALLOCATE(Sp%Rayleigh%rayleigh_coeff)
    ALLOCATE(Sp%Rayleigh%rayleigh_coeff(Sp%Dim%nd_band))

!   Obtain the file containing the H2 refractive index data.
    CALL get_free_unit(ierr, iu_n_h2)
    CALL open_file_in(ierr, iu_n_h2, & 
      'Enter the name of the file containing the H2 refractive index data.')
    IF (ierr /= i_normal) RETURN

!   Read first to find the number of points in the spectrum.
    Refract_H2%n_points = 0
    ios = 0
    DO WHILE (ios == 0)
      READ(iu_n_h2,'(f22.15)',IOSTAT=ios)
      Refract_H2%n_points = Refract_H2%n_points + 1
    END DO
    Refract_H2%n_points = Refract_H2%n_points - 2

    CALL allocate_refract(Refract_H2)

!   Read in the file.
    REWIND(iu_n_h2)
    READ(iu_n_h2, '(A)', IOSTAT=ios) ! Skip header line
    DO i = 1, Refract_H2%n_points
      READ(iu_n_h2, *, IOSTAT=ios) &
        Refract_H2%wavelength(i), &
        Refract_H2%re_part(i)
      Refract_H2%wavelength(i) = &
        Refract_H2%wavelength(i)*1e-6_RealK ! Convert to metre
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/A)') '*** Error: Corrupt refractive index data.'
        ierr=i_err_fatal
        RETURN
      ENDIF
    ENDDO

!   Check that reading was carried out: failure to read the data
!   will cause n_points to be 0.
    IF (Refract_H2%n_points == 0) THEN
      WRITE(iu_err, '(/a)') &
        '*** Error: No data were read. ' // &
        'Check format of file of refractive index data.'
      ierr=i_err_fatal
      RETURN
    ENDIF

    l_h2he_atm = .TRUE.

  ELSE IF (char_atm == 'C' .OR. char_atm == 'c') THEN
!   Rayleigh scattering coefficients will be tabulated independently
!   for each gas

!   Find the spectral indices of each gas in the file
    DO i = 1, Sp%Gas%n_absorb
      type_index(Sp%Gas%type_absorb(i)) = i
    END DO

    IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff_gas)) &
        DEALLOCATE(Sp%Rayleigh%rayleigh_coeff_gas)
    ALLOCATE(Sp%Rayleigh%rayleigh_coeff_gas(Sp%Dim%nd_species,Sp%Dim%nd_band))

!   If block is already present in spectral file, and the Rayleigh scattering
!   scheme is custom, ask if Rayleigh scattering coefficients should be
!   tabulated for the same set of gases
    IF (Sp%Basic%l_present(3) .AND. &
        Sp%Rayleigh%i_rayleigh_scheme == ip_rayleigh_custom) THEN

      WRITE(*, '(/A/A)') 'Presently Rayleigh scattering coefficients are', &
          'tabulated for the following gases:'
      DO i = 1, Sp%Rayleigh%n_gas_rayleigh
        type_gas = Sp%Gas%type_absorb(Sp%Rayleigh%index_rayleigh(i))
        WRITE(*, '(A)') TRIM(name_absorb(type_gas))
      END DO
      WRITE(*, '(/A/A)') 'Are these the gases for which you would like to', &
          'tabulate Rayleigh scattering coefficients? (Y/N)'
      DO
        READ(*, *, IOSTAT=ios) char_yn
        IF (ios /= 0 .OR. .NOT. (char_yn == 'Y' .OR. char_yn == 'y' .OR. &
            char_yn == 'N' .OR. char_yn == 'n')) THEN
          WRITE(iu_err, '(A)') '***Error: unrecognized response'
          WRITE(iu_err, '(A)') 'Please re-enter.'
        ELSE IF (char_yn == 'Y' .OR. char_yn == 'y') THEN
          l_retain_gases = .TRUE.
          EXIT
        ELSE
          l_retain_gases = .FALSE.
          EXIT
        ENDIF
      ENDDO

    ELSE
      l_retain_gases = .FALSE.
      Sp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_custom
    END IF

    IF (.NOT. l_retain_gases) THEN
!     Obtain the number of gases for which the Rayleigh scattering should be
!     tabulated
      WRITE(*, '(/A)') 'Enter the number of gases for which Rayleigh ' // &
          'scattering coefficients should be tabulated.'
      WRITE(*, '(A/)') 'Enter A to include all gases listed in block 0 ' // &
          'in spectral file.'
      DO
        READ(*, *, IOSTAT=ios) line
        IF (ios /= 0) THEN
          WRITE(iu_err, '(A)') '***Error: unrecognized response'
          WRITE(iu_err, '(A)') 'Please re-enter.'
        ELSE
!         Check for valid response
          IF (TRIM(line) == 'A' .OR. TRIM(line) == 'a') THEN
            Sp%Rayleigh%n_gas_rayleigh = Sp%Gas%n_absorb
            EXIT
          ELSE
            READ(line, *, IOSTAT=ios) Sp%Rayleigh%n_gas_rayleigh
            IF (ios /= 0) THEN
              WRITE(iu_err, '(A)') '***Error: unrecognized response'
              WRITE(iu_err, '(A)') 'Please re-enter.'
            ELSE IF (Sp%Rayleigh%n_gas_rayleigh > Sp%Gas%n_absorb) THEN
              WRITE(iu_err, '(A)') '***Error: Number cannot be greater ' // &
                  'than number of absorbers in spectral file.'
              WRITE(iu_err, '(A)') 'Please re-enter.'
            ELSE IF (Sp%Rayleigh%n_gas_rayleigh <= 0) THEN
              WRITE(iu_err, '(A)') '***Error: Rayleigh scattering ' // &
                  'coefficients must be included for at least one gas.'
              WRITE(iu_err, '(A)') 'Please re-enter.'
            ELSE
              EXIT
            END IF
          END IF
        ENDIF
      ENDDO

      IF (ALLOCATED(Sp%Rayleigh%index_rayleigh)) &
          DEALLOCATE(Sp%Rayleigh%index_rayleigh)
      ALLOCATE(Sp%Rayleigh%index_rayleigh(Sp%Dim%nd_species))
      IF (Sp%Rayleigh%n_gas_rayleigh == Sp%Gas%n_absorb) THEN
        WRITE(*, '(A)') 'Including Rayleigh scattering coefficients for ' // &
            'all gases listed in block 0 in spectral file.'
        DO i = 1, Sp%Gas%n_absorb
          Sp%Rayleigh%index_rayleigh(i) = i
        END DO
      ELSE
!       Get list of gases for which to include Rayleigh scattering coefficients
        DO i = 1, Sp%Rayleigh%n_gas_rayleigh
          WRITE(*, '(A, 1X, I2)') &
              'Enter the identifier for Rayleigh gas', i
          DO
            READ(iu_stdin, *, IOSTAT=ios) type_gas
            IF (ios == 0 .AND. type_gas > 0 .AND. type_gas <= npd_gases) THEN
              IF (type_index(type_gas) < 1) THEN
                WRITE(iu_err, '(A)') '***Error: gas not in spectral file'
                WRITE(iu_err, '(A)') 'Please re-enter.'
              ELSE
                Sp%Rayleigh%index_rayleigh(i) = type_index(type_gas)
                EXIT
              END IF
            ELSE
              WRITE(iu_err, '(A)') '***Error: unrecognized response'
              WRITE(iu_err, '(A)') 'Please re-enter.'
            END IF
          ENDDO
        ENDDO
      END IF
    END IF

  ELSE
!   Rayleigh scattering by air
    Sp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_total
    IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
        DEALLOCATE(Sp%Rayleigh%rayleigh_coeff)
    ALLOCATE(Sp%Rayleigh%rayleigh_coeff(Sp%Dim%nd_band))

  END IF

  CALL make_block_3_1(Sp, Sol, Refract_H2, l_h2he_atm)
  IF (l_h2he_atm) CALL deallocate_refract(Refract_H2)

  ! Treat single custom gas as total for compatibility with sub-bands
  IF (Sp%Rayleigh%i_rayleigh_scheme == ip_rayleigh_custom .AND. &
      Sp%Rayleigh%n_gas_rayleigh == 1) THEN
    WRITE(*, '(A)') &
      'Assume entire atmosphere is this gas for Rayleigh scattering? (Y/N)'
    READ(*, *, IOSTAT=ios) char_yn
    IF (char_yn == 'Y' .OR. char_yn == 'y') THEN
      Sp%Rayleigh%i_rayleigh_scheme = ip_rayleigh_total
      IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
          DEALLOCATE(Sp%Rayleigh%rayleigh_coeff)
      ALLOCATE(Sp%Rayleigh%rayleigh_coeff(Sp%Dim%nd_band))
      DO i=1, Sp%Basic%n_band
        Sp%Rayleigh%rayleigh_coeff(i) = Sp%Rayleigh%rayleigh_coeff_gas(1,i)
      END DO
      DEALLOCATE(Sp%Rayleigh%rayleigh_coeff_gas)
      DEALLOCATE(Sp%Rayleigh%index_rayleigh)
    END IF
  END IF

  Sp%Basic%l_present(3)=.TRUE.

END SUBROUTINE make_block_3

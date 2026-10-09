! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to construct a spectral file.
!
! Description:
!   This program is used to interactively contruct a spectral file,
!   either entirely from scratch or based on an existing file.
!
! Method:
!   For a new file the user is asked to supply the limits on the
!   spectral bands and the gases and aerosols present. Alternatively,
!   a spectral file is read in. The user then selects a block of data
!   to add to the file. The resulting file is then written out.
!
!- ---------------------------------------------------------------------
PROGRAM prep_spec

  USE realtype_rd
  USE def_spectrum
  USE def_solarspec
  USE rad_pcf
  USE gas_list_pcf
  USE dimensions_pp_ucf
  USE dimensions_spec_ucf

  IMPLICIT NONE


  INTEGER :: ierr = 0
!   Error flag
  LOGICAL :: l_interactive
!   Flag for interactive operation
  LOGICAL :: l_exclude
!   Flag for the presence of non-contiguous bands
  LOGICAL :: l_exist
!   Existence flag for file
  INTEGER :: i, i_block
!   Loop variable
  CHARACTER (LEN=132) :: file_spectral
!   Name of the spectral file
  CHARACTER (LEN=1) :: char_in
!   Character response variable
  INTEGER, Dimension(npd_gases) :: type_index = &
                                (/ (-1, i = 1, npd_gases) /)
!   The spectral indices of each gas in the file

  TYPE(StrSpecData) :: Spectrum
!   The spectral configuration to be defined
  TYPE(StrSolarSpec) :: SolarSpec
!   The solar spectral irradiance data

! Functions called:
  LOGICAL, EXTERNAL :: set_interactive
!   Function to determine whether operation is interactive


! Set the flag for interactive operation
  l_interactive=set_interactive()

! Read in the spectral file.
  WRITE(*, "(a)") "Enter the name of the spectral file."
  READ(*, "(a)") file_spectral
  INQUIRE(FILE=file_spectral, EXIST=l_exist)
  IF (.NOT.l_exist) THEN

    WRITE(*, "(a)") "A new spectral file will be created."

!   Obtain the summary information.
    CALL make_block_0(Spectrum, type_index, l_interactive, ierr)
    IF (ierr /= i_normal) STOP
    
!   Obtain the limits on the wavelengths.
    CALL make_block_1(Spectrum, l_interactive, ierr)
    IF (ierr /= i_normal) STOP
    
!   Set the types of absorber in each band
    CALL make_block_4(Spectrum, type_index, l_interactive, ierr)
    IF (ierr /= i_normal) STOP
    
!   Set the continuum types.
    CALL make_block_8(Spectrum, l_interactive, ierr)
    IF (ierr /= i_normal) STOP
    
!   Exclude regions from bands if required.
    CALL make_block_14(Spectrum, l_exclude, l_interactive, ierr)
    IF (ierr /= i_normal) STOP

!   Set the types of continua in each band
    IF (Spectrum%ContGen%n_cont > 0) THEN
      CALL make_block_18(Spectrum, l_interactive, ierr)
      IF (ierr /= i_normal) STOP
    END IF

!   Set the index number of water vapour for use with the continuum.
    Spectrum%Cont%index_water = 0
    DO i=1, Spectrum%Gas%n_absorb
      IF (Spectrum%Gas%type_absorb(i) == IP_h2o) THEN
        Spectrum%Cont%index_water = i
      END IF
    END DO
    
!   Setting up of the basic spectrum is now complete.
    Spectrum%Dim%nd_type = npd_type
    Spectrum%Dim%nd_sub_band_gas = 1
    ALLOCATE(Spectrum%Basic%l_present(0:Spectrum%Dim%nd_type))
    Spectrum%Basic%l_present(:)   = .FALSE.
    Spectrum%Basic%l_present(0:1) = .TRUE.
    Spectrum%Basic%l_present(4)   = .TRUE.
    Spectrum%Basic%l_present(8)   = .TRUE.
    Spectrum%Basic%l_present(14)  = l_exclude
    Spectrum%Basic%l_present(18)  = Spectrum%ContGen%n_cont > 0
    Spectrum%Planck%l_planck_tbl  = .FALSE.

  ELSE

    CALL read_spectrum(file_spectral, Spectrum)
    WRITE(*, '(a)') 'Type "a" to append data to the existing file;'
    WRITE(*, '(a)') '  or "n" to create a new file.'
    READ(*, '(a)') char_in
    IF (char_in /= 'A' .AND. char_in /= 'a') THEN
      WRITE(*, '(a)') 'Enter the name of the new file.'
      DO
        READ(*, '(a)') file_spectral
        INQUIRE(FILE=file_spectral, EXIST=l_exist)
        IF (l_exist) THEN
          WRITE(*, '(a)') 'This file already exists: '// &
            'do you wish to overwrite? (y/n)'
          READ(*, '(a)') char_in
          IF (char_in == 'Y' .OR. char_in == 'y') EXIT
          WRITE(*, '(a)') 'Please specify another file.'
        ELSE
          EXIT
        END IF
      END DO
    END IF

  END IF

  DO
!   Now decide which blocks are to be written.
    WRITE(*, '(/a/,15(6x, a/),/)')                                      &
      'Select from the following types of data:',                       &
      '0.   Block 0: Change number of bands.',                          &
      '2.   Block 2: Solar spectrum in each band.',                     &
      '3.   Block 3: Rayleigh scattering in each band.',                &
      '5.   Block 5: k-terms and p, T scaling data.',                   &
      '6.   Block 6: Thermal source function in each band.',            &
      '9.   Block 9: Continuum extinction and scaling data.',           &
      '10.  Block 10: Droplet parameters in each band.',                &
      '11.  Block 11: Aerosol parameters in each band.',                &
      '12.  Block 12: Ice crystal parameters in each band.',            &
      '15.  Block 15: Parameters for aerosol optical depths.',          &
      '17.  Block 17: Spectral variability data in sub-bands.',         &
      '19.  Block 19: Continuum k-terms and T scaling data.',           &
      '20.  Block 20: Photolysis pathways.',                            &
      '-1.  To write spectral file and exit.',                          &
      '-2.  To quit without writing spectral file.'
    READ(*, *) i_block

!   Check whether the block is present and request confirmation
!   for overwriting.
    IF ( (i_block /= -1).AND.(i_block /= -2).AND.(i_block /= 10).AND.   &
         (i_block /= 11).AND.(i_block /= 12) ) THEN
      IF (Spectrum%Basic%l_present(i_block)) THEN
        WRITE(*, '(/a)') 'This block already exists'
        WRITE(*, '(a/)') 'Continue ? (y/n)'
        READ(*, '(a)') char_in
          IF ( (char_in /= 'Y').AND.(char_in /= 'y') ) CYCLE
      END IF
    END IF

!   For each valid type of block call the appropriate routine.
    SELECT CASE (i_block)
    CASE (-2)
      EXIT
    CASE (-1)
!     Write out the spectral file.
      CALL out_spectrum(file_spectral, Spectrum, ierr)
      EXIT
    CASE(0)
      CALL change_block_0(Spectrum)
    CASE(2)
      CALL make_block_2(Spectrum, SolarSpec, ierr)
    CASE(3)
      CALL make_block_3(Spectrum, SolarSpec, ierr)
    CASE(5)
      CALL make_block_5(Spectrum, ierr)
    CASE(6)
      CALL make_block_6(Spectrum, ierr)
    CASE(9)
      CALL make_block_9(Spectrum, ierr)
    CASE(10)
      CALL make_block_10(Spectrum, ierr)
    CASE(11)
      CALL make_block_11(Spectrum, ierr)
    CASE(12)
      CALL make_block_12(Spectrum, ierr)
    CASE(15)
      CALL make_block_15(Spectrum, ierr)
    CASE(17)
      CALL make_block_17(Spectrum, SolarSpec, ierr)
    CASE(19)
      CALL make_block_19(Spectrum, ierr)
    CASE(20)
      CALL make_block_20(Spectrum, ierr)
    CASE DEFAULT
      WRITE(*, '(a)') '+++ Invalid block number.'
    END SELECT
  END DO

END PROGRAM prep_spec

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 0.
!
SUBROUTINE make_block_0 &
!
(Spectrum, type_index, l_interactive, ierr)
!
! Description:
!   This routine constructs the basic summmary information for
!   a spectral file.
!
! Method:
!      Straightforward.
!
! Modules used
  USE realtype_rd
  USE def_spectrum
  USE def_std_io_icf
  USE gas_list_pcf
  USE rad_pcf
  USE dimensions_spec_ucf
  USE dimensions_pp_ucf
!
!
  IMPLICIT NONE
!
!
! Arguments:
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  TYPE (StrSpecData), Intent(OUT) :: Spectrum
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(INOUT) :: type_index(npd_gases)
!   Indices of types of absorbers
!
!
! Local variables:
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i
!   Loop variable
  INTEGER :: i_type_1, i_type_2
!   Types of continuum gases
!
!- End of Header


  WRITE(iu_stdout, '(/A/)') 'Enter number of spectral bands.' 
  DO
    READ(iu_stdin, *, IOSTAT=ios) Spectrum%Basic%n_band
    IF ( (ios == 0) .AND. &
         (Spectrum%Basic%n_band > 0) ) EXIT
    IF (l_interactive) THEN
      WRITE(iu_err, '(A)') &
        '+++ Invalid input: Please re-enter.'
    ELSE
      WRITE(iu_err, '(A)') &
        '*** Error: Invalid number of bands.'
      ierr = i_err_fatal
      RETURN
    ENDIF
  ENDDO
  Spectrum%Dim%nd_band = Spectrum%Basic%n_band


  WRITE(iu_stdout, '(/A/)') 'Enter number of absorbing gases.' 
  DO
    READ(iu_stdin, *, IOSTAT=ios) Spectrum%Gas%n_absorb
    IF ( (ios == 0) .AND. &
         (Spectrum%Gas%n_absorb >= 0) .AND. &
         (Spectrum%Gas%n_absorb <= npd_gases) ) EXIT
    IF (l_interactive) THEN
      WRITE(iu_err, '(A)') &
        '+++ Invalid input: Please re-enter.'
    ELSE
      WRITE(iu_err, '(A)') &
        '*** Error: Invalid number of gaseous absorbers.'
      ierr = i_err_fatal
      RETURN
    ENDIF
  ENDDO
  Spectrum%Dim%nd_species = MAX(Spectrum%Gas%n_absorb, 1)

! Prepare the list of gaeous absorbers.
  ALLOCATE(Spectrum%Gas%type_absorb(Spectrum%Dim%nd_species))
  IF (Spectrum%Gas%n_absorb > 0) THEN
    WRITE(iu_stdout, '(/A)') &
      'Enter the physical types of absorber.'
    DO i = 1, Spectrum%Gas%n_absorb
      WRITE(iu_stdout, '(A29, 1X, I5)') &
          'Enter the identifier for gas ',i
      DO
        READ(iu_stdin, *, IOSTAT=ios) Spectrum%Gas%type_absorb(i)
        IF ( (ios == 0)         .AND. &
           (Spectrum%Gas%type_absorb(i) > 0) .AND. &
           (Spectrum%Gas%type_absorb(i) <= npd_gases) ) THEN
!         Set the index for this type for future use.
          type_index(Spectrum%Gas%type_absorb(i)) = i
          EXIT
        ELSE IF (l_interactive) THEN
          WRITE(iu_err, '(A)') &
            '+++ Invalid input: Please re-enter.'
        ELSE
          WRITE(iu_err, '(A)') &
            '*** Error: Invalid gaseous identifier.'
          ierr = i_err_fatal
          RETURN
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  WRITE(*, '(/A/)') &
    'Enter the number of generalised continua to be included.'
  DO
    READ(iu_stdin, *, IOSTAT=ios) Spectrum%ContGen%n_cont
    IF (ios == 0 .AND. Spectrum%ContGen%n_cont >= 0) THEN
      EXIT
    ELSE IF (l_interactive) THEN
      WRITE(iu_err, '(A)') &
        '+++ Invalid input: Please re-enter.'
    ELSE
      WRITE(iu_err, '(A)') &
        '*** Error: Invalid number of continua.'
      ierr = i_err_fatal
      RETURN
    END IF
  END DO
  Spectrum%Dim%nd_cont = Spectrum%ContGen%n_cont

! Prepare list of continuum gas pairs.
  ALLOCATE(Spectrum%ContGen%index_cont_gas_1(Spectrum%ContGen%n_cont))
  ALLOCATE(Spectrum%ContGen%index_cont_gas_2(Spectrum%ContGen%n_cont))
  IF (Spectrum%ContGen%n_cont > 0) THEN
    WRITE(iu_stdout, '(/A)') &
      'Enter the physical types of the continuum gases.'
    DO i = 1, Spectrum%ContGen%n_cont
      WRITE(*,'(A,I5)') &
        'Enter the gas identifiers for continuum ', i
      DO
        READ(iu_stdin, *, IOSTAT=ios) &
          i_type_1, &
          i_type_2
        IF (ios == 0 .AND. &
            ANY(i_type_1 == &
                Spectrum%Gas%type_absorb(1:Spectrum%Gas%n_absorb)) .AND. &
            ANY(i_type_2 == &
                Spectrum%Gas%type_absorb(1:Spectrum%Gas%n_absorb))) THEN
          Spectrum%ContGen%index_cont_gas_1(i) = type_index(i_type_1)
          Spectrum%ContGen%index_cont_gas_2(i) = type_index(i_type_2)
          EXIT
        ELSE IF (l_interactive) THEN
            WRITE(iu_err, '(A)') &
              '+++ Invalid input: Please re-enter.'
        ELSE
          WRITE(iu_err, '(A)') &
            '*** Error: Invalid gas identifiers.'
          ierr = i_err_fatal
          RETURN
        END IF
      END DO
    END DO
  END IF

  WRITE(iu_stdout, '(/A/)') 'Enter number of aerosols.' 
  DO
    READ(iu_stdin, *, IOSTAT=ios) Spectrum%Aerosol%n_aerosol
    IF ( (ios == 0)       .AND. &
         (Spectrum%Aerosol%n_aerosol >= 0) .AND. &
         (Spectrum%Aerosol%n_aerosol <= npd_aerosol_component) ) EXIT
    IF (l_interactive) THEN
      WRITE(iu_err, '(A)') &
        '+++ Invalid input: Please re-enter.'
    ELSE
      WRITE(iu_err, '(A)') &
        '*** Error: Invalid number of aerosols.'
      ierr = i_err_fatal
      RETURN
    ENDIF
  ENDDO

! Prepare the list of aerosols.
  Spectrum%Dim%nd_aerosol_species = Spectrum%Aerosol%n_aerosol
  ALLOCATE(Spectrum%Aerosol%type_aerosol(Spectrum%Dim%nd_aerosol_species))
  IF (Spectrum%Aerosol%n_aerosol > 0) THEN
    WRITE(iu_stdout, '(/A)') &
      'Enter the physical types of aerosol.'
    DO i = 1, Spectrum%Aerosol%n_aerosol
      WRITE(iu_stdout, '(A33, 1X, I5)') &
        'Enter the identifier for aerosol ',i
      DO
        READ(iu_stdin, *, IOSTAT=ios) Spectrum%Aerosol%type_aerosol(i)
        IF ( (ios == 0)         .AND. &
           (Spectrum%Aerosol%type_aerosol(i) > 0) .AND. &
           (Spectrum%Aerosol%type_aerosol(i) <= npd_aerosol_component) ) THEN
          EXIT
        ELSE IF (l_interactive) THEN
          WRITE(iu_err, '(A)') &
            '+++ Invalid input: Please re-enter.'
        ELSE
          WRITE(iu_err, '(A)') &
            '*** Error: Invalid aerosol identifier.'
          ierr = i_err_fatal
          RETURN
        ENDIF
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE make_block_0

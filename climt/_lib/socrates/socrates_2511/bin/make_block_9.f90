! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 9.
!- ---------------------------------------------------------------------
SUBROUTINE make_block_9(Spectrum, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE file_type_pcf

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local arguments.
  CHARACTER (LEN=80) :: line
!   Input line
  INTEGER :: iu_esft
!   Unit number for file of ESFT/k-term data
  INTEGER :: ios
!   IO status
  INTEGER :: i_input_type
!   Type of input file
  INTEGER :: i_index
!   Index number of absorber
  INTEGER :: i_band
!   Number of band
  INTEGER :: i, j, l
!   Loop variables
  LOGICAL :: l_index_band(Spectrum%Dim%nd_band, Spectrum%Dim%nd_continuum)
!   Continua present
  REAL (RealK) :: dummy
!   Dummy reading variable

  LOGICAL, EXTERNAL :: non_blank
!   Function to detect blank lines

! Pointers to dimensions: used to shorten declarations later
  INTEGER, POINTER :: nd_band
!   Size allocated for spectral bands
  INTEGER, POINTER :: nd_continuum
!   Size allocated for continua
  INTEGER, POINTER :: nd_scale_variable
!   Size allocated for scaling variables

! Alias pointers to dimensions to the actual structure.
  nd_band            => Spectrum%Dim%nd_band
  nd_continuum       => Spectrum%Dim%nd_continuum
  nd_scale_variable  => Spectrum%Dim%nd_scale_variable

! If the block does not exist it is filled with grey null fits.
  IF (.NOT.Spectrum%Basic%l_present(9)) THEN
!   Allocate space for the variables.
    IF (ALLOCATED(Spectrum%Cont%k_cont)) &
        DEALLOCATE(Spectrum%Cont%k_cont)
    ALLOCATE(Spectrum%Cont%k_cont(nd_band, nd_continuum))
    IF (ALLOCATED(Spectrum%Cont%i_scale_fnc_cont)) &
        DEALLOCATE(Spectrum%Cont%i_scale_fnc_cont)
    ALLOCATE(Spectrum%Cont%i_scale_fnc_cont(nd_band, nd_continuum))
    IF (ALLOCATED(Spectrum%Cont%scale_cont)) &
        DEALLOCATE(Spectrum%Cont%scale_cont)
    ALLOCATE(Spectrum%Cont%scale_cont(nd_scale_variable, &
      nd_band, nd_continuum))
    IF (ALLOCATED(Spectrum%Cont%t_ref_cont)) &
        DEALLOCATE(Spectrum%Cont%t_ref_cont)
    ALLOCATE(Spectrum%Cont%t_ref_cont(nd_continuum, nd_band))
    IF (ALLOCATED(Spectrum%Cont%p_ref_cont)) &
        DEALLOCATE(Spectrum%Cont%p_ref_cont)
    ALLOCATE(Spectrum%Cont%p_ref_cont(nd_continuum, nd_band))
    DO i=1, Spectrum%Basic%n_band
      DO j=1, Spectrum%Cont%n_band_continuum(i)
        Spectrum%Cont%k_cont(i, j)=0.0_RealK
        Spectrum%Cont%i_scale_fnc_cont(i, j)=IP_scale_fnc_null
        Spectrum%Cont%p_ref_cont(j, i)=1.0_RealK
        Spectrum%Cont%t_ref_cont(j, i)=200.0_RealK
      END DO
    END DO
  END IF

! Obtain the band data from the prepared file of ESFT terms.
  CALL get_free_unit(ierr, iu_esft)
  CALL open_file_in(ierr, iu_esft, &
    'enter the name of the file containing the continuum fit.')
  DO
    READ(iu_esft, '(A)', IOSTAT=ios) line
    IF (ios /= 0) THEN
      WRITE(*, '(/a)') '***error: file type not found.'
      ierr = i_err_fatal
      RETURN
    END IF
    IF (line(1:10) == '*FILE TYPE') THEN
      BACKSPACE iu_esft
      EXIT
    END IF
  END DO
  IF (ierr /= i_normal) RETURN

! Assemble the list of indexing numbers.
  DO i=1, Spectrum%Basic%n_band
    l_index_band(i, :)=.FALSE.
    DO j=1, Spectrum%Cont%n_band_continuum(i)
      l_index_band(i, Spectrum%Cont%index_continuum(i, j))=.TRUE.
    END DO
  END DO

  outer: DO
    inner: DO
      READ(iu_esft, '(A)', IOSTAT=ios) line
      IF (ios < 0) EXIT outer
      IF (line(1:10) == '*FILE TYPE') THEN
        BACKSPACE iu_esft
        EXIT inner
      END IF
    END DO inner
    READ(iu_esft, '(15x, i5, //)', IOSTAT=ios) i_input_type
    IF (ios < 0) EXIT
    IF (i_input_type /= it_file_cont_fit) THEN
      WRITE(*, '(/a)') &
        '*** error: the input file is of an invalid type.'
      ierr=i_err_fatal
      RETURN
    END IF
    READ(iu_esft, '(14x, i5, 21x, i5)') i_band, i_index

!   Find the position of this datum in the array of continua.
    IF (.NOT.l_index_band(i_band, i_index)) THEN
      WRITE(*, '(/a, i5)') 'Adding continuum to band',i_band
      l_index_band(i_band, i_index) = .TRUE.
      Spectrum%Cont%n_band_continuum(i_band) =                          &
        Spectrum%Cont%n_band_continuum(i_band)+1
      Spectrum%Cont%index_continuum(                                    &
        Spectrum%Cont%n_band_continuum(i_band), i_band) = i_index
    END IF

    READ(iu_esft, '(18x, 1pe10.3, 21x, 1pe10.3)')                       &
      Spectrum%Cont%p_ref_cont(i_index, i_band),                        &
      Spectrum%Cont%t_ref_cont(i_index, i_band)
    READ(iu_esft, '(//)')
!   Read over the transmission data.
    DO
      READ(iu_esft, '(a)') line
      IF (.NOT.non_blank(line)) EXIT
    END DO
    READ(iu_esft, '(/, 74x, i5, //)')                                   &
      Spectrum%Cont%i_scale_fnc_cont(i_band, i_index)
    READ(iu_esft, '(2(3x, 1pe16.9), (t39, 2(3x, 1pe16.9)))')            &
      Spectrum%Cont%k_cont(i_band, i_index), dummy,                     &
      (Spectrum%Cont%scale_cont(l, i_band, i_index), l=1,               &
      n_scale_variable(Spectrum%Cont%i_scale_fnc_cont(i_band, i_index)))
    READ(iu_esft, '(/)')
  END DO outer

  Spectrum%Basic%l_present(9)=.TRUE.
  CLOSE(iu_esft)

END SUBROUTINE make_block_9

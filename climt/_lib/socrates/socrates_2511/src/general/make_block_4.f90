! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 4.
!
! Description:
!   This routine assigns gaseous absorbers to specific bands.
!
!- ----------------------------------------------------------------------------
SUBROUTINE make_block_4(Spectrum, type_index, l_interactive, ierr)

  USE realtype_rd
  USE def_spectrum
  USE def_std_io_icf
  USE gas_list_pcf, ONLY : npd_gases
  USE error_pcf

  IMPLICIT NONE


! Dummy arguments
  INTEGER, Intent(IN) :: type_index(npd_gases)
!   Indices of each gas type
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  TYPE (StrSpecData), Intent(INOUT), Target :: Spectrum
!   Spectral data
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables
  CHARACTER (LEN=35), Parameter :: prompt = &
    'identifiers of absorbers in band '
!   Prompt for input
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: i_list(Spectrum%Dim%nd_species)
!   Indices for one band
!
  TYPE (StrSpecDim), Pointer :: SpDim
!   Pointer to dimensions within the spectrum
  TYPE (StrSpecBasic), Pointer :: SpBasic
!   Pointer to basic components of the spectrum
  TYPE (StrSpecGas), Pointer :: SpGas
!   Pointer to gaseous within the spectrum


  SpDim   => Spectrum%Dim
  ALLOCATE(Spectrum%Gas%n_band_absorb(SpDim%nd_band))
  ALLOCATE(Spectrum%Gas%i_overlap(SpDim%nd_band))
  Spectrum%Gas%i_overlap = 0
  ALLOCATE(Spectrum%Gas%index_absorb(SpDim%nd_species, SpDim%nd_band))
  SpBasic => Spectrum%Basic
  SpGas   => Spectrum%Gas

! Obtain the input data.
  WRITE(iu_stdout, '(//A, /A)') &
    'For each band specify the type numbers of the absorbers', &
    'active in that band.'
  WRITE(iu_stdout, '(A, /A/)') &
    'To continue input on the next line terminate the line with an &.', &
    'enter ''0'' if there is no gaseous absorption in the band.'

  i = 0
  DO
    i = i + 1
    DO
      CALL read_line(ierr, i, prompt, SpDim%nd_species, &
        SpGas%n_band_absorb(i), i_list)
      IF (ierr /= i_normal) RETURN
      IF ((SpGas%n_band_absorb(i) == 0).AND.(i_list(1) > 0)) THEN
        SpGas%n_band_absorb(i+1:i+i_list(1)-1) = 0
        i = i + i_list(1) - 1
        EXIT
      END IF

      j=1
      DO
        IF (j > SpGas%n_band_absorb(i)) EXIT
        SpGas%index_absorb(j, i)=type_index(i_list(j))
        IF (SpGas%index_absorb(j, i) == -1) then
          WRITE(iu_err, '(/A, /A)') &
            '+++ Erroneous input:  An absorber has been specified', &
            'which has not prevoiusly been declared.'
          IF (l_interactive) THEN
            WRITE(iu_err, '(A/)') 'Please re-specify.'
            EXIT
          ELSE
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
        j = j + 1
      ENDDO

      IF (j > SpGas%n_band_absorb(i)) EXIT
    ENDDO
    IF (i >= SpBasic%n_band) EXIT
  ENDDO

END SUBROUTINE make_block_4

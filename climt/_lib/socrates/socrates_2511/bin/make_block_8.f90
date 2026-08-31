! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 8.
!
! Description:
!   This routine defines continua present in the spectral file.
!
!- ----------------------------------------------------------------------------
SUBROUTINE make_block_8(Spectrum, l_interactive, ierr)

  USE realtype_rd
  USE def_spectrum
  USE dimensions_spec_ucf
  USE def_std_io_icf
  USE rad_pcf

  IMPLICIT NONE


! Dummy arguments
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive use
  TYPE (StrSpecData), Intent(INOUT), Target :: Spectrum
!   Spectral data
  integer, Intent(INOUT) :: ierr
!   Error flag

! Local variables
  CHARACTER (LEN=33) :: prompt = 'type numbers of continua in band'
!   Prompt for input
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: i_list(npd_continuum)
!   Indices for one band
!
  TYPE (StrSpecDim), Pointer :: SpDim
!   Pointer to dimensions within the spectrum
  TYPE (StrSpecBasic), Pointer :: SpBasic
!   Pointer to basic components of the spectrum
  TYPE (StrSpecCont), Pointer :: SpCont
!   Pointer to continuum data within the spectrum


  SpDim   => Spectrum%Dim
  ALLOCATE(Spectrum%Cont%n_band_continuum(SpDim%nd_band))
  ALLOCATE(Spectrum%Cont%index_continuum(SpDim%nd_band, npd_continuum))
  SpBasic => Spectrum%Basic
  SpCont  => Spectrum%Cont

! Obtain the input data.
  WRITE(iu_stdout, '(//a, /A)') &
    'For each band, specify the types of the continuum absorbers ', &
    'active in that band.'
  WRITE(iu_stdout, '(a/, a/)') &
   'To continue input on the next line terminate the line with an &.'
  WRITE(iu_stdout, '(/a/)') &
    'Enter "0" if there is no continuum absorption in the band.'

  i = 0
  DO
    i = i + 1
    DO
      CALL read_line(ierr, i, prompt, npd_continuum, &
        SpCont%n_band_continuum(i), i_list)
      IF (ierr /= i_normal) RETURN
      IF ((SpCont%n_band_continuum(i) == 0).AND.(i_list(1) > 0)) THEN
        SpCont%n_band_continuum(i+1:i+i_list(1)-1) = 0
        i = i + i_list(1) - 1
        EXIT
      END IF

      j=1
      DO
        IF (j > SpCont%n_band_continuum(i)) EXIT
        SpCont%index_continuum(i, j)=i_list(j)
        IF (SpCont%index_continuum(i, j) == -1) then
          WRITE(iu_err, '(/a, /a)') &
            '+++ Erroneous input:  a continuum absorber has been ', &
            'specified which is not known.'
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

      IF (j > SpCont%n_band_continuum(i)) EXIT
    ENDDO
    IF (i >= SpBasic%n_band) EXIT
  ENDDO

END SUBROUTINE make_block_8

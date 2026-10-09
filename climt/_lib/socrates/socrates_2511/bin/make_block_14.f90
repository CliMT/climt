! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 14.
!
SUBROUTINE make_block_14 &
!
(Spectrum, l_exclude, l_interactive, ierr)
!
! Description:
!   This determines exclusions to deal with non-contiguous bands.
!
! Method:
!   Straightforward.
!
!- End of header
!
!
! Modules to set types of variables:
  use realtype_rd
  use def_spectrum
  USE def_std_io_icf
  USE error_pcf
!
!
  implicit none
!
!
!
! Dummy arguments
  INTEGER, Intent(INOUT) :: ierr
!       Error flag
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  LOGICAL, Intent(OUT) :: l_exclude
!       Flag for exclusion of bands
  TYPE (StrSpecData), Intent(INOUT), Target :: Spectrum
!   Spectral data
!
!
! Local variables
  character (LEN=25) :: prompt='bands excluded from band '
!       Prompt for input
  CHARACTER (LEN=1) :: char_yn
!       Character response variable
  INTEGER :: i
!       Loop variable
  INTEGER :: j
!       Loop variable
  INTEGER :: index_exclude(Spectrum%Dim%nd_band, Spectrum%Dim%nd_band)
!       Array of indices for absorbers
  TYPE (StrSpecDim), Pointer :: SpDim
!   Pointer to dimensions within the spectrum
  TYPE (StrSpecBasic), Pointer :: SpBasic
!   Pointer to basic components of the spectrum
!
!
  SpDim   => Spectrum%Dim
  ALLOCATE(Spectrum%Basic%n_band_exclude(SpDim%nd_band))
  SpBasic => Spectrum%Basic
!
! Determine whether exclusion will occur.
  l_exclude=.FALSE.
  write(iu_stdout, '(//a)') &
    'Do you wish to exclude regions from particular bands? (y/n)'
  read(iu_stdin, '(a)') char_yn
  if ( (char_yn.ne.'y').and.(char_yn.ne.'Y') ) then
    SpDim%nd_exclude = 1
    ALLOCATE(SpBasic%index_exclude(SpDim%nd_exclude, SpDim%nd_band))
    return
  end if
!
  l_exclude=.TRUE.
  SpDim%nd_exclude=0
  write(iu_stdout, '(a)') &
    'For each band enter the list of bands to be excluded therefrom.'
  write(iu_stdout, '(a)') &
    'To continue input on the next line terminate the line with an &.'
  write(iu_stdout, '(/a/)') &
    'Enter ''0'' if no bands are to be excluded.'
!
  DO i=1, SpBasic%n_band
    DO
      call read_line(ierr, i, prompt, SpDim%nd_band, &
        SpBasic%n_band_exclude(i), index_exclude(1, i))
      if (ierr /= i_normal) return
      if ( (SpBasic%n_band_exclude(i) > SpBasic%n_band) ) then
        write(iu_err, '(/a)') &
          '*** error:  the list of excluded bands is too large.'
        IF (l_interactive) then
          WRITE(iu_err, '(A)') 'Please reenter.'
        ELSE
          ierr=i_err_fatal
          return
        ENDIF
      ELSE
        SpDim%nd_exclude=MAX(SpDim%nd_exclude, SpBasic%n_band_exclude(i))
        EXIT
      ENDIF
    ENDDO
!
    j=1
    DO
      IF (j > SpBasic%n_band_exclude(i)) EXIT
      if ( (index_exclude(j, i) < 1) .OR. &
           (index_exclude(j, i) > SpBasic%n_band) ) then
        write(iu_err, '(/a)') &
          '+++ Erroneous input:  an excluded band is out of range.'
        if (l_interactive) then
          write(iu_err, '(a/)') 'Please re-specify.'
        else
          ierr=i_err_fatal
          return
        endif
      endif
      j=j+1
    ENDDO
  enddo
!
! Transfer to the proper copy.
  ALLOCATE(SpBasic%index_exclude(SpDim%nd_exclude, SpDim%nd_band))
  DO i=1, SpBasic%n_band
    DO j=1, SpBasic%n_band_exclude(i)
      SpBasic%index_exclude(j, i)=index_exclude(j, i)
    ENDDO
  ENDDO
!
!
!
  return
END SUBROUTINE make_block_14

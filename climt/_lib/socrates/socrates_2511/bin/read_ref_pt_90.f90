! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read gaseous reference p and T from a file.
!
SUBROUTINE read_ref_pt_90 &
!
(i_gas, &
 n_selected_band, list_band, &
 p_ref, t_ref, &
 ierr &
)
!
! Method:
!   The name of the file is requested. The file is opened
!   and read. Directives encoutered are processed.
!
!
! Modules used:
  USE realtype_rd
  USE dimensions_spec_ucf
  USE dimensions_pp_ucf
  USE def_std_io_icf
  USE gas_list_pcf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(IN) :: i_gas
!   Types of absorber
  INTEGER, Intent(IN) :: n_selected_band
!   Number of Spectral bands selected
  INTEGER, Intent(IN) :: list_band(:)
!   List of bands in which index falls
  REAL  (RealK), Intent(OUT) :: p_ref(:)
!   Reference pressures for scaling
  REAL  (RealK), Intent(OUT) :: t_ref(:)
!   Reference pressures for scaling
!
!
! Local variables.
  CHARACTER (LEN=80) :: line
!   Line read from input file
  CHARACTER (LEN=16) :: word
!   Substring of input line
  INTEGER :: iu_file_in
!   Unit number for input from file
  INTEGER :: i_begin
!   Beginning of substring
  INTEGER :: i_end
!   End of substring
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i_first_band
!   Fisrt band in list
  INTEGER :: i_last_band
!   Fisrt band in list
  INTEGER :: i_gas_in
!   Type of gas read from input file
  INTEGER :: k
!   Loop variable
  REAL  (RealK) :: p0
!   Pressure read from file
  REAL  (RealK) :: t0
!   Temperature read from file
!
! Subroutines called:
  EXTERNAL &
      open_file_in, read_word
!
!
!
!
  CALL get_free_unit(ierr, iu_file_in)
  CALL open_file_in(ierr, iu_file_in, &
    'Enter the name of the file setting the reference conditions.')
!
! Set the reference conditions in the required bands to negative
! values to flag the absence of data.
  DO k=1, n_selected_band
    p_ref(list_band(k))=-1.0_RealK
    t_ref(list_band(k))=-1.0_RealK
  ENDDO
! Read the file looking for the directive "*REF".
  Input: DO
    READ(iu_file_in, '(a)', IOSTAT=ios) line
!
    IF (ios /= 0) EXIT
!
!   Lines of valid data begin with the directive *REF.
    IF (line(1:4) == '*REF') THEN
!     Process the remainder of the line to extract the data.
!
!     Extract the first band of the range.
      i_begin=6
      CALL read_word(ierr, line, 80, i_begin, i_end, word)
      IF (ierr /= i_normal) RETURN
        READ(word, "(i16)", IOSTAT=ios) i_first_band
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') &
          '*** Error: first band is incorrectly specified in line:'
        WRITE(iu_err, "(a)") line
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Extract the last band in the range.
      i_begin=i_end+1
      CALL read_word(ierr, line, 80, i_begin, i_end, word)
      IF (ierr /= i_normal) RETURN
        READ(word, '(i16)', IOSTAT=ios) i_last_band
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') &
          '*** Error: last band is incorrectly specified in line:'
        WRITE(iu_err, '(a)') line
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Extract the type of the gas. This is in part an inheritance from
!     the earlier version of this routine, where several gases were
!     treated at once, but it also allows all reference conditions to
!     be placed in one file.
      i_begin=i_end+1
      CALL read_word(ierr, line, 80, i_begin, i_end, word)
      IF (ierr /= i_normal) RETURN
      READ(word, '(i16)', IOSTAT=ios) i_gas_in
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') &
          '*** Error: gas is incorrectly specified in line:'
        WRITE(iu_err, '(a)') line
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Extract the pressure.
      i_begin=i_end+1
      CALL read_word(ierr, line, 80, i_begin, i_end, word)
      IF (ierr /= i_normal) RETURN
!     Real values are read with ".0" in the format. A decimal
!     point, if present, will overide the "0"; but if no point
!     is present, no characters are taken as lying to the right
!     thereof.
      READ(word, '(e16.0)', IOSTAT=ios) p0
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') &
          '*** Error: Pressure is incorrectly specified in line:'
        WRITE(iu_err, '(a)') line
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Extract the temperature.
      i_begin=i_end+1
      CALL read_word(ierr, line, 80, i_begin, i_end, word)
      IF (ierr /= i_normal) return
        READ(word, '(e16.0)', IOSTAT=ios) t0
      IF (ios /= 0) THEN
        WRITE(iu_err, '(/a)') &
          '*** Error: Temperature is incorrectly specified in line:'
        WRITE(iu_err, '(a)') line
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      IF (i_gas_in == i_gas) THEN
!       The gas is valid.
        IF (i_last_band > SIZE(p_ref)) THEN
          WRITE(*, '(/a)') &
            'Warning: too many bands in reference P/T file.'
          i_last_band=SIZE(p_ref)
        END IF
        DO k=i_first_band, i_last_band
           p_ref(k)=p0
           t_ref(k)=t0
         ENDDO
      ENDIF
!
    ENDIF
  ENDDO Input
!
  CLOSE(iu_file_in)
!
! Check that reference conditions are set in all the required bands.
  DO k=1, n_selected_band
    IF ( (p_ref(list_band(k)) < 0.0_RealK) .OR.  &
         (t_ref(list_band(k)) < 0.0_RealK) ) THEN
      WRITE(iu_err, '(/a, i5)') &
        '*** Error: reference conditions are missing for  band ', &
        list_band(k)
      ierr=i_err_fatal
      RETURN
    ENDIF
  ENDDO
!
!
!
  RETURN
END

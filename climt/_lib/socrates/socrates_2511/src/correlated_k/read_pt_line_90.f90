! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read and process a line giving p and T.
!
SUBROUTINE read_pt_line_90 &
!
(ierr, l_interactive, l_file, iu_file_in, &
 nd_pt, n_pt_pair, p, t, l_finish)
!
! Method:
!  Straightforward.
!
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  LOGICAL, Intent(IN) :: l_interactive
!   Flag for interactive operation
  LOGICAL, Intent(IN) :: l_file
!   Flag set as true when reading from a file
  INTEGER, Intent(IN) :: iu_file_in
!   Unit assigned for input from a file
  INTEGER, Intent(IN) :: nd_pt
!   Size allocated for number of pairs of p and T
  INTEGER, Intent(INOUT) :: n_pt_pair
!   Number of pairs
  LOGICAL, Intent(INOUT) :: l_finish
!   Termination flag
  REAL  (RealK), Intent(OUT) :: p(nd_pt)
!   Values of pressures.
  REAL  (RealK), Intent(OUT) :: t(nd_pt)
!   Values of temperatures.
!
!
! Local variables
  INTEGER :: length
!   Length of input line
  INTEGER :: n_word
!   Number of words in line
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: begin
!   Beginning of string
  INTEGER :: end
!   End of string
  INTEGER :: iu_pt
!   Unit number assigned for input of p and T
  INTEGER :: ios
!   I/O error flag
  INTEGER,PARAMETER :: n_char_line=1024
!   Maximum number of characters in a line
  INTEGER,PARAMETER :: n_char_list=2048
!   Maximum number of characters in a list
  INTEGER,PARAMETER :: n_char_word=20
!   Maximum number of characters in a word
  CHARACTER (LEN=n_char_line) :: line
!   Line of input data
  CHARACTER (LEN=n_char_list) :: list
!   Assembled input list
  CHARACTER (LEN=n_char_word) :: word
!   Input number
  LOGICAL :: l_next
!   Flag to read next line
  LOGICAL :: l_reread
!   Flag to reread the line
  REAL  (RealK), Dimension(nd_pt) :: x
!   Words read in
!
!
!
! Set unit numbers for I/O.
  IF (l_file) THEN
    iu_pt = iu_file_in
  ELSE
    iu_pt = iu_stdin
  ENDIF
!
  DO
!
!   Initialize input variables.
    list(1:1)=' '
    begin=2
!
    IF (.NOT.l_file) THEN
      WRITE(iu_stdout, '(a)') &
        'Specify pressure and corresponding temperatures (*END to finish)'
    ENDIF
!
    l_next= .TRUE. 
!   The next line will be read until l_next is false.
    Input: DO
      IF (.NOT.l_next) EXIT
      READ(iu_pt, '(a)', IOSTAT=ios) line
      IF (ios /= 0) THEN
        WRITE(iu_err, "(a)") "Erroneous input"
        IF ( l_interactive .AND. (.NOT.l_file) ) THEN
          WRITE(iu_stdout, '(a)') "Please re-enter this line."
          CYCLE
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
!
!     Check whether concatenation is required.
      end=n_char_line
      DO
        IF (end == 1) THEN
          EXIT
        ELSE IF (line(end:end) == ' ') THEN
          end=end-1
        ELSE
          EXIT
        ENDIF
      ENDDO
      IF (line(end:end) == '&') THEN
!       Set the return flag.
        l_next= .TRUE.
        end=end-1
      ELSE
        l_next= .FALSE.
      ENDIF
!
!     Add this line on to the long list.
      IF (end >= 1) THEN
        list(begin: begin+end)=' '//line(1:end)
!       Advance the free counter.
        begin=begin+end+1
      ENDIF
!
    ENDDO Input
!
!
!   Determine the length of the input line and process to determine
!   the numbers contained therein.
    list(begin:begin)=' '
    l_reread=.FALSE.
    length=begin
    n_word=0
    j=2
    Process: DO
!
      IF (j >= length) EXIT
!
      IF ( (list(j-1:j-1) == ' ') .AND. (list(j:j) /= ' ') ) THEN
!       Beginning of word found.
        begin=j
        n_word=n_word+1
!       Check for termination of input: a line begins with 'F' 
!       or 'f' or the directive "*END"
        IF (n_word == 1) THEN
          IF ( (list(j:j) == 'f') .OR. &
               (list(j:j) == 'f') .OR.  &
               (list(j:j+3) == '*END') ) THEN  
            l_finish= .TRUE. 
            RETURN
          ENDIF
        ENDIF
      ENDIF
!
      IF ( (list(j+1:j+1) == ' ') .AND. (list(j:j) /= ' ') ) THEN
!       End of word found.
        end=j
!
!       Check and perform read.
        IF ( (end-begin) > n_char_word-1 ) THEN
          WRITE(iu_err, '(/a)') &
            '+++ Erroneous input: a number is too large: '
          IF (l_file .OR. (.NOT.l_interactive) ) THEN
            ierr=i_err_fatal
            RETURN
          ELSE
            WRITE(iu_stdout, '(a)') 'Please re-enter the list.'
            l_reread=.TRUE.
            EXIT
          ENDIF
        ENDIF
        word(1:n_char_word) = '            '
        word(n_char_word+begin-end:n_char_word)=list(begin:end)
!
!       Read in the word.
        READ(word, FMT=*, IOSTAT=ios) x(n_word)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(/a/)') '+++ erroneous input:'
          IF (l_file .OR. (.NOT.l_interactive) ) THEN
            ierr=i_err_fatal
            RETURN
          ELSE
            WRITE(iu_stdout, '(/a/)') 'Please re-enter'
            l_reread=.TRUE.
            EXIT
          ENDIF
        ENDIF
      ENDIF
      j=j+1
    ENDDO Process
!
    IF (l_reread) CYCLE
!
!   Transfer the numbers read in to the arrays of P and T.
!   check array sizes.
    IF (n_word-1+n_pt_pair > nd_pt) THEN
      WRITE(iu_err, '(/a22, a)') &
        '+++ Too many pairs of pressure and temperature. '
      WRITE(iu_err, '(a, i5)') &
        'Maximum number allowed = ', nd_pt
      IF (l_file .OR. (.NOT.l_interactive) ) THEN
        ierr=i_err_fatal
        RETURN
      ELSE
        WRITE(iu_stdout, '(a)') 'Please re-specify this last line.'
      ENDIF
    ELSE IF (n_word == 1) THEN
      WRITE(iu_err, '(/a)') '+++ No temperatures given:'
      IF (l_file .OR. (.NOT.l_interactive) ) THEN
        ierr=i_err_fatal
        RETURN
      ELSE
        WRITE(iu_stdout, '(/a)') 'Re-enter the list.'
      ENDIF
    ELSE
      DO i=2, n_word
        n_pt_pair=n_pt_pair+1
        p(n_pt_pair)=x(1)
        t(n_pt_pair)=x(i)
      ENDDO
      EXIT
    ENDIF
!
  ENDDO
!
!
!
  RETURN
END

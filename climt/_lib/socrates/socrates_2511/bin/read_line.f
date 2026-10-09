! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a line giving the gases in a band.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_line(ierr, i, prompt, max_size
     &  , n_abs_band, i_abs_band
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
      CHARACTER !, Intent(IN)
     &    prompt*(*)
!           Prompt to user
      INTEGER, Intent(IN) ::
     &    i
!           Number of band
     &  , max_size
!           Maximum array size
!
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , n_abs_band
!           Number of absorbers in band
     &  , i_abs_band(max_size)
!           Indexing numbers
!
!     Local varaibles
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      INTEGER
     &    length
!           Length of input line
     &  , j
!           Loop variable
     &  , k
!           Loop variable
     &  , begin
!           Beginning of string
     &  , END
!           End of string
     &  , ios
!           Io error flag
      CHARACTER
     &    line*80
!           Line of input data
     &  , list*1024
!           Assembled input list
     &  , word*5
!           Input number
      LOGICAL
     &    l_next
!           Flag to read next line
!
!
!
!     Initialize input variables.
1     list(1:1)=' '
      begin=2
!
      WRITE(*, '(a, 1x, i5)')
     &  prompt, i
2     read(iu_stdin, '(a)', iostat=ios) line
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ erroneous input '
        IF (lock_code(.true.)) THEN
          ierr=i_err_fatal
          RETURN
        ELSE
          WRITE(*, '(a)') ' please re-enter this line.'
          goto 2
        ENDIF
      ENDIF

      IF (line(1:2).eq.'0*') THEN
        n_abs_band=0
        READ(line(3:), fmt='(i5)', iostat=ios) i_abs_band(1)
        RETURN
      END IF

!     Check whether concatenation is required.
      END=80
3     if (line(end:end).eq.' ') then
        IF (end > 1) THEN
          END=end-1
          goto 3
        ELSE
          WRITE(iu_err, '(a)') '+++ erroneous input '
          IF (lock_code(.true.)) THEN
          ierr=i_err_fatal
            RETURN
          ELSE
            WRITE(*, '(a)') ' please re-enter this line.'
            goto 2
          ENDIF
        ENDIF
      ENDIF
      IF (line(end:end) == '&') THEN
!       Set the return flag.
        l_next=.true.
        END=end-1
      ELSE
        l_next=.false.
      ENDIF
!
!     Add this line on to the long list.
      list(begin: begin+end)=' '//line(1:end)
!     Advance the free counter.
      begin=begin+end+1
!     Read the next line if necessary.
      IF (l_next) goto 2
!
!     Determine the length of the input line and process to determine
!     the numbers contained therein.
      list(begin:begin)=' '
      length=begin
      n_abs_band=0
      j=2
4     continue
        IF ( (list(j-1:j-1) == ' ').AND.(list(j:j) /= ' ') ) THEN
!         Beginning of word found.
          begin=j
          n_abs_band=n_abs_band+1
!
5         continue
            IF ( (list(j+1:j+1) == ' ').AND.
     &           (list(j:j) /= ' ') ) THEN
!             End of word found.
              END=j
!
!             Check and perform read.
              IF ( (end-begin) > 4 ) THEN
                WRITE(iu_err, '(/a)')
     &            '+++ erroneous input: a number is too large: '
                IF (lock_code(.true.)) THEN
                  ierr=i_err_fatal
                  RETURN
                ELSE
                  WRITE(iu_err, '(a/)')
     &              'please re-enter the data for this band.'
                  goto 1
                ENDIF
              ENDIF
              DO k=1, 5
                word(k:k)=' '
              ENDDO
              DO k=5+begin-end, 5
                word(k:k)=list(end-5+k:end-5+k)
              ENDDO
!
!             Read in the word.
              READ(word, fmt='(i5)', iostat=ios)
     &          i_abs_band(n_abs_band)
              IF (ios /= 0) THEN
                WRITE(iu_err, '(/a)')
     &            '+++ erroneous input:'
                IF (lock_code(.true.)) THEN
                  ierr=i_err_fatal
                  RETURN
                ELSE
                  WRITE(iu_err, '(a/)')
     &              'please re-enter the data for this band.'
                  goto 1
                ENDIF
              ENDIF
              j=j+1
              goto 4
!
          ELSE
            j=j+1
            goto 5
          ENDIF
        ELSE IF (j < length) THEN
          j=j+1
          goto 4
        ELSE
          continue
        ENDIF
!
!     Reset N_ABS_BAND if there are actually no absorbers active in
!     the band.
      IF (i_abs_band(1) == 0) THEN
        n_abs_band=0
      ENDIF
!
!
!
      RETURN
      END

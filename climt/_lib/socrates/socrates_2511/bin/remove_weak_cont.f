! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to remove weakly absorbing comntinua from bands.
!
! Method:
!   A tolerance for negelecting continuum absorption is obtained.
!   A typical amount of each continuum in a column is given: for
!   the self-broadened continuum this must be multiplied by a
!   broadening density. For each band a test is made to see 
!   whether this amount of continuum absorber gives a transmission 
!   close to 1 as defined by the tolerance. If so, the continuum is
!   removed from this spectral band.
!
!- ---------------------------------------------------------------------
      SUBROUTINE REMOVE_WEAK_CONT(IERR, N_BAND
     &  , N_BAND_CONTINUUM, INDEX_CONTINUUM, K_CONTINUUM
     &  )
!
!
!
!     Modules to set types of variables:
      USE REALTYPE_RD
      USE rad_pcf
      USE dimensions_spec_ucf
      USE def_std_io_icf


      IMPLICIT NONE


!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::
     &    N_BAND
!           Number of bands
      INTEGER, INTENT(INOUT) ::
     &    N_BAND_CONTINUUM(NPD_BAND)
!           Number of continua in bands
     &  , INDEX_CONTINUUM(NPD_BAND, NPD_SPECIES)
!           List of active continua
      REAL  (RealK), INTENT(IN) ::
     &    K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)
!           Continuum extinction coefficients
!
!     Local arguments.
!
      LOGICAL
     &    LOCK_CODE
!           Logical to forbid interactive looping
      EXTERNAL
     &    LOCK_CODE
!
      INTEGER
     &    I_CONTINUUM
!           Index number of continuum
     &  , I
!           Loop variable
     &  , J
!           Loop variable
     &  , K
!           Loop variable
     &  , N_BAND_CONTINUUM_TEMP
!           Temporary number of continua
     &  , IOS
!           I/O error flag
      LOGICAL
     &    L_SET_CONTINUUM
!           Flag for setting amounts
      REAL  (RealK) ::
     &    COLUMN(NPD_CONTINUUM)
!           Column amounts for continua
     &  , TRANS_NEGLECT
!           Transmission for neglecting
     &  , TRANS_COLUMN
!           Transmission of column
!
!
!     Obtain column amounts of continua.
!     Self-broadened continuum of water vapour.
      L_SET_CONTINUUM=.FALSE.
      DO J=1, N_BAND
        DO K=1, N_BAND_CONTINUUM(J)
          IF ( (INDEX_CONTINUUM(J, K).EQ.IP_SELF_CONTINUUM).AND.
     &         (.NOT.L_SET_CONTINUUM) ) THEN
            WRITE(IU_USER, '(/A, /A, A)') 'ENTER AMOUNT OF WATER '
     &        //'VAPOUR TIMES MOLAR DENSITY OF WATER VAPOUR'
     &        , 'TO TEST FOR NEGLECTING THE SELF-BROADENED '
     &        , 'CONTINUUM IN A BAND.'
1           READ(IU_STDIN, *, IOSTAT=IOS) COLUMN(IP_SELF_CONTINUUM)
            IF (IOS.NE.0) THEN
              WRITE(IU_ERR, '(A)') '+++ ERRONEOUS RESPONSE:'
              IF (LOCK_CODE(.TRUE.)) THEN
                INCLUDE 'batch_error_main.finc'
              ELSE
                WRITE(IU_USER, '(A)') 'PLEASE RE-TYPE.'
                GOTO 1
              ENDIF
            ENDIF
            L_SET_CONTINUUM=.TRUE.
          ENDIF
        ENDDO
      ENDDO
!
!     Foreign broadened continuum of water vapour.
      L_SET_CONTINUUM=.FALSE.
      DO J=1, N_BAND
        DO K=1, N_BAND_CONTINUUM(J)
          IF ( (INDEX_CONTINUUM(J, K).EQ.IP_FRN_CONTINUUM).AND.
     &         (.NOT.L_SET_CONTINUUM) ) THEN
            WRITE(IU_USER, '(/A, A, /A)') 'ENTER AMOUNT OF WATER '
     &        , 'VAPOUR TO TEST FOR NEGLECTING THE '
     &        , 'FOREIGN-BROADENED CONTINUUM IN A BAND.'
2           READ(IU_STDIN, *, IOSTAT=IOS) COLUMN(IP_FRN_CONTINUUM)
            IF (IOS.NE.0) THEN
              WRITE(IU_ERR, '(A)') '+++ ERRONEOUS RESPONSE:'
              IF (LOCK_CODE(.TRUE.)) THEN
                INCLUDE 'batch_error_main.finc'
              ELSE
                WRITE(IU_USER, '(A)') 'PLEASE RE-TYPE.'
                GOTO 2
              ENDIF
            ENDIF
            L_SET_CONTINUUM=.TRUE.
          ENDIF
        ENDDO
      ENDDO
!
      WRITE(IU_USER, '(/A)')
     &  'ENTER THE TRANSMISSION FOR NEGLECTING CONTINUA.'
4     READ(IU_STDIN, *, IOSTAT=IOS) TRANS_NEGLECT
      IF (IOS.NE.0) THEN
        WRITE(IU_ERR, '(A)') '+++ ERRONEOUS RESPONSE:'
        IF (LOCK_CODE(.TRUE.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(IU_USER, '(A)') 'PLEASE RE-TYPE.'
          GOTO 4
        ENDIF
      ENDIF
!
!     Go through the bands removing continua which are too weak.
      DO I=1, N_BAND
        N_BAND_CONTINUUM_TEMP=N_BAND_CONTINUUM(I)
        N_BAND_CONTINUUM(I)=0
        DO J=1, N_BAND_CONTINUUM_TEMP
          I_CONTINUUM=INDEX_CONTINUUM(I, J)
          TRANS_COLUMN=EXP(-K_CONTINUUM(I, I_CONTINUUM)
     &      *COLUMN(I_CONTINUUM))
          IF (TRANS_COLUMN.LT.TRANS_NEGLECT) THEN
!           The continuum is included.
            N_BAND_CONTINUUM(I)=N_BAND_CONTINUUM(I)+1
            INDEX_CONTINUUM(I, N_BAND_CONTINUUM(I))=I_CONTINUUM
          ENDIF
        ENDDO
      ENDDO
!
!
!
      RETURN
      END

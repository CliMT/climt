! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 2.
!
! Method:
!   Perform integration using Simpsons Rule.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SIMPSONS_RULE(N,DELTA,INPUT_ARRAY,RESULT)


      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     &     N
!             Number of data points

      REAL(KIND=KIND(1.0D0)) , INTENT(IN) ::
     &     DELTA
!             Size of Interval
     &  ,  INPUT_ARRAY(0:N)
!             Array being summed over

      REAL(KIND=KIND(1.0D0)) , INTENT(OUT) ::
     &     Result
!             The final sum of all array entries

!
! Local Variables
!
      INTEGER
     &   I
!            Looping variable

      REAL (KIND=KIND(1.0D0)) ::
     &    SUM
!            The sum of the array entries


      SUM=0.0
      DO I=0,N
         IF ((I.EQ.0).OR.(I.EQ.N)) THEN
            SUM=SUM+INPUT_ARRAY(I)
         ELSE IF (MOD(I,2).EQ.0) THEN
            SUM=SUM+2.0*INPUT_ARRAY(I)
         ELSE IF (MOD(I,2).EQ.1) THEN
            SUM=SUM+4.0*INPUT_ARRAY(I)
         ENDIF
      ENDDO

      RESULT=SUM*DELTA/3.0

      RETURN
      END

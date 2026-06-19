! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to define a map for sorting.
!
FUNCTION map_shell &
!
(a) &
!
RESULT (map)
!
! Description:
!   This routine receives an input array and returns a key sorted in 
!   increasing order of the input array.
!
! Method:
!   The standard shell sorting algorithm is used with 
!   the gap changing by a factor of 3.
!
! Modules used:
  USE realtype_rd
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments
  REAL  (RealK), Intent(IN), Dimension(:) :: a
!   Array to define the sorting.
!
  INTEGER, Pointer, Dimension(:) :: map
!   Mapping defining the sorted array.
!
!
! Local variables
  INTEGER :: gap
!   Gap between elemnts in list to be sorted
  INTEGER :: j
!   Loop variable
  INTEGER :: k
!   Loop variable
  INTEGER :: itemp
!   Temporary storage for swapping
!
!
!
! Initialize the mapping.
  ALLOCATE(map(SIZE(a)))
  DO j=1, SIZE(a)
    map(j)=j
  ENDDO
!
!
  gap=1
  DO 
    gap=3*gap+1
    IF(gap > SIZE(a)) EXIT
  END DO
!
  DO 
    gap=gap/3
    IF (gap < 1) EXIT
    DO j=gap, SIZE(a)-1
      DO k=j-gap+1, 1, -gap
        IF ( a(map(k)) > a(map(k+gap)) ) THEN
          itemp=map(k)
          map(k)=map(k+gap)
          map(k+gap)=itemp
        ENDIF
      ENDDO
    ENDDO
  ENDDO
!
!
!
  RETURN
END FUNCTION map_shell

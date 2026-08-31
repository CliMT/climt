! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to define a map for sorting.
!
SUBROUTINE map_heap_func &
!
(arr, map)
!
! Description:
!   This routine sorts an array arr(1:n) by index into ascending
!   numerical order using the Heapsort algorithm. 
!   The output is an array of indexes giving the order of the sorted
!   array elements.
!
! Method:
!   A standard heapsort sorting algorithm is adapted to
!   be a function that gives an array of indeces showing the
!   order of the sorted data in the input array.
!   To get the the actual array sorted print arr(map)
!   
!
! Modules used:
  USE realtype_rd
!
  IMPLICIT NONE
!
!
! Dummy arguments
  REAL (RealK), Intent(IN), Dimension(:) :: arr
!
  INTEGER, Intent(OUT), Dimension(:) :: map
!   Mapping defining the sorted array.
!
! Local variables
  INTEGER :: i, j
!   Loop indeces
  INTEGER :: n
!   The size of the input array
  REAL (RealK), Allocatable, Dimension(:) :: arr_copy 
!   A copy of the input array to change and do testing on.
!
! Initialize the size variable.
  n=size(arr)
!
! Initialize the data array.
  ALLOCATE(arr_copy(n))
  arr_copy = arr
!   
! Initialize the mapping by indexing the array data.
  DO j=1, n
    map(j) = j
  ENDDO


  DO i =n/2, 1, -1
    ! The index i, which here determines the "left" range of the sift-down,
    ! i.e., the element to be sifted down, is decremented from n/2 down to 1 
    ! during the "hiring" (heap creation) phase.

    CALL sift_down(i,n)
  END DO

  DO i=n,1,-1
    ! Here the "right" range of the sift down is decremented from n-1 down
    ! to 1 during the "retirement-and-promotion" (heap selection phase).
    
    CALL swap(arr_copy(1),arr_copy(i))      ! Clear a space at the end of the array, and
    CALL sift_down(1,i-1)                   ! retire the top of the heap into it.
  END DO

  DEALLOCATE(arr_copy)


  CONTAINS


  SUBROUTINE sift_down(l,r)
    INTEGER, INTENT(IN) :: l,r
      !  Carry out the sift down on element arr_copy(l) to maintain the heap structure.
    Integer :: j, jold
    REAL (RealK) :: a
    INTEGER      :: saved_map_index
  
    a               = arr_copy(l)
    saved_map_index = map(l)
    jold = l
    j    = l+l

    DO
      IF (j > r) EXIT
      IF (j < r) THEN
        IF (arr_copy(j) < arr_copy(j+1)) THEN
          j = j+1                       ! Compare to better underling.
        ENDIF
      ENDIF

      IF(a >= arr_copy(j)) EXIT         ! Found a's level. Terminate the sift-down.

      arr_copy(jold) = arr_copy(j)
      map(jold)      = map(j)           ! Otherwise, demote a and continue.

      jold      = j
      j         = j+j

    ENDDO

    arr_copy(jold) = a                  ! Put the data into its slot.
    map(jold)      = saved_map_index    ! Make the corresponding change to the index.

  END SUBROUTINE sift_down
  
  SUBROUTINE swap(a,b)
    ! Swap the contents of a and b.
    REAL (RealK) ,INTENT(INOUT) :: a,b
    REAL :: dum
    INTEGER :: map_dum

    dum=a
    a=b
    b=dum

    map_dum = map(1) 
    map(1)  = map(i)
    map(i)  = map_dum

  END SUBROUTINE swap

END SUBROUTINE map_heap_func



! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Standard quicksort routine with  
! 
! Method:
!   Uses the quicksort algorithm to sort an array a_array. The
!   corresponding rearrangement is also performed on the input array
!   b_array. Insertion sort is used for small subarrays to increase
!   efficiency.
!
!- ---------------------------------------------------------------------
MODULE quicksort_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE quicksort(n_elem, a_array, b_array)

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_normal, i_err_fatal
  USE errormessagelength_mod, ONLY:errormessagelength
  USE ereport_mod, ONLY:ereport

  IMPLICIT NONE

  INTEGER,INTENT(IN) ::                                                        &
      n_elem
!     Number of elements in arrays
  REAL(RealK),INTENT(INOUT) ::                                                 &
      a_array(n_elem)                                                          &
!     Array to be sorted in increasing order
    , b_array(n_elem)
!     The rearrangement of b_array is also performed on this array

! Local variables
  INTEGER,PARAMETER ::                                                         &
      n_elem_isort = 7                                                         &
!     Maximum length of subarray for which to use insertion sort
    , n_elem_stack = 50
!     Stack size

! Temporary variables
  INTEGER ::                                                                   &
      i, ir, j, k, l                                                           &
!     Array and loop indices
    , j_stack, i_stack(n_elem_stack)
!     Stack index and array
  REAL(RealK) ::                                                               &
      a_elem, b_elem, temp_elem
!     Array elements

  INTEGER                            :: ierr = i_normal
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'QUICKSORT'

! Initialise indices
  j_stack = 0
  l = 1
  ir = n_elem

  Outer_loop: DO

!   Use insertion sort for small subarrays
    IF (ir-l < n_elem_isort) THEN
      DO j=l+1,ir
        a_elem = a_array(j)
        b_elem = b_array(j)

        i = j - 1
        DO WHILE (i >= l .AND. a_array(i) > a_elem)
          a_array(i+1) = a_array(i)
          b_array(i+1) = b_array(i)
          i = i - 1
        END DO
        a_array(i+1) = a_elem
        b_array(i+1) = b_elem
      ENDDO

!     Check if array is sorted
      IF (j_stack == 0) EXIT Outer_loop

      ir = i_stack(j_stack)
      l = i_stack(j_stack-1)
      j_stack = j_stack-2

    ELSE

!     Choose partitioning element a_elem and perform sort so that
!     a_array(l) <= a_array(l+1) <= a_array(ir).
      k = (l+ir)/2
      temp_elem = a_array(k)
      a_array(k) = a_array(l+1)
      a_array(l+1) = temp_elem
      temp_elem = b_array(k)
      b_array(k) = b_array(l+1)
      b_array(l+1) = temp_elem
      IF (a_array(l) > a_array(ir)) THEN
        temp_elem = a_array(l)
        a_array(l) = a_array(ir)
        a_array(ir) = temp_elem
        temp_elem = b_array(l)
        b_array(l) = b_array(ir)
        b_array(ir) = temp_elem
      ENDIF
      IF (a_array(l+1) > a_array(ir)) THEN
        temp_elem = a_array(l+1)
        a_array(l+1) = a_array(ir)
        a_array(ir) = temp_elem
        temp_elem = b_array(l+1)
        b_array(l+1) = b_array(ir)
        b_array(ir) = temp_elem
      ENDIF
      IF (a_array(l) > a_array(l+1)) THEN
        temp_elem = a_array(l)
        a_array(l) = a_array(l+1)
        a_array(l+1) = temp_elem
        temp_elem = b_array(l)
        b_array(l) = b_array(l+1)
        b_array(l+1) = temp_elem
      ENDIF

!     Set start values for indices before partitioning
      i = l+1
      j = ir

!     Start partitioning
      a_elem = a_array(l+1)
      b_elem = b_array(l+1)

      Inner_loop: DO

!       Search upwards to locate element grater than a_elem
        i = i+1
        DO WHILE (a_array(i) < a_elem)
          i = i+1
        END DO

!       Search downwards to find element smaller than a_elem
        j = j-1
        DO WHILE (a_array(j) > a_elem)
          j = j-1
        END DO

!       Check if partitioning is complete
        IF (j < i) EXIT Inner_loop

!       Swap elements in arrays
        temp_elem = a_array(i)
        a_array(i) = a_array(j)
        a_array(j) = temp_elem
        temp_elem = b_array(i)
        b_array(i) = b_array(j)
        b_array(j) = temp_elem

      END DO Inner_loop

!     Put partitioning elements a_elem and b_elem into arrays
      a_array(l+1) = a_array(j)
      a_array(j) = a_elem
      b_array(l+1) = b_array(j)
      b_array(j) = b_elem

!     Iterate indices for the larger subarray on stack
      j_stack = j_stack+2

!     Check if stack size has been exceeded
      IF (j_stack > n_elem_stack) THEN
        WRITE(cmessage,'(A,I3,A)') 'Stack size is too small.' //  &
          ' Please increase n_elem_stack to at least ', j_stack, '.'
        ierr = i_err_fatal
        EXIT Outer_loop
      END IF

      IF (ir-i+1 >= j-l) THEN
        i_stack(j_stack) = ir
        i_stack(j_stack-1) = i
        ir = j-1
      ELSE
        i_stack(j_stack) = j-1
        i_stack(j_stack-1) = l
        l = i
      ENDIF
    ENDIF

  END DO Outer_loop

  IF (ierr /= i_normal) THEN
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

END SUBROUTINE quicksort
END MODULE quicksort_mod

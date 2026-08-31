! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read an instrumental response from function.
!
SUBROUTINE read_instrument_response_90 &
!
(filter, ierr)
!
! Description:
!
! Reads file containing instrument spectral response data
! 
! Method:
!
! Straightforward
!
! Modules used:    
  USE def_inst_flt
  USE def_std_io_icf
  USE unit_list_pcf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
! Subroutine arguments
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  TYPE  (StrFiltResp), Intent(OUT) :: filter
!   Instrumental response function
!
!
! Local scalars:
  INTEGER :: iu_flt
!   Unit number for input of the filter function
  CHARACTER  (LEN=20) :: ch_unit
!   Character string holding the definition of the unit
  INTEGER :: iunit
!   Index of the unit used for the wavelength or frequency 
  INTEGER :: i
!   Loop variable
!
! Local arrays
  INTEGER, Allocatable, Dimension(:) :: map
!   Map defing the list of points in increasing order or wavenumber
  REAL  (RealK), Allocatable, Dimension(:) :: dum
!   Dummy array for reordering loops
!
! Functions called
!
  INTERFACE
!
   SUBROUTINE map_heap_func(a, map)
!
      USE realtype_rd
!
      REAL  (RealK), Intent(IN), Dimension(:) :: a
!
      INTEGER, Intent(OUT), Dimension(:) :: map
!
    END SUBROUTINE map_heap_func
!
  END INTERFACE
!
!- End of header  
!
!
!
! Open the file containing the instrumental response.
  CALL get_free_unit(ierr, iu_flt)
  IF (ierr /= i_normal) RETURN
  CALL open_file_in(ierr, iu_flt, "Give the instrumental response.")
  IF (ierr /= i_normal) RETURN
!
!
   READ(iu_flt,*)
   READ(iu_flt,"(a)") filter%satellite
   READ(iu_flt,*)
   READ(iu_flt,"(a)") filter%channel
   READ(iu_flt,*)
   READ(iu_flt,"(a)") ch_unit
   IF (ch_unit(1:6) == "*METRE") THEN
     iunit=IP_unit_metre
   ELSE IF (ch_unit(1:7) == "*INV_CM") THEN
     iunit=IP_unit_inverse_cm
   ELSE IF (ch_unit(1:7) == "*MICRON") THEN
     iunit=IP_unit_micron
   ELSE IF (ch_unit(1:10) == "*NANOMETRE") THEN
     iunit=IP_unit_nanometre
   ELSE
     WRITE(iu_err, "(/a)") "*** Error: Invalid unit for frequency."
     ierr=i_err_fatal
     return
   ENDIF
   READ(iu_flt,*)
!
   READ(iu_flt,*) filter%n_pts
   ALLOCATE(filter%wavenumber(filter%n_pts))
   ALLOCATE(filter%response(filter%n_pts))
   ALLOCATE(filter%d2_response(filter%n_pts))
!
   DO i = 1, filter%n_pts
     READ(iu_flt, *) filter%wavenumber(i),  filter%response(i)
!    Convert to m-1 in SI.
     SELECT CASE(iunit)
       CASE(IP_unit_metre)
         filter%wavenumber(i) = 1.0_RealK / filter%wavenumber(i)
       CASE(IP_unit_inverse_cm)
         filter%wavenumber(i) = 1.0E+02_RealK * filter%wavenumber(i)
       CASE(IP_unit_micron)
         filter%wavenumber(i) = 1.0E+06_RealK / filter%wavenumber(i)
       CASE(IP_unit_nanometre)
         filter%wavenumber(i) = 1.0E+09_RealK / filter%wavenumber(i)
     END SELECT
   ENDDO
!
!  Rearrange the filter function in order of increasing wavenumber.
   ALLOCATE(dum(filter%n_pts))
   ALLOCATE(map(filter%n_pts))
   CALL map_heap_func(filter%wavenumber(1:filter%n_pts), map)
   DO i=1, filter%n_pts
     dum(i)=filter%wavenumber(map(i))
   ENDDO
   filter%wavenumber(1:filter%n_pts)=dum(1:filter%n_pts)
   DO i=1, filter%n_pts
     dum(i)=filter%response(map(i))
   ENDDO
   filter%response(1:filter%n_pts)=dum(1:filter%n_pts)
   DEALLOCATE(dum)
   DEALLOCATE(map)
!
!  Zero the second derivative for now: eventually, we may
!  move to spline fitting.
   filter%d2_response(1:filter%n_pts)=0.0_RealK
!
  CLOSE (iu_flt)
!
!
!
  RETURN
END SUBROUTINE read_instrument_response_90

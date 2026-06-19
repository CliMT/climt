! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for MCICA data.
!
! Description:
!   This module contains the declaration of the structure
!   used to store MCICA data in the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_mcica

USE realtype_rd, ONLY: RealK

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DEF_MCICA'

INTEGER, PARAMETER :: ip_mcica_full_sampling = 0
! Each k-term "sees" every sub-column
INTEGER, PARAMETER :: ip_mcica_single_sampling = 1
! Each k-term "sees" a single sub-column
INTEGER, PARAMETER :: ip_mcica_optimal_sampling = 2
! Each k-term "sees" an optimal number of sub-columns

TYPE StrMcica
  INTEGER :: n_subcol_gen = 0
!   Number of sub-columns to generate
  INTEGER :: n_subcol_req_single = 0
!   Number of cloudy sub-columns required for single sampling
  INTEGER :: n_subcol_req_optimal = 0
!   Number of cloudy sub-columns required for optimal sampling
  INTEGER :: ipph = 0
!   Plane-parallel homogeneous flag

  INTEGER, ALLOCATABLE :: lw_subcol_reorder_single(:)
!   LW order of sub-columns (for single sampling)
  INTEGER, ALLOCATABLE :: lw_subcol_reorder_optimal(:)
!   LW order of sub-columns (for optimal sampling)
  INTEGER, ALLOCATABLE :: sw_subcol_k(:, :)
!   Number of subcolumns for each k_term in each band.
  INTEGER, ALLOCATABLE :: lw_subcol_k(:, :)
!   Number of subcolumns for each k_term in each band.

  INTEGER :: n1 = 0
  INTEGER :: n2 = 0
!   Dimensions of xcw array:
!     Cumulative probability (n1)
!     Relative standard deviation (n2)
  REAL (RealK), ALLOCATABLE :: xcw(:, :)
!   Distribution of normalised condensate amount as a function of
!   cumulative probability and relative standard deviation.
END TYPE StrMcica


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE read_mcica_data(mcica_data, filename, sp_sw, sp_lw)

USE def_spectrum, ONLY: StrSpecData
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport

IMPLICIT NONE

TYPE (StrMcica), INTENT(INOUT) :: mcica_data
CHARACTER (LEN=*),  INTENT(IN) :: filename
TYPE (StrSpecData), INTENT(IN), OPTIONAL :: sp_sw
TYPE (StrSpecData), INTENT(IN), OPTIONAL :: sp_lw

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'READ_MCICA_DATA'
INTEGER, PARAMETER :: iu_mcd = 80
INTEGER            :: icode
CHARACTER (LEN=80) :: cmessage
CHARACTER (LEN=80) :: line
INTEGER            :: band, k, subcol

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

OPEN(UNIT=iu_mcd, FILE=filename, IOSTAT=icode, STATUS='OLD')
IF (icode /= 0) THEN
  cmessage='McICA data file could not be opened.'
  CALL ereport(RoutineName, icode, cmessage)
END IF

DO
  ! Read header until data block is found
  READ(iu_mcd, '(a80)', IOSTAT=icode) line
  IF (line(1:5) == '*DATA') EXIT
  IF (icode /= 0) THEN
    cmessage = 'No *DATA block present in McICA data file'
    CALL ereport(RoutineName, icode, cmessage)
  END IF
END DO

DO
  ! Read variables from data block
  READ(iu_mcd, '(a80)', IOSTAT=icode) line
  IF (line(1:4) == '*END') EXIT

  SELECT CASE (line)

  CASE ('tot_subcol_gen','n_subcol_gen')
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%n_subcol_gen
  CASE ('subcol_need_single','n_subcol_req_single')
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%n_subcol_req_single
  CASE ('subcol_need_optimal','n_subcol_req_optimal')
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%n_subcol_req_optimal
  CASE ('ipph')
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%ipph

  CASE ('lw_subcol_reorder_single')
    ALLOCATE(mcica_data%lw_subcol_reorder_single( &
      mcica_data%n_subcol_req_single), STAT=icode)
    IF (icode /= 0) THEN
      cmessage = 'Cannot allocate array: lw_subcol_reorder_single'
      CALL ereport(RoutineName, icode, cmessage)
    END IF
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%lw_subcol_reorder_single

  CASE ('lw_subcol_reorder_optimal')
    ALLOCATE(mcica_data%lw_subcol_reorder_optimal( &
      mcica_data%n_subcol_req_optimal), STAT=icode)
    IF (icode /= 0) THEN
      cmessage = 'Cannot allocate array: lw_subcol_reorder_optimal'
      CALL ereport(RoutineName, icode, cmessage)
    END IF
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%lw_subcol_reorder_optimal

  CASE ('sw_subcol_k')
    IF (PRESENT(sp_sw)) THEN
      ALLOCATE(mcica_data%sw_subcol_k(sp_sw%dim%nd_band, &
                                      sp_sw%dim%nd_k_term), STAT=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: sw_subcol_k'
        CALL ereport(RoutineName, icode, cmessage)
      END IF
      mcica_data%sw_subcol_k=1
      DO
        READ(iu_mcd, '(3i4)', IOSTAT=icode) band, k, subcol
        IF (band == -99) EXIT
        mcica_data%sw_subcol_k(band, k)=subcol
        IF (icode /= 0) THEN
          cmessage = 'Error reading data for sw_subcol_k'
          CALL ereport(RoutineName, icode, cmessage)
        END IF
      END DO
    END IF

  CASE ('lw_subcol_k')
    IF (PRESENT(sp_lw)) THEN
      ALLOCATE(mcica_data%lw_subcol_k(sp_lw%dim%nd_band, &
                                      sp_lw%dim%nd_k_term), STAT=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: lw_subcol_k'
        CALL ereport(RoutineName, icode, cmessage)
      END IF
      mcica_data%lw_subcol_k=1
      DO
        READ(iu_mcd, '(3i4)', IOSTAT=icode) band, k, subcol
        IF (band == -99) EXIT
        mcica_data%lw_subcol_k(band, k)=subcol
        IF (icode /= 0) THEN
          cmessage = 'Error reading data for lw_subcol_k'
          CALL ereport(RoutineName, icode, cmessage)
        END IF
      END DO
    END IF

  CASE ('n1')
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%n1
  CASE ('n2')
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%n2
  CASE ('xcw')
    ALLOCATE(mcica_data%xcw(mcica_data%n1, mcica_data%n2), STAT=icode)
    IF (icode /= 0) THEN
      cmessage = 'Cannot allocate array: xcw'
      CALL ereport(RoutineName, icode, cmessage)
    END IF
    READ(iu_mcd, *, IOSTAT=icode) mcica_data%xcw

  END SELECT

  IF (icode /= 0) THEN
    cmessage = 'Error reading data from McICA data file'
    CALL ereport(RoutineName, icode, cmessage)
  END IF
END DO

CLOSE(iu_mcd)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_mcica_data
!------------------------------------------------------------------------------
SUBROUTINE deallocate_mcica_data(mcica_data)

IMPLICIT NONE

TYPE (StrMcica), INTENT(INOUT) :: mcica_data

IF (ALLOCATED(mcica_data%xcw)) &
   DEALLOCATE(mcica_data%xcw)
IF (ALLOCATED(mcica_data%lw_subcol_k)) &
   DEALLOCATE(mcica_data%lw_subcol_k)
IF (ALLOCATED(mcica_data%sw_subcol_k)) &
   DEALLOCATE(mcica_data%sw_subcol_k)
IF (ALLOCATED(mcica_data%lw_subcol_reorder_optimal)) &
   DEALLOCATE(mcica_data%lw_subcol_reorder_optimal)
IF (ALLOCATED(mcica_data%lw_subcol_reorder_single)) &
   DEALLOCATE(mcica_data%lw_subcol_reorder_single)

END SUBROUTINE deallocate_mcica_data
!------------------------------------------------------------------------------

END MODULE def_mcica

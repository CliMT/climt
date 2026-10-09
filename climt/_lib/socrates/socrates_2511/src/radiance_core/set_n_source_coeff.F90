! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE set_n_source_coeff_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_N_SOURCE_COEFF_MOD'
CONTAINS

! Function to set number of source coefficients.
!
! Method:
!   The two-stream approximation is examined and the number
!   of coefficients is set accordingly.
!- ---------------------------------------------------------------------
FUNCTION set_n_source_coeff(isolir, l_ir_source_quad)

  USE rad_pcf,  ONLY: ip_solar
  USE yomhook,  ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


  INTEGER, INTENT(IN) :: isolir
!   Spectral region
  LOGICAL, INTENT(IN) :: l_ir_source_quad
!   Flag for quadratic infra-red source

  INTEGER :: set_n_source_coeff
!   Returned number of source coefficients

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_N_SOURCE_COEFF'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (isolir == ip_solar) THEN
    set_n_source_coeff=2
  ELSE
    IF (l_ir_source_quad) THEN
      set_n_source_coeff=2
    ELSE
      set_n_source_coeff=1
    END IF
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION set_n_source_coeff
END MODULE set_n_source_coeff_mod

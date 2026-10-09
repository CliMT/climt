! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define the arrays used in the CKD continuum model.
!
MODULE ckd_extern_data
!
! Description:
!   This module defines arrays used by the CKD continuum model.
!   These arrays would usually be defined with the code, but
!   are too big for compilation on all machines, so a 
!   process of reading in is adopted.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd
!
!
  IMPLICIT NONE
!
!
!
  TYPE StrCKDFld
!
   INTEGER :: n_freq
!    Number of frequencies in the table
   REAL  (RealK) :: table_start
!    First frequency in table
   REAL  (RealK) :: table_end
!    Last frequency in table
   REAL  (RealK) :: table_inc
!    Increment along table
!
   REAL  (RealK), Pointer, Dimension(:) :: c => NULL()
!    Continuum coefficients
   REAL  (RealK), Pointer, Dimension(:) :: texp => NULL()
!    Temperature exponent
!
  END TYPE StrCKDFld
!
!
!
  TYPE (StrCKDFld) :: c_foreign_h2o_296
!   Foreign-broadened coefficients for the CKD continuum
  TYPE (StrCKDFld) :: c_self_h2o_260
!   Self-broadened coefficients for the CKD continuum at 260 K
  TYPE (StrCKDFld) :: c_self_h2o_296
!   Self-broadened coefficients for the CKD continuum at 296 K
!
END MODULE ckd_extern_data

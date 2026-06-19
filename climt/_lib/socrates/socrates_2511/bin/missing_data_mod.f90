! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

! Data module containing missing data indicies

! Migrated from include file c_mdi.h

! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!

MODULE missing_data_mod

IMPLICIT NONE

     ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768

END MODULE missing_data_mod

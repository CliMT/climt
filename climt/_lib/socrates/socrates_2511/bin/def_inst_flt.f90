! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to declare a structure for instrumental filter fiunctions.

MODULE def_inst_flt

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE


  TYPE StrFiltResp

    CHARACTER (LEN=32) :: satellite = "                                "
!     Name of the satellite
    CHARACTER (LEN=32) :: instrument = "                                "
!     Name of the instrument
    CHARACTER (LEN=32) :: channel = "                                "
!     Name of the channel

    INTEGER :: n_pts
!     Number of defining the instrumental response

    REAL  (RealK), Pointer :: wavenumber(:)
!     Wavenumbers where the filter response is given
    REAL  (RealK), Pointer :: response(:)
!     Values of the reponse function at these wavenumbers
    REAL  (RealK), Pointer :: d2_response(:)
!     Second derivative of the reponse function at these wavenumbers

  END TYPE StrFiltResp

END MODULE def_inst_flt

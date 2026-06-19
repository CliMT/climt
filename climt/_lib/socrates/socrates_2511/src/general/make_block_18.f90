! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 18.
!
! Description:
!   Block 18 contains the indexing numbers for the generalised continua
!   active in each band. As this information is automatically generated
!   when reading in the continuum data, there is no need for the user
!   to specify this when creating the spectral file. This block is
!   therefore initially created without user input with zero continua
!   in each band initially.
!
!------------------------------------------------------------------------------
SUBROUTINE make_block_18(Spectrum, l_interactive, ierr)

  USE def_spectrum

  IMPLICIT NONE


! Dummy arguments
  LOGICAL, INTENT(IN) :: l_interactive
!   Flag for interactive operation
  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral data
  INTEGER, INTENT(INOUT) :: ierr
!   Error flag

! Local variables
  TYPE (StrSpecDim), Pointer :: SpDim
!   Pointer to dimensions within the spectrum

  SpDim   => Spectrum%Dim
  ALLOCATE(Spectrum%ContGen%n_band_cont(SpDim%nd_band))
  ALLOCATE(Spectrum%ContGen%index_cont(SpDim%nd_cont, SpDim%nd_band))
  ALLOCATE(Spectrum%ContGen%l_cont_major(SpDim%nd_band))

  Spectrum%ContGen%n_band_cont = 0
  Spectrum%ContGen%index_cont = 0
  Spectrum%ContGen%l_cont_major = .FALSE.

END SUBROUTINE make_block_18

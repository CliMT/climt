! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module holds tools to generate appropriate masks for COSP2 diagnostics
!
!------------------------------------------------------------------------------

module cosp2_diagnostics_mod

use cosp_kinds, only: wp
use realtype_rd, only: RealExt

implicit none

contains

elemental subroutine create_mask(mdi,x,y)
implicit none
real(wp),intent(in) :: mdi
real(wp),intent(inout) :: x
real(RealExt),intent(inout) :: y

if (x == mdi) then
  x = 0.0_wp
  y = 0.0_RealExt
else
  y = 1.0_RealExt
end if
end subroutine create_mask

end module cosp2_diagnostics_mod

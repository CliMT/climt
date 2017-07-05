module phy_finalize
! stub that provides do nothing finalize_phy routine.

 use phy_data, only: destroy_phydata
 implicit none
 private
 public :: finalize_phy

 contains

 subroutine finalize_phy()
    call destroy_phydata()
    return
 end subroutine finalize_phy


end module phy_finalize

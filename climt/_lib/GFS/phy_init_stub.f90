module phy_init
! stub that provides do nother init_phy routine.

 use phy_data, only: init_phydata
 implicit none
 private
 public :: init_phy

 contains

 subroutine init_phy() bind(c,name='gfs_init_physics')
 end subroutine init_phy

end module phy_init

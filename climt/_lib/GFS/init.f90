module init_mod
! initialization module (provides init driver subroutine)
private
!public :: init
contains

!subroutine init() bind(c,name='dycoreInit')
!!  use params, only: read_namelist
!  use dyn_init ,only: init_dyn
!  use phy_init ,only: init_phy
!
!  implicit none
!  print *, 'Im here!!!'
!!  call read_namelist ()
!  call init_dyn ()
!  call init_phy ()
!  
!  return
!end subroutine init
end module init_mod

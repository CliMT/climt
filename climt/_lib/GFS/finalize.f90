module finalize_mod
! finalization module (provides finalize driver subroutine)
private
public :: finalize
contains

subroutine finalize() bind(c,name='gfs_finalise');
  use dyn_finalize ,only: finalize_dyn
  use phy_finalize ,only: finalize_phy
  !JOY this was not being called at all
  use shtns, only: shtns_destroy

  implicit none

  call finalize_dyn ()
  call finalize_phy ()
  !JOY this was not being called at all
  call shtns_destroy()
  
  return
end subroutine finalize
end module finalize_mod

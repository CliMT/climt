program driver

  ! driver for gfs dynamical core

  use init_mod, only: init
  use run_mod, only: run
  use finalize_mod, only: finalize 

  implicit none

  ! read in namelist, initial conditions, 
  ! initialize physics and dynamics, allocate arrays.
  call init() 
  ! run model, write out at specified intervals.
  call run()
  ! deallocate arrays.
  call finalize()

end program driver

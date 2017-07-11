module dyn_finalize
! finalize dynamics (provides finalize_dyn subroutine)

 use shtns, only: shtns_destroy
 use spectral_data, only: destroy_specdata
 use pressure_data, only: destroy_pressdata
 use grid_data, only: destroy_griddata
! JOY removing stochastic routines 
! use stoch_data, only: destroy_stochdata
 use semimp_data, only: destroy_semimpdata
 use params, only: explicit, iau
! JOY remove Incremental Analysis Unit (IAU) forcing                          
! use iau_module, only: destroy_iau

 implicit none
 private
 public :: finalize_dyn

 contains

 subroutine finalize_dyn()
    ! call routines to deallocate arrays
    call destroy_specdata()
    call destroy_griddata()
    call destroy_pressdata()
    if (.not. explicit) call destroy_semimpdata()
! JOY removing stochastic routines 
!    call destroy_stochdata()
    call shtns_destroy()
! JOY remove Incremental Analysis Unit (IAU) forcing                          
!    if (iau) call destroy_iau()
 end subroutine finalize_dyn


end module dyn_finalize

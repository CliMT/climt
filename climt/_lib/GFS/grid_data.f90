module grid_data
! gaussian grid model state variable data.
! public subroutines:
! init_griddata: allocate arrays.
! destroy_griddata: deallocate arrays.
 use kinds, only: r_kind
 use params, only: nlons,nlats,nlevs,ntrac
 use iso_c_binding, only: c_double
 implicit none
 private
 public :: init_griddata, destroy_griddata
 real(r_kind), pointer, public, dimension(:,:,:) :: ug,vg,vrtg,divg,&
 virtempg,etadot,dlnpdtg
 real(r_kind), pointer, public, dimension(:,:,:,:) :: tracerg
! (nlons,nlons,nlevs) arrays (bottom to top unless otherwise noted)
! they are transformed to the grid from spectral space in subroutine
! getdyntend in module dyn_run.
! ug: zonal wind
! vg: meridional wind
! vrtg: vorticity
! divg: divergence
! tracerg: tracers (first one is specific humidity)
! virtempg: virtual temperature
! the following last two are computed in subroutine omega from module dyn_run
! (called by dyntend):
! etadot: vertical motion in hybrid sigma-pressure (top to bottom)
! dlnpdtg: d(lnp)/dt = omega/p (compute in subroutine omega from module dyn_run)
! dlnpsdt: local tendency of ln(ps)
 real(r_kind), pointer, public, dimension(:,:) :: &
 dlnpsdt,lnpsg,dlnpdst,dphisdx,dphisdy,phis

 contains

     subroutine assign_grid_arrays(py_ug, py_vg, py_virt_tempg, py_lnpsg, &
             py_tracerg, py_vortg, py_divg) bind(c, name='gfs_assign_grid_arrays')

         real(c_double), target, intent(inout):: py_ug(nlons, nlats, nlevs)
         real(c_double), target, intent(inout):: py_vg(nlons, nlats, nlevs)
         real(c_double), target, intent(inout):: py_virt_tempg(nlons, nlats, nlevs)
         real(c_double), target, intent(inout):: py_vortg(nlons, nlats, nlevs)
         real(c_double), target, intent(inout):: py_divg(nlons, nlats, nlevs)

         real(c_double), target, intent(inout):: py_lnpsg(nlons, nlats)
         real(c_double), target, intent(inout):: py_tracerg(nlons, nlats, nlevs, ntrac)

         ug => py_ug
         vg => py_vg
         vrtg => py_vortg
         divg => py_divg
         virtempg => py_virt_tempg
         tracerg => py_tracerg
         lnpsg => py_lnpsg

     end subroutine

     subroutine init_griddata(pyDlnpdtg, pyEtaDotg, pyPhis, &
             pyDPhisdx, pyDPhisdy, pyDlnpsdt) bind(c,name='gfs_initialise_grid_arrays')

         real(c_double), target, intent(inout) :: pyDlnpdtg(nlons, nlats, nlevs)

         real(c_double), target, intent(inout) :: pyEtaDotg(nlons,nlats,nlevs+1)

         real(c_double), target, intent(inout) :: pyPhis(nlons,nlats), pyDPhisdx(nlons,nlats), &
             pyDPhisdy(nlons,nlats), pyDlnpsdt(nlons,nlats)

         !print *, 'Initialising grid data'

         !print *,  shape(pyDlnpdtg);
         dlnpdtg => pyDlnpdtg
         etadot => pyEtaDotg
         phis => pyPhis
         dphisdy => pyDPhisdy
         dphisdx => pyDPhisdx
         dlnpsdt => pyDlnpsdt
         !print *, shape(dlnpdtg);
         !print *, ' Done Initialising grid data'
     end subroutine init_griddata
     
     subroutine destroy_griddata()
         !print *, 'Deallocating grid data'
         nullify(ug,vg,vrtg,divg,virtempg,dlnpdtg,etadot)
         nullify(tracerg)
         nullify(lnpsg,phis,dphisdx,dphisdy,dlnpsdt)
     end subroutine destroy_griddata

 end module grid_data

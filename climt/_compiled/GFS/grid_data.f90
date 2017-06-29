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

 !public :: init_griddata, destroy_griddata
 !real(r_kind), allocatable, public, dimension(:,:,:) :: ug,vg,vrtg,divg,&
 !virtempg,etadot,dlnpdtg
 !real(r_kind), allocatable, public, dimension(:,:,:,:) :: tracerg

 !real(r_kind), allocatable, public, dimension(:,:) :: &
 !dlnpsdt,lnpsg,dlnpdst,dphisdx,dphisdy,phis
 contains
 subroutine init_griddata(pyUg, pyVg, pyVrtg, pyDivg, pyVirtTempg, &
         pyTracerg, pyDlnpdtg, pyEtaDotg, pyLnPsg, pyPhis, &
         pyDPhisdx, pyDPhisdy, pyDlnpsdt) bind(c,name='initialiseGridArrays')

    real(c_double), target, intent(inout) :: pyUg(nlons,nlats,nlevs), &
    pyVg(nlons,nlats,nlevs),pyVrtg(nlons,nlats,nlevs),pyDivg(nlons,nlats,nlevs),&
    pyVirtTempg(nlons,nlats,nlevs),pyDlnpdtg(nlons,nlats,nlevs)

    real(c_double), target, intent(inout) :: pyEtaDotg(nlons,nlats,nlevs+1),&
    pyTracerg(nlons,nlats,nlevs,ntrac)

    real(c_double), target, intent(inout) :: pyLnPsg(nlons,nlats), &
    pyPhis(nlons,nlats), pyDPhisdx(nlons,nlats), pyDPhisdy(nlons,nlats), &
    pyDlnpsdt(nlons,nlats)

    print *, 'Initialising grid data'

    print *, shape(pyUg), shape(pyDlnpdtg);
    ug => pyUg
    vg => pyVg
    vrtg => pyVrtg
    divg => pyDivg
    virtempg => pyVirtTempg
    tracerg => pyTracerg
    dlnpdtg => pyDlnpdtg
    etadot => pyEtaDotg
    lnpsg => pyLnPsg
    phis => pyPhis
    dphisdy => pyDPhisdy
    dphisdx => pyDPhisdx
    dlnpsdt => pyDlnpsdt
    print *, shape(ug), shape(dlnpdtg);
 end subroutine init_griddata
 subroutine destroy_griddata()
    print *, 'Deallocating grid data'
    nullify(ug,vg,vrtg,divg,virtempg,dlnpdtg,etadot)
    nullify(tracerg)
    nullify(lnpsg,phis,dphisdx,dphisdy,dlnpsdt)
 end subroutine destroy_griddata

end module grid_data

!    allocate(ug(nlons,nlats,nlevs))
!    allocate(vg(nlons,nlats,nlevs))
!    allocate(vrtg(nlons,nlats,nlevs))
!    allocate(divg(nlons,nlats,nlevs))
!    allocate(virtempg(nlons,nlats,nlevs))
!    allocate(tracerg(nlons,nlats,nlevs,ntrac))
!    allocate(dlnpdtg(nlons,nlats,nlevs))
!    allocate(etadot(nlons,nlats,nlevs+1))
!    allocate(lnpsg(nlons,nlats))
!    allocate(phis(nlons,nlats))
!    allocate(dphisdy(nlons,nlats))
!    allocate(dphisdx(nlons,nlats))
!    allocate(dlnpsdt(nlons,nlats))
!    deallocate(ug,vg,vrtg,divg,virtempg,dlnpdtg,etadot)
!    deallocate(tracerg)
!    deallocate(lnpsg,phis,dphisdx,dphisdy,dlnpsdt)
!

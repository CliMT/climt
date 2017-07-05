module pressure_data
! data for model level pressures and pressure-related
! quantities.
! init_pressdata: allocates arrays.
! calc_pressdata: computes pressure related vars given ln(psg).
! (called by getdyntend in module dyn_run).
! destroy_pressdata: deallocate arrays.
 use params, only: nlons,nlats,nlevs
 use kinds, only: r_kind
 use physcons, only: con_rd,con_cp,rk => con_rocp
 implicit none
 private
 public :: ak,bk,ck,dbk,bkl,psg,pk,alfa,rlnp,dpk,&
  si,sl,prs,init_pressdata, calc_pressdata, destroy_pressdata
! ak,bk: definition of hybrid sigma-pressure coordinates
! si,sl: equivalent values of sigma coordinate at interfaces, mid-levels.
! dbk(k) = bk(k+1)-bk(k)
! bkl(k) = 0.5*(bk(k+1)+bk(k))
! ck(k)  = ak(k+1)*bk(k)-ak(k)*bk(k+1)
! above 1-d arrays computed in init_dyn, defined from top to bottom
! of model.
! the following are definied by call to calc_pressdata:
! psg: surface pressure (Pa)
! pk: model interface pressures = ak(k) + bk(k)*psg(:,:) (k=1,nlevs+1)
! dpk = pk(:,:,k+1)-pk(:,:,k)
! prs: model layer pressures 
! rlnp = log(pk(:,:,k+1)/pk(:,:,k))
! alfa = 1.-(pk(:,:,k)/dpk(:,:,k))*rlnp(:,:,k) (for k=1 alfa(:,:,1)=log(2))
! all of the above arrays are defined top to bottom, except
! prs which is bottom to top.
 real(r_kind), dimension(:), allocatable :: ak,bk,ck,dbk,si,sl,bkl
! real(r_kind), dimension(:,:), allocatable :: psg
 real(r_kind), dimension(:,:,:), allocatable :: alfa,rlnp,dpk
!JOY converted prs, pk and psg to pointer, memory allocated by python
 real(r_kind), dimension(:,:,:), pointer :: prs => null()
 real(r_kind), dimension(:,:,:), pointer :: pk => null()
 real(r_kind), dimension(:,:), pointer :: psg => null()
 
 contains

 subroutine assignPressureArrays(pySurfPressure, pyPressGrid, pyInterfacePressure) bind(c, name='initialisePressureArrays')

    real(r_kind), intent(inout), target :: pyPressGrid(nlons,nlats,nlevs)
    real(r_kind), intent(inout), target :: pyInterfacePressure(nlons,nlats,nlevs)
    real(r_kind), intent(inout), target :: pySurfPressure(nlons,nlats)

    print *, 'Initialising pressure arrays'

    prs => pyPressGrid
    psg => pySurfPressure
    pk => pyInterfacePressure

end subroutine assignPressureArrays

subroutine init_pressdata()

    if (allocated(ak)) then
        deallocate(ak)
    endif
    allocate(ak(nlevs+1))

    if (allocated(bk)) then
        deallocate(bk)
    endif
    allocate(bk(nlevs+1))

    if (allocated(si)) then
        deallocate(si)
    endif
    allocate(si(nlevs+1))

    if (allocated(sl)) then
        deallocate(sl)
    endif
    allocate(sl(nlevs))

    if (allocated(ck)) then
        deallocate(ck)
    endif
    allocate(ck(nlevs))

    if (allocated(dbk)) then
        deallocate(dbk)
    endif
    allocate(dbk(nlevs))

    if (allocated(bkl)) then
        deallocate(bkl)
    endif
    allocate(bkl(nlevs))
    !JOY allocated in python
    !allocate(psg(nlons,nlats))
    if (allocated(alfa)) then
        deallocate(alfa)
    endif
    allocate(alfa(nlons,nlats,nlevs))

    if (allocated(rlnp)) then
        deallocate(rlnp)
    endif
    allocate(rlnp(nlons,nlats,nlevs))

    if (allocated(dpk)) then
        deallocate(dpk)
    endif
    allocate(dpk(nlons,nlats,nlevs))
    !JOY allocated in python
   ! allocate(pk(nlons,nlats,nlevs+1))
   ! allocate(prs(nlons,nlats,nlevs))
 end subroutine init_pressdata

 subroutine destroy_pressdata()
    print *, 'Deallocating pressure data'
    deallocate(ak,bk,ck,dbk,si,sl,bkl)
    deallocate(alfa,rlnp,dpk)
    !JOY these are just pointers, no deallocation needed
    nullify(prs,psg,pk)
 end subroutine destroy_pressdata

 subroutine updatePressure() bind(c, name='gfsCalcPressure')
    use grid_data, only: lnpsg
    call calc_pressdata(lnpsg)

 end subroutine updatePressure

 subroutine calc_pressdata(lnpsg) 
    real(r_kind),  intent(in) :: lnpsg(nlons,nlats)
! update pressure related variables using latest estimate of lnps
    integer k
! compute sigma coordinate quantities (bottom to top).
    do k=1,nlevs+1
       si(nlevs+2-k)= ak(k)/101300.0+bk(k) ! si are now sigmas
    enddo
    do k=1,nlevs
       sl(k) = 0.5*(si(k)+si(k+1))
       !sl(k) = ((si(k)**(rk+1.) - si(k+1)**(rk+1.))/&
       !        ((rk+1.)*(si(k)-si(k+1))))**(1./rk)
    enddo
    ! surface press from log(ps)
    psg = exp(lnpsg)
    ! interface pressure (psg is ps)
    ! ak,bk go top to bottom, so does pk
!$omp parallel
!$omp do private(k)
    do k=1,nlevs+1
      pk(:,:,k)=ak(k) + bk(k)*psg(:,:)
    enddo
!$omp end do
!$omp do private(k)
    do k=1,nlevs
       ! layer pressure thickness
       dpk(:,:,k)=    pk(:,:,k+1) - pk(:,:,k)
       ! sela's layer pressure from hyb2press.f
       ! (goes from bottom to top, unlike pk)
       prs(:,:,nlevs-k+1) = ((pk(:,:,k+1)**rk*pk(:,:,k+1) - pk(:,:,k)**rk*pk(:,:,k))/&
                    ((rk+1.)*dpk(:,:,k))) ** (1./rk)
    enddo
!$omp end do
!$omp end parallel  
    alfa(:,:,1)=log(2.)
    rlnp(:,:,1)=99999.99 !doesn't matter, should never be used.
!$omp parallel do private(k)
    do k=2,nlevs
      rlnp(:,:,k)= log( pk(:,:,k+1)/pk(:,:,k) )
      alfa(:,:,k)= 1.-( pk(:,:,k)/dpk(:,:,k) )*rlnp(:,:,k)
    enddo
!$omp end parallel do 
 end subroutine calc_pressdata

end module pressure_data

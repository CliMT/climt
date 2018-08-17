module pressure_data
! data for model level pressures and pressure-related
! quantities.
! init_pressdata: allocates arrays.
! calc_pressdata: computes pressure related vars given ln(psg).
! (called by getdyntend in module dyn_run).
! destroy_pressdata: deallocate arrays.
 use params, only: nlons,nlats,nlevs,toa_pressure
 use kinds, only: r_kind
 use physcons, only: con_rd,con_cp,rk => con_rocp
 use iso_c_binding
 implicit none
 private
 public :: ak,bk,ck,dbk,bkl,psg,pk,alfa,rlnp,dpk,&
  si,sl,prs, calc_pressdata, destroy_pressdata
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
 real(c_double), dimension(:), pointer :: ak => null()
 real(c_double), dimension(:), pointer :: bk => null()
 real(c_double), dimension(:), pointer :: ck => null()
 real(c_double), dimension(:), pointer :: dbk => null()
 real(c_double), dimension(:), pointer :: si => null()
 real(c_double), dimension(:), pointer :: sl => null()
 real(c_double), dimension(:), pointer :: bkl => null()
! real(c_double), dimension(:,:), allocatable :: psg
 real(c_double), dimension(:,:,:), pointer :: alfa => null()
 real(c_double), dimension(:,:,:), pointer :: rlnp => null()
 real(c_double), dimension(:,:,:), pointer :: dpk => null()
!JOY converted prs, pk and psg to pointer, memory allocated by python
 real(c_double), dimension(:,:,:), pointer :: prs => null()
 real(c_double), dimension(:,:,:), pointer :: pk => null()
 real(c_double), dimension(:,:), pointer :: psg => null()

 contains

     subroutine assign_coordinate_arrays(py_ak, py_bk, py_ck, py_dbk,&
             py_si, py_sl, py_bkl, py_alfa, py_rlnp, py_dpk)bind(c, name='gfs_init_coordinate_arrays')

         real(c_double), intent(inout), target:: py_ak(nlevs+1), py_bk(nlevs+1), py_ck(nlevs),&
             py_dbk(nlevs), py_si(nlevs+1), py_sl(nlevs), py_bkl(nlevs)
         real(c_double), intent(inout), target:: py_alfa(nlons, nlats, nlevs), &
             py_rlnp(nlons, nlats, nlevs), py_dpk(nlons, nlats, nlevs)

         ak => py_ak
         bk => py_bk
         ck => py_ck
         dbk => py_dbk
         si => py_si
         sl => py_sl
         bkl => py_bkl

         alfa => py_alfa
         rlnp => py_rlnp
         dpk => py_dpk

     end subroutine

     subroutine assignPressureArrays(pySurfPressure, pyPressGrid, pyInterfacePressure) bind(c, name='gfs_assign_pressure_arrays')

         real(c_double), intent(inout), target :: pyPressGrid(nlons,nlats,nlevs)
         real(c_double), intent(inout), target :: pyInterfacePressure(nlons,nlats,nlevs)
         real(c_double), intent(inout), target :: pySurfPressure(nlons,nlats)

         prs => pyPressGrid
         pk => pyInterfacePressure
         psg => pySurfPressure

     end subroutine assignPressureArrays


     subroutine destroy_pressdata()
         !print *, 'Deallocating pressure data'
         nullify(ak,bk,ck,dbk,si,sl,bkl)
         nullify(alfa,rlnp,dpk)
         !JOY these are just pointers, no deallocation needed
         nullify(prs,psg,pk)
     end subroutine destroy_pressdata

     subroutine updatePressure() bind(c, name='gfs_calculate_pressure')
         use grid_data, only: lnpsg

         call calc_pressdata(lnpsg)

     end subroutine updatePressure

     subroutine calc_pressdata(lnpsg)
         real(r_kind),  intent(in) :: lnpsg(nlons,nlats)
         ! update pressure related variables using latest estimate of lnps
         integer k
         ! surface press from log(ps)
         psg = exp(lnpsg)
         ! interface pressure (psg is ps)
         ! ak,bk go top to bottom, so does pk
         !$omp parallel
         !$omp do private(k)
         do k=1,nlevs+1
         pk(:,:,k)=ak(k) + bk(k)*(psg(:,:) - toa_pressure)
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

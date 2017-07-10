 module phy_run
! compute physics tendencies for simplifed held-suarez forcing
! getphytend: compute physics tendencies.
! Public subroutines:
! getphytend: compute newtonian damping of temperature, drag terms
! in vort and div eqns.

 use params, only: nlevs,nlons,nlats,ndimspec,ntrac,pdryini
 use kinds, only: r_kind, r_double
 use shtns, only: grdtospec, lats, areawts
 !use grid_data, only: vrtg,divg,virtempg
 !use pressure_data, only:  prs,psg
 !use spectral_data, only: lnpsspec

 use phy_data, only: spectral_vort_tend, spectral_div_tend, spectral_virtemp_tend, spectral_tracer_tend, spectral_lnps_tend
 use,intrinsic  :: iso_c_binding


 implicit none
 private


 ! Define interface of call-back routine.
   !    ABSTRACT INTERFACE
   !      SUBROUTINE callback (x)
 !         USE, INTRINSIC :: ISO_C_BINDING
 !          REAL(KIND=C_DOUBLE), INTENT(IN), VALUE :: x
!         END SUBROUTINE callback
!       END INTERFACE
 ! Define interface of call-back routine.
       ABSTRACT INTERFACE
         SUBROUTINE callback ()
           USE, INTRINSIC :: ISO_C_BINDING
         END SUBROUTINE callback
       END INTERFACE

 public :: getphytend, gfsRegisterPhysicsCallback, calculate_heldsuarez, set_tendencies

 procedure(callback), pointer :: pyPhysicsCallback => null()

 contains

 subroutine gfsRegisterPhysicsCallback(callback) bind(c,name='gfs_register_physics_callback')

    type(c_funptr), intent(in), value :: callback

    print *, 'Registering callback'
    print *, callback
    call c_f_procpointer(callback, pyPhysicsCallback)

    call pyPhysicsCallback()

 end subroutine gfsRegisterPhysicsCallback

 subroutine getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt,t,dt)
     complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
         dvrtspecdt,ddivspecdt,dvirtempspecdt
     complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
         dtracerspecdt
     real(r_double), intent(in) :: t,dt
     complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
     integer k, trac


     !$omp parallel do private(k)
     do k=1,nlevs
     dvrtspecdt(:,k) = spectral_vort_tend(:,k)
     ddivspecdt(:,k) = spectral_div_tend(:,k)
     dvirtempspecdt(:,k) = spectral_virtemp_tend(:,k)
     do trac=1,ntrac
         dtracerspecdt(:,k,trac) = spectral_tracer_tend(:,k,trac)
     enddo
     enddo
     !$omp end parallel do 
     dlnpsspecdt(:) = spectral_lnps_tend(:)


 end subroutine getphytend

 subroutine calculate_heldsuarez(vrtg, divg, virtempg, prs, psg, tracerg, t, dt) bind(c, name='gfs_calculate_tendencies')

     real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: vrtg, divg, virtempg, prs
     real(r_kind), intent(in), dimension(nlons,nlats) :: psg 
     real(r_kind), intent(in), dimension(nlons,nlats,nlevs,ntrac) :: tracerg
     real(r_kind) p0,sigbot,tempstrat,delthz,deltmp,&
                kdrag,krada,kradb,pdry,pcorr
     real(r_kind), dimension(:,:,:),allocatable :: blprof,radequiltemp,forcingg
     complex(r_kind), dimension(:,:),allocatable :: forcingspec
     real(r_double), intent(in) :: t,dt
     integer k

     allocate(blprof(nlons,nlats,nlevs),radequiltemp(nlons,nlats,nlevs))
     allocate(forcingg(nlons,nlats,nlevs),forcingspec(ndimspec,nlevs))


     sigbot = 0.7
     delthz = 10.
     tempstrat = 200.
     kdrag = 1./(1.*86400.)
     krada = 1./(40.*86400.)
     kradb = 1./(4.*86400. )
     p0 = 1.e5
     deltmp = 60.
     !$omp parallel do private(k)
     do k=1,nlevs
     blprof(:,:,k) = prs(:,:,k)/psg
     radequiltemp(:,:,k) = (prs(:,:,k)/p0)**(2./7.)*&
         (315.-deltmp*sin(lats)**2-delthz*log(prs(:,:,k)/p0)*cos(lats)**2)
     enddo
     !$omp end parallel do 
     blprof = (blprof-sigbot)/(1.-sigbot)
     where (blprof < 0) 
         blprof = 0
     endwhere
     where (radequiltemp < tempstrat)
         radequiltemp=tempstrat
     endwhere
     !$omp parallel do private(k)
     do k=1,nlevs
     forcingg(:,:,k)=(krada+(kradb-krada)*blprof(:,:,k)*cos(lats)**4)*&
         (radequiltemp(:,:,k)-virtempg(:,:,k))
     call grdtospec(forcingg(:,:,k), spectral_virtemp_tend(:,k))
     !dvirtempspecdt(:,k) = forcingspec(:,k)
     forcingg(:,:,k) = -(blprof(:,:,k)*kdrag)*vrtg(:,:,k)
     call grdtospec(forcingg(:,:,k), spectral_vort_tend(:,k))
     !dvrtspecdt(:,k) = forcingspec(:,k)
     forcingg(:,:,k) = -(blprof(:,:,k)*kdrag)*divg(:,:,k)
     call grdtospec(forcingg(:,:,k), spectral_div_tend(:,k))
     !ddivspecdt(:,k) = forcingspec(:,k)
     enddo
     !$omp end parallel do 
     spectral_lnps_tend = 0.

     !print *, 'Maximum tendency in divergence: '
     !print *, maxval(forcingg)

     deallocate(blprof,radequiltemp)
     deallocate(forcingg,forcingspec)

     return

end subroutine calculate_heldsuarez

subroutine set_tendencies(dvrtgdt,ddivgdt,dvirtempgdt,dlnpsgdt,dtracergdt) bind(c, name='gfs_set_tendencies')

    real(r_kind), intent(in), dimension(nlons,nlats,nlevs) :: &
        dvrtgdt,ddivgdt,dvirtempgdt
    real(r_kind), intent(in), dimension(nlons,nlats,nlevs,ntrac) :: &
        dtracergdt
    real(r_kind), intent(in), dimension(nlons,nlats) :: dlnpsgdt
    !complex(r_kind), dimension(:,:),allocatable :: forcingspec
    integer k,trac

    !allocate(forcingspec(ndimspec,nlevs))

    !print *, 'In phy_run: '
    !print *, maxval(abs(dvrtgdt(:,:,9))), maxval(abs(ddivgdt)), maxval(abs(dvirtempgdt(:,:,9)))
    !$omp parallel do private(k)
    do k=1,nlevs
    call grdtospec(dvrtgdt(:,:,k), spectral_vort_tend(:,k))
    call grdtospec(ddivgdt(:,:,k), spectral_div_tend(:,k))
    call grdtospec(dvirtempgdt(:,:,k), spectral_virtemp_tend(:,k))
    !spectral_vort_tend(:,k) = dvrtspecdt(:,k)
    !spectral_div_tend(:,k) = ddivspecdt(:,k) 
    !spectral_virtemp_tend(:,k) = dvirtempspecdt(:,k)
    do trac=1,ntrac
        call grdtospec(dtracergdt(:,:,k,ntrac), spectral_tracer_tend(:,k,ntrac))
    enddo
    enddo
    !$omp end parallel do 
    call grdtospec(dlnpsgdt(:,:), spectral_lnps_tend(:))
    !spectral_lnps_tend(:) = dlnpsspecdt(:)

    !deallocate(spectral_lnps_tend)
    !deallocate(spectral_vort_tend,spectral_div_tend,spectral_virtemp_tend)
    !deallocate(forcingspec)
end subroutine set_tendencies


! subroutine getphytend_old(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt,t,dt)
!   ! compute physics tendencies for held-suarez test case.
!   ! http://dx.doi.org/10.1175/1520-0477(1994)075%3C1825:APFTIO%3E2.0.CO;2
!   complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
!   dvrtspecdt,ddivspecdt,dvirtempspecdt
!   complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
!   dtracerspecdt
!   real(r_double), intent(in) :: t,dt
!   complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
!   real(r_kind) p0,sigbot,tempstrat,delthz,deltmp,&
!                kdrag,krada,kradb,pdry,pcorr
!   real(r_kind), dimension(:,:,:),allocatable :: blprof,radequiltemp,forcingg
!   complex(r_kind), dimension(:,:),allocatable :: forcingspec
!   real(r_kind), allocatable, dimension(:,:) :: wrkg
!   complex(r_kind), allocatable, dimension(:) :: wrkspec
!   integer k
!
!   allocate(blprof(nlons,nlats,nlevs),radequiltemp(nlons,nlats,nlevs))
!   allocate(forcingg(nlons,nlats,nlevs),forcingspec(ndimspec,nlevs))
!
!   sigbot = 0.7
!   delthz = 10.
!   tempstrat = 200.
!   kdrag = 1./(1.*86400.)
!   krada = 1./(40.*86400.)
!   kradb = 1./(4.*86400. )
!   p0 = 1.e5
!   deltmp = 60.
!!$omp parallel do private(k)
!   do k=1,nlevs
!      blprof(:,:,k) = prs(:,:,k)/psg
!      radequiltemp(:,:,k) = (prs(:,:,k)/p0)**(2./7.)*&
!                     (315.-deltmp*sin(lats)**2-delthz*log(prs(:,:,k)/p0)*cos(lats)**2)
!   enddo
!!$omp end parallel do 
!   blprof = (blprof-sigbot)/(1.-sigbot)
!   where (blprof < 0) 
!     blprof = 0
!   endwhere
!   where (radequiltemp < tempstrat)
!      radequiltemp=tempstrat
!   endwhere
!!$omp parallel do private(k)
!   do k=1,nlevs
!      forcingg(:,:,k)=(krada+(kradb-krada)*blprof(:,:,k)*cos(lats)**4)*&
!                      (radequiltemp(:,:,k)-virtempg(:,:,k))
!      call grdtospec(forcingg(:,:,k), forcingspec(:,k))
!      dvirtempspecdt(:,k) = forcingspec(:,k)
!      forcingg(:,:,k) = -(blprof(:,:,k)*kdrag)*vrtg(:,:,k)
!      call grdtospec(forcingg(:,:,k), forcingspec(:,k))
!      dvrtspecdt(:,k) = forcingspec(:,k)
!      forcingg(:,:,k) = -(blprof(:,:,k)*kdrag)*divg(:,:,k)
!      call grdtospec(forcingg(:,:,k), forcingspec(:,k))
!      ddivspecdt(:,k) = forcingspec(:,k)
!   enddo
!!$omp end parallel do 
!   dlnpsspecdt = 0.
!
!   deallocate(blprof,radequiltemp)
!   deallocate(forcingg,forcingspec)
!
!   return
! end subroutine getphytend_old
!

end module phy_run

module run_mod
! time step loop for model run.
! Public subroutines:
! run: main time step loop (advances model state, writes out
! data at specified intervals).

use kinds, only: r_kind,r_double
use params, only: ndimspec, nlevs, ntmax, tstart, dt, nlons, nlats, nlevs,&
  fhzer,ntrac,ntrunc,ntout, explicit, adiabatic, ntrac, iau,&
  massfix,gfsio_out, sigio_out, sfcinitfile, ntdfi, shum, svc, sppt,&
  addnoise, addnoise_vfilt,addnoise_kenorm, addnoise_lscale, addnoise_dissfact,&
  addnoise_vrtonly
use shtns, only: grdtospec, spectogrd, lats, areawts, lap, invlap, degree, getuv
use dyn_run, only: getdyntend, dry_mass_fixer
use phy_run, only: getphytend
use phy_data, only: wrtout_sfc, wrtout_flx, init_phydata, pwat
!! JOY removing calls to write data
!use dyn_init, only: wrtout_sig, wrtout_gfsgrb
use spectral_data, only:  lnpsspec, vrtspec, divspec, virtempspec,&
                          tracerspec, disspec, dmp_prof, diff_prof
! JOY remove Incremental Analysis Unit (IAU) forcing                          
!use iau_module, only: getiauforcing, init_iau, init_iauialized
! these arrays used to print diagnostics after each time step.
use pressure_data, only: psg
use grid_data, only: ug,vg,dlnpsdt
use physcons, only: rerth => con_rerth

use iso_c_binding
private

real(r_double) :: t
bind(c) :: t
public :: run,takeOneStep, t

contains

subroutine takeOneStep() bind(c,name='gfsTakeOneStep')

    integer nt,ntstart
    real(r_kind) fh,pstendmean,spdmax
    real(8) t1,t2
    real(r_kind), dimension(nlons,nlats) :: pstend
    integer(8) count, count_rate, count_max

    call system_clock(count, count_rate, count_max)
    t1 = count*1.d0/count_rate
    ! advance solution with RK
    call advance(t)
    t = t + dt ! absolute forecast time.
    fh = t/3600.
    call system_clock(count, count_rate, count_max)
    t2 = count*1.d0/count_rate
    spdmax = maxval(sqrt(ug**2+vg**2)) ! max wind speed
    pstend = (36.*psg*dlnpsdt)**2 ! ps tend variance (mb/hr)**2
    pstendmean = sqrt(sum(pstend*areawts))
    if (modulo(fh, 10.) .eq. 0.) then
        write(6,8998) fh,spdmax,minval(psg/100.),maxval(psg/100.),pstendmean,t2-t1
        8998 format('t = ',f0.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', pstend = ',f0.3,', cpu time = ',f0.3)
    endif


end subroutine takeOneStep

subroutine advance(t)
    ! advance model state to next time step using 3-stage 2nd-order additive
    ! (i.e. implicit-explicit (IMEX) or semi-implicit) runge-kutta scheme.
!
! Double Butcher tableau looks like this:
!  *explicit part*                *implicit part* 
!
!   0    | 0                       0  | 0
!  a21   | a21 0             aa21+aa22| aa21 aa22
! a31+a32| a31 a32 0    aa31+aa32+aa33| aa31 aa32 aa33
!   1    | b1  b2  b3 0            1  | bb1  bb2   bb3 bb4
! ---------------------  ----------------------------------
!          b1  b2  b3 0                 bb1  bb2   bb3 bb4
!
! the values of the coefficients are defined in semimp_data.f90
!
! progression in time: to time t0+a21*dt, then t0+(a31+a32)*dt, then dt.
! for progression in time to be the same in implicit and
! explicit parts, a21=aa21+a22 and a31+a32=aa31+aa32+aa33
! consistency condition:
! b1 + b2 + b3 = bb1 + bb2 + bb3 + bb4 = 1
! accuracy condition (2nd order):
! a21*b2+(a31+a32)*b3 = 1/2
! (aa21+aa22)*bb2+(aa31+aa32+aa33)*bb3+bb4 = 1/2

! JOY remove patterngenerator and stoch references    
!  use patterngenerator, only: patterngenerator_advance
!  use stoch_data, only: rpattern_svc,rpattern_sppt,&
!     spec_svc,spec_sppt,grd_svc,grd_sppt,&
!     spec_shum,grd_shum,rpattern_shum,getstochforcing,&
!     rpattern_addnoise,specpsi_addnoise,vfact_addnoise
  use semimp_data, only: amhyb,bmhyb,svhyb,d_hyb_m,tor_hyb,&
                         a21,a31,a32,aa21,aa22,aa31,aa32,aa33,b1,b2,b3,bb1,bb2,bb3,bb4
  real(r_double), intent(in) :: t
  complex(r_kind),dimension(ndimspec) :: &
      lnpsspec_orig, dlnpsdtlin1, dlnpsdtlin2,&
      dlnpsspecdt_orig,dlnpsdtlin_orig,&
      dlnpsspecdt1, dlnpsspecdt2
  complex(r_kind), dimension(ndimspec,nlevs) :: &
      vrtspec_orig,divspec_orig,virtempspec_orig,&
      dvrtspecdt_orig,ddivspecdt_orig,dvirtempspecdt_orig,&
      dvrtspecdt1,ddivspecdt1,dvirtempspecdt1,&
      dvrtspecdt2,ddivspecdt2,dvirtempspecdt2,&
      psiforcing,tvforcing,ddivdtlin_orig, dtvdtlin_orig,&
      ddivdtlin1, ddivdtlin2, dtvdtlin1, dtvdtlin2
  complex(r_kind), dimension(ndimspec,nlevs,ntrac) :: &
      dtracerspecdt_orig,&
      tracerspec_orig,dtracerspecdt1,dtracerspecdt2
  complex(r_kind) :: rhs(nlevs)
  complex(r_kind), allocatable, dimension(:,:) :: dvrtspecdt_iau,ddivspecdt_iau,dvirtempspecdt_iau
  complex(r_kind), allocatable, dimension(:,:,:) :: dtracerspecdt_iau
  complex(r_kind), allocatable, dimension(:)  :: dlnpsspecdt_iau
  integer nt, k, n
  logical :: profile = .false. ! print out timing stats
  integer(8) count, count_rate, count_max
  real(8) t1,t2,t0,t4

  ! if iau, allocate space for iau tendencies.
! JOY remove Incremental Analysis Unit (IAU) forcing                          
!  if (iau) then
!     if (.not. init_iauialized) call init_iau()
!     allocate(dvrtspecdt_iau(ndimspec,nlevs))
!     allocate(ddivspecdt_iau(ndimspec,nlevs))
!     allocate(dvirtempspecdt_iau(ndimspec,nlevs))
!     allocate(dtracerspecdt_iau(ndimspec,nlevs,ntrac))
!     allocate(dlnpsspecdt_iau(ndimspec))
!     call getiauforcing(dvrtspecdt_iau,ddivspecdt_iau,dvirtempspecdt_iau,dtracerspecdt_iau,dlnpsspecdt_iau,t)
!  endif

  ! save original fields.
!$omp workshare
  vrtspec_orig = vrtspec
  divspec_orig = divspec
  virtempspec_orig = virtempspec
  tracerspec_orig = tracerspec
!$omp end workshare
  lnpsspec_orig = lnpsspec
  call system_clock(count, count_rate, count_max)
  t0 = count*1.d0/count_rate

  ! compute random patters for stochastic physics to be used this timestep.
  ! JOY removing code which calls patterngenerator_advance
!  if (svc > tiny(svc)) then
!     call patterngenerator_advance(spec_svc,rpattern_svc)
!     call spectogrd(spec_svc,grd_svc)
!  endif
!  if (sppt > tiny(sppt)) then
!     call patterngenerator_advance(spec_sppt,rpattern_sppt)
!     call spectogrd(spec_sppt,grd_sppt)
     ! logit transform to prevent changing the sign of tendencies.
!     grd_sppt = (2./(1.+exp(grd_sppt)))-1.
!  endif
!  if (shum > tiny(shum)) then
!     call patterngenerator_advance(spec_shum,rpattern_shum)
!     call spectogrd(spec_shum,grd_shum)
!  endif
  ! additive noise perturbation.
!  if (addnoise > tiny(addnoise)) then
!     do k=1,nlevs
!        call patterngenerator_advance(specpsi_addnoise(:,k),rpattern_addnoise)
!     enddo
!  endif

  ! stage 1
  ! compute dynamics tendencies.
  call system_clock(count, count_rate, count_max)
  t1 = count*1.d0/count_rate
  call getdyntend(dvrtspecdt_orig,ddivspecdt_orig,dvirtempspecdt_orig,dtracerspecdt_orig,dlnpsspecdt_orig,1)
  call system_clock(count, count_rate, count_max)
  t2 = count*1.d0/count_rate
  ! add IAU contribution (constant over RK sub-steps).
  if (iau) then
!$omp workshare
     dvrtspecdt_orig = dvrtspecdt_orig + dvrtspecdt_iau
     ddivspecdt_orig = ddivspecdt_orig + ddivspecdt_iau
     dvirtempspecdt_orig = dvirtempspecdt_orig + dvirtempspecdt_iau
     dtracerspecdt_orig = dtracerspecdt_orig + dtracerspecdt_iau
     dlnpsspecdt_orig = dlnpsspecdt_orig + dlnpsspecdt_iau
!$omp end workshare
  endif
  !if (profile) print *,'time in getdyntend (stage 1) =',t2-t1
  !print *, dt, a21, aa21
  ! update vorticity and tracers (always explicit)
!$omp workshare
  vrtspec=vrtspec_orig+a21*dt*dvrtspecdt_orig
  tracerspec=tracerspec_orig+a21*dt*dtracerspecdt_orig
!$omp end workshare
  if (explicit) then
!$omp workshare
     divspec=divspec_orig+a21*dt*ddivspecdt_orig
     virtempspec=virtempspec_orig+a21*dt*dvirtempspecdt_orig
     lnpsspec=lnpsspec_orig+a21*dt*dlnpsspecdt_orig
!$omp end workshare
  else
! solve for updated divergence.
! back subsitution to get updated virt temp, lnps.
!$omp parallel do private(n,rhs)
  do n=1,ndimspec
! first remove linear terms from computed tendencies.
      ddivdtlin_orig(n,:) = -lap(n)*& 
      (matmul(amhyb,virtempspec_orig(n,:)) + tor_hyb(:)*lnpsspec_orig(n))
      ddivspecdt_orig(n,:) = ddivspecdt_orig(n,:) - ddivdtlin_orig(n,:)
      dtvdtlin_orig(n,:) = -matmul(bmhyb,divspec_orig(n,:))
      dvirtempspecdt_orig(n,:) = dvirtempspecdt_orig(n,:) - dtvdtlin_orig(n,:)
      dlnpsdtlin_orig(n) = -sum(svhyb(:)*divspec_orig(n,:))
      dlnpsspecdt_orig(n) = dlnpsspecdt_orig(n) - dlnpsdtlin_orig(n)
! compute explicit part of update
      divspec(n,:) = divspec_orig(n,:) + dt*(a21*ddivspecdt_orig(n,:)  +&
                                             aa21*ddivdtlin_orig(n,:))
      virtempspec(n,:) = virtempspec_orig(n,:) + dt*(a21*dvirtempspecdt_orig(n,:)  +&
                                                     aa21*dtvdtlin_orig(n,:))
      lnpsspec(n) = lnpsspec_orig(n) + dt*(a21*dlnpsspecdt_orig(n) +&
                                           aa21*dlnpsdtlin_orig(n))
! solve implicit part for updated divergence.
      rhs(:)   = divspec(n,:) - aa22*lap(n)*dt*&
                 (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
      divspec(n,:) = matmul(d_hyb_m(:,:,degree(n)+1,1),rhs(:))
! back substitute to get implicit contribution to virtual temp and lnps update.
      virtempspec(n,:) = virtempspec(n,:) - aa22*dt*matmul(bmhyb,divspec(n,:))
      lnpsspec(n) = lnpsspec(n) - aa22*dt*sum(svhyb(:)*divspec(n,:))
  enddo
!$omp end parallel do 
  endif

  ! stage 2
  ! compute dynamics tendencies.
  call system_clock(count, count_rate, count_max)
  t1 = count*1.d0/count_rate
  call getdyntend(dvrtspecdt1,ddivspecdt1,dvirtempspecdt1,dtracerspecdt1,dlnpsspecdt1,2)
  call system_clock(count, count_rate, count_max)
  t2 = count*1.d0/count_rate
  !if (profile) print *,'time in getdyntend (stage 2) =',t2-t1
  ! add IAU contribution (constant over RK sub-steps).
  if (iau) then
!$omp workshare
     dvrtspecdt1 = dvrtspecdt1 + dvrtspecdt_iau
     ddivspecdt1 = ddivspecdt1 + ddivspecdt_iau
     dvirtempspecdt1 = dvirtempspecdt1 + dvirtempspecdt_iau
     dtracerspecdt1 = dtracerspecdt1 + dtracerspecdt_iau
     dlnpsspecdt1 = dlnpsspecdt1 + dlnpsspecdt_iau
!$omp end workshare
  endif
  ! update vorticity and tracers (always explicit)
!$omp workshare
  vrtspec=vrtspec_orig+dt*(a31*dvrtspecdt_orig+a32*dvrtspecdt1)
  tracerspec=tracerspec_orig+dt*(a31*dtracerspecdt_orig+a32*dtracerspecdt1)
!$omp end workshare
  if (explicit) then
!$omp workshare
     divspec=divspec_orig+dt*(a31*ddivspecdt_orig+a32*ddivspecdt1)
     virtempspec=virtempspec_orig+dt*(a31*dvirtempspecdt_orig+a32*dvirtempspecdt1)
     lnpsspec=lnpsspec_orig+dt*(a31*dlnpsspecdt_orig+a32*dlnpsspecdt1)
!$omp end workshare
  else
! solve for updated divergence.
! back subsitution to get updated virt temp, lnps.
!$omp parallel do private(n,rhs)
  do n=1,ndimspec
! remove linear terms from computed tendencies.
      ddivdtlin1(n,:) = -lap(n)*& 
      (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
      ddivspecdt1(n,:) = ddivspecdt1(n,:) - ddivdtlin1(n,:)
      dtvdtlin1(n,:) = -matmul(bmhyb,divspec(n,:))
      dvirtempspecdt1(n,:) = dvirtempspecdt1(n,:) - dtvdtlin1(n,:)
      dlnpsdtlin1(n) = -sum(svhyb(:)*divspec(n,:))
      dlnpsspecdt1(n) = dlnpsspecdt1(n) - dlnpsdtlin1(n)
! compute explicit part of update
      divspec(n,:) = divspec_orig(n,:) + dt*(a31*ddivspecdt_orig(n,:)  +&
                                             aa31*ddivdtlin_orig(n,:)  +&
                                             a32*ddivspecdt1(n,:)      +&
                                             aa32*ddivdtlin1(n,:))   
      virtempspec(n,:) = virtempspec_orig(n,:) + dt*(a31*dvirtempspecdt_orig(n,:)  +&
                                                     aa31*dtvdtlin_orig(n,:)       +&
                                                     a32*dvirtempspecdt1(n,:)      +&
                                                     aa32*dtvdtlin1(n,:))
      lnpsspec(n) = lnpsspec_orig(n) + dt*(a31*dlnpsspecdt_orig(n)     +&
                                           aa31*dlnpsdtlin_orig(n)     +&
                                           a32*dlnpsspecdt1(n)         +&
                                           aa32*dlnpsdtlin1(n))     
! solve implicit part for updated divergence.
      rhs(:)   = divspec(n,:) - aa33*lap(n)*dt*&
                 (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
      divspec(n,:) = matmul(d_hyb_m(:,:,degree(n)+1,2),rhs(:))
! back substitute to get implicit contribution to virtual temp and lnps update.
      virtempspec(n,:) = virtempspec(n,:) - aa33*dt*matmul(bmhyb,divspec(n,:))
      lnpsspec(n) = lnpsspec(n) - aa33*dt*sum(svhyb(:)*divspec(n,:))
  enddo
!$omp end parallel do 
  endif

  ! stage 3
  ! compute dynamics tendencies.
  call system_clock(count, count_rate, count_max)
  t1 = count*1.d0/count_rate
  call getdyntend(dvrtspecdt2,ddivspecdt2,dvirtempspecdt2,dtracerspecdt2,dlnpsspecdt2,3)
  call system_clock(count, count_rate, count_max)
  t2 = count*1.d0/count_rate
  !if (profile) print *,'time in getdyntend (stage 3) =',t2-t1
  ! add IAU contribution (constant over RK sub-steps).
  if (iau) then
!$omp workshare
     dvrtspecdt2 = dvrtspecdt2 + dvrtspecdt_iau
     ddivspecdt2 = ddivspecdt2 + ddivspecdt_iau
     dvirtempspecdt2 = dvirtempspecdt2 + dvirtempspecdt_iau
     dtracerspecdt2 = dtracerspecdt2 + dtracerspecdt_iau
     dlnpsspecdt2 = dlnpsspecdt2 + dlnpsspecdt_iau
!$omp end workshare
  endif
  ! final update of vorticity and tracers (always explicit)
!$omp workshare
  vrtspec=vrtspec_orig+dt*(b1*dvrtspecdt_orig+b2*dvrtspecdt1+b3*dvrtspecdt2)
  tracerspec=tracerspec_orig+dt*(b1*dtracerspecdt_orig+b2*dtracerspecdt1+b3*dtracerspecdt2)
!$omp end workshare
  if (explicit) then
!$omp workshare
     divspec=divspec_orig+dt*(b1*ddivspecdt_orig+b2*ddivspecdt1+b3*ddivspecdt2)
     virtempspec=virtempspec_orig+dt*(b1*dvirtempspecdt_orig+b2*dvirtempspecdt1+b3*dvirtempspecdt2)
     lnpsspec=lnpsspec_orig+dt*(b1*dlnpsspecdt_orig+b2*dlnpsspecdt1+b3*dlnpsspecdt2)
!$omp end workshare
  else
!$omp parallel do private(n,rhs)
  do n=1,ndimspec
! remove linear terms from computed tendencies.
      ddivdtlin2(n,:) = -lap(n)*& 
      (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
      ddivspecdt2(n,:) = ddivspecdt2(n,:) - ddivdtlin2(n,:)
      dtvdtlin2(n,:) = -matmul(bmhyb,divspec(n,:))
      dvirtempspecdt2(n,:) = dvirtempspecdt2(n,:) - dtvdtlin2(n,:)
      dlnpsdtlin2(n) = -sum(svhyb(:)*divspec(n,:))
      dlnpsspecdt2(n) = dlnpsspecdt2(n) - dlnpsdtlin2(n)
! compute final update
      divspec(n,:) = divspec_orig(n,:) + dt*(b1*ddivspecdt_orig(n,:)  +&
                                             bb1*ddivdtlin_orig(n,:)  +& 
                                             b2*ddivspecdt1(n,:)      +&
                                             bb2*ddivdtlin1(n,:)      +& 
                                             b3*ddivspecdt2(n,:)      +&
                                             bb3*ddivdtlin2(n,:))       
      virtempspec(n,:) = virtempspec_orig(n,:) + dt*(b1*dvirtempspecdt_orig(n,:)  +&
                                                     bb1*dtvdtlin_orig(n,:)       +& 
                                                     b2*dvirtempspecdt1(n,:)      +&
                                                     bb2*dtvdtlin1(n,:)           +& 
                                                     b3*dvirtempspecdt2(n,:)      +&
                                                     bb3*dtvdtlin2(n,:))       
      lnpsspec(n) = lnpsspec_orig(n) + dt*(b1*dlnpsspecdt_orig(n) +&
                                           bb1*dlnpsdtlin_orig(n) +&
                                           b2*dlnpsspecdt1(n)     +&
                                           bb2*dlnpsdtlin1(n)     +&
                                           b3*dlnpsspecdt2(n)     +&
                                           bb3*dlnpsdtlin2(n))
! solve implicit part for updated divergence.
      if (abs(bb4) .gt. 1.e-5) then
         rhs(:)   = divspec(n,:) - bb4*lap(n)*dt*&
                    (matmul(amhyb,virtempspec(n,:)) + tor_hyb(:)*lnpsspec(n))
         divspec(n,:) = matmul(d_hyb_m(:,:,degree(n)+1,3),rhs(:))
! back substitute to get implicit contribution to virtual temp and lnps update.
         virtempspec(n,:) = virtempspec(n,:) - bb4*dt*matmul(bmhyb,divspec(n,:))
         lnpsspec(n) = lnpsspec(n) - bb4*dt*sum(svhyb(:)*divspec(n,:))
      endif
  enddo
!$omp end parallel do 
  endif

  ! forward implicit treatment of linear damping/diffusion
  !$omp parallel do private(k,nt)
  do k=1,nlevs
     vrtspec(:,k) = vrtspec(:,k)/(1. - (disspec(:)*diff_prof(k) - dmp_prof(k))*dt)
     divspec(:,k) = divspec(:,k)/(1. - (disspec(:)*diff_prof(k) - dmp_prof(k))*dt)
     virtempspec(:,k) = virtempspec(:,k)/(1. - disspec(:)*diff_prof(k)*dt)
     do nt=1,ntrac
        tracerspec(:,k,nt) = tracerspec(:,k,nt)/(1. - disspec(:)*diff_prof(k)*dt)
     enddo
  enddo
  !$omp end parallel do

  call system_clock(count, count_rate, count_max)
  t4 = count*1.d0/count_rate
  !if (profile) print *,'time in dynamics update=',t4-t0
  ! apply physics parameterizations as an adjustment to fields updated by dynamics.
  ! (time-split approach)
  if (.not. adiabatic) then
     ! update variables on grid.
     call getdyntend(dvrtspecdt1,ddivspecdt1,dvirtempspecdt1,&
          dtracerspecdt1,dlnpsspecdt1,3,.true.) ! <- .true. for spectogrd only
     ! compute stochastic backscatter forcing, if needed
     ! if addnoise_dissfact is true,
     ! streamfunctione noise forcing by smoothed estimate of ke dissipation
     ! rate.
     ! Joy removing stochastic forcing routines
     !if (addnoise > tiny(addnoise)) &
     !call getstochforcing(specpsi_addnoise, psiforcing, tvforcing)
     ! compute physics tendencies, apply as an adjustment to updated state
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     call getphytend(dvrtspecdt1,ddivspecdt1,dvirtempspecdt1,dtracerspecdt1,dlnpsspecdt1,t,dt)
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
     ! if (profile) print *,'time in getphytend=',t2-t1
!$omp workshare
     vrtspec=vrtspec+dt*dvrtspecdt1
     divspec=divspec+dt*ddivspecdt1
     virtempspec=virtempspec+dt*dvirtempspecdt1
     tracerspec=tracerspec+dt*dtracerspecdt1
!$omp end workshare
! modify spectral lnps tendency to include contribution
! from dry mass 'fixer' (adjusts dry surface pressure to
! remain equal to pdryini).
     if (massfix) then
        call dry_mass_fixer(psg,pwat,dlnpsspecdt1,dt)
     endif
     lnpsspec=lnpsspec+dt*dlnpsspecdt1 ! only needed for dry mass fixer.
     ! add stochastic backsscatter forcing in vorticity eqn.
     ! psiforcing has units of streamfunction tendency
     if (addnoise > tiny(addnoise)) then
         !$omp parallel do private(k)
         do k=1,nlevs
            if (addnoise_kenorm .and. addnoise_vrtonly) then
               !vrtspec(:,k) = vrtspec(:,k) + dt*(lap/rerth**2)*sqrt(-invlap*rerth**2)*psiforcing(:,k)
               vrtspec(:,k) = vrtspec(:,k) - dt*(sqrt(-lap)/rerth)*psiforcing(:,k)
            else
               vrtspec(:,k) = vrtspec(:,k) + dt*(lap/rerth**2)*psiforcing(:,k)
            endif
            virtempspec(:,k) = virtempspec(:,k) + dt*tvforcing(:,k)
         enddo
         !$omp end parallel do 
     end if
  end if

  ! free-up temporary storage.
  if (iau) then
     deallocate(dvrtspecdt_iau)
     deallocate(ddivspecdt_iau)
     deallocate(dvirtempspecdt_iau)
     deallocate(dtracerspecdt_iau)
     deallocate(dlnpsspecdt_iau)
  endif

end subroutine advance

subroutine set_dfi_wts(dfi_wts)
 ! set Lanczos filter weights for digital filter 
 real(r_kind), intent(out) :: dfi_wts(0:2*ntdfi)
 real(r_kind) totsum,sx,wx
 integer kstep
 dfi_wts = 1.
 totsum = 0.
 do kstep=0,2*ntdfi
    sx     = acos(-1.)*(kstep-ntdfi)/ntdfi
    wx     = acos(-1.)*(kstep-ntdfi)/(ntdfi+1)
    if (kstep .NE. ntdfi) then
       dfi_wts(kstep) = sin(wx)/wx*sin(sx)/sx
    endif
    totsum = totsum + dfi_wts(kstep)
 enddo
 dfi_wts = dfi_wts/totsum
 print *,'lanczos dfi wts:'
 do kstep=0,2*ntdfi
    print *,kstep-ntdfi,dfi_wts(kstep)
 enddo
end subroutine set_dfi_wts

subroutine run()
  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
  implicit none
  integer nt,ntstart,my_id
  real(r_kind) fh,fha,pstendmean,spdmax
  real(r_double) ta,t
  real(8) t1,t2
  real(r_kind), dimension(nlons,nlats) :: pstend
  complex(r_kind), dimension(:,:), allocatable :: &
  vrtspec_dfi,divspec_dfi,virtempspec_dfi
  complex(r_kind), dimension(:,:,:), allocatable :: tracerspec_dfi
  complex(r_kind), dimension(:), allocatable :: lnpsspec_dfi
  real(r_kind), dimension(:), allocatable :: dfi_wts
  integer(8) count, count_rate, count_max
  character(len=500) filename,filename_save

!$omp parallel
  my_id = omp_get_thread_num()
  if (my_id .eq. 0) print *,'running with',omp_get_num_threads(),' threads'
!$omp end parallel

  t = tstart
  ta = 0.
  ntstart = 1

  !JOY not doing the DFI loop
  ! digital filter loop.
!  if (ntdfi > 0) then
!     print *,'in dfi time step loop...'
     ! allocate work space
!     allocate(vrtspec_dfi(ndimspec,nlevs),divspec_dfi(ndimspec,nlevs))
 !    allocate(virtempspec_dfi(ndimspec,nlevs),lnpsspec_dfi(ndimspec))
 !    allocate(tracerspec_dfi(ndimspec,nlevs,ntrac),dfi_wts(0:2*ntdfi))
     ! compute dfi weights
!     call set_dfi_wts(dfi_wts)
     ! intialize weighted time averages.
!!$omp workshare
!     vrtspec_dfi = dfi_wts(0)*vrtspec
!     divspec_dfi = dfi_wts(0)*divspec
!     virtempspec_dfi = dfi_wts(0)*virtempspec
!     tracerspec_dfi = dfi_wts(0)*tracerspec
!!$omp end workshare
!     lnpsspec_dfi = dfi_wts(0)*lnpsspec
!     do nt=1,2*ntdfi
!        call system_clock(count, count_rate, count_max)
!        t1 = count*1.d0/count_rate
!        call advance(t)
!        t = t + dt ! absolute forecast time.
!        ta = ta + dt ! absolute forecast time.
!        fh = t/3600.
!!$omp workshare
!        vrtspec_dfi = vrtspec_dfi + dfi_wts(nt)*vrtspec
!        divspec_dfi = divspec_dfi + dfi_wts(nt)*divspec
!        virtempspec_dfi = virtempspec_dfi + dfi_wts(nt)*virtempspec
!        tracerspec_dfi = tracerspec_dfi + dfi_wts(nt)*tracerspec
!!$omp end workshare
!        lnpsspec_dfi = lnpsspec_dfi + dfi_wts(nt)*lnpsspec
!        call system_clock(count, count_rate, count_max)
!        t2 = count*1.d0/count_rate
!        spdmax = maxval(sqrt(ug**2+vg**2)) ! max wind speed
!        pstend = (36.*psg*dlnpsdt)**2 ! ps tend variance (mb/hr)**2
!        pstendmean = sqrt(sum(pstend*areawts))
!        write(6,8998) fh,spdmax,minval(psg/100.),maxval(psg/100.),pstendmean,t2-t1
        ! write out surface and flux data in middle of dfi window.
!        if (nt .eq. ntdfi) then
!           write(filename_save,9000) nint(fh)
!           print *,'writing to ',trim(filename_save),' fh=',fh
!           call wrtout_sfc(fh,filename_save)
!           write(filename,9001) nint(fh)
!           print *,'writing to ',trim(filename),' fh=',fh
!           call wrtout_flx(fh,ta,filename)
!        ! write first time step output
!        else if (nt .eq. 1) then
!           if (sigio_out) then
!              write(filename,8999) nint(fh)
!              print *,'writing to ',trim(filename),' fh=',fh
!              call wrtout_sig(fh,filename)
!           endif
!           write(filename,9000) nint(fh)
!           print *,'writing to ',trim(filename),' fh=',fh
!           call wrtout_sfc(fh,filename)
!           write(filename,9001) nint(fh)
!           print *,'writing to ',trim(filename),' fh=',fh
!           call wrtout_flx(fh,ta,filename)
!           if (gfsio_out) then
!              write(filename,9002) nint(fh)
!              print *,'writing to ',trim(filename),' fh=',fh
!              call wrtout_gfsgrb(fh,filename)
!           endif
!        end if
!     enddo
!     print *,'done with dfi loop, resetting fields and restarting...'
!     ! reset model state to weighted time average.
!     vrtspec = vrtspec_dfi; divspec = divspec_dfi; virtempspec = virtempspec_dfi
!     lnpsspec = lnpsspec_dfi; tracerspec = tracerspec_dfi
!     ! deallocate work space.
!     deallocate(vrtspec_dfi,divspec_dfi,virtempspec_dfi,lnpsspec_dfi,tracerspec_dfi,dfi_wts)
!     ! reset surface data to values at middle of window (also zeros flux arrays).
!     sfcinitfile = filename_save; call init_phydata(); ta = 0.
     ! reset time.
!     t = tstart + ntdfi*dt; ntstart = ntdfi+1
     ! write out spectral data after dfi.
!     fh = t/3600.
!     if (sigio_out) then
!        write(filename,8999) nint(fh)
!        print *,'writing to ',trim(filename),' fh=',fh
!        call wrtout_sig(fh,filename)
!     endif
!     ! write out gaussian grid data after DFI
!     if (gfsio_out) then
!        write(filename,9002) nint(fh)
!        print *,'writing to ',trim(filename),' fh=',fh
!        call wrtout_gfsgrb(fh,filename)
!     endif
!  endif

  ! main time step loop
  do nt=ntstart,ntmax
     call system_clock(count, count_rate, count_max)
     t1 = count*1.d0/count_rate
     ! advance solution with RK
     call advance(t)
     t = t + dt ! absolute forecast time.
     ta = ta + dt ! time in accumulaton interval.
     fh = t/3600.
     fha = ta/3600.
     call system_clock(count, count_rate, count_max)
     t2 = count*1.d0/count_rate
     spdmax = maxval(sqrt(ug**2+vg**2)) ! max wind speed
     pstend = (36.*psg*dlnpsdt)**2 ! ps tend variance (mb/hr)**2
     pstendmean = sqrt(sum(pstend*areawts))
     write(6,8998) fh,spdmax,minval(psg/100.),maxval(psg/100.),pstendmean,t2-t1
8998 format('t = ',f0.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', pstend = ',f0.3,', cpu time = ',f0.3)
     ! write out data at specified intervals.
     ! data always written at first time step.
!     if (nt .eq. 1 .or. (ntout .ne. 0 .and. mod(nt,ntout) .eq. 0)) then
!        if (sigio_out) then
!           write(filename,8999) nint(fh)
!8999       format('SIG.F',i0.2) ! at least three digits used
!           print *,'writing to ',trim(filename),' fh=',fh
!           call wrtout_sig(fh,filename)
!        endif    
!        ! write out boundary and flux files if using gfs physics.
!        write(filename,9000) nint(fh)
!9000    format('SFC.F',i0.2) ! at least three digits used
!        print *,'writing to ',trim(filename),' fh=',fh
!        call wrtout_sfc(fh,filename)
!        write(filename,9001) nint(fh)
!9001    format('FLX.F',i0.2) ! at least three digits used
!        print *,'writing to ',trim(filename),' fh=',fh
!        call wrtout_flx(fh,ta,filename)
!        if (gfsio_out) then
!           write(filename,9002) nint(fh)
!9002       format('GFS.F',i0.2) ! at least three digits used
!           print *,'writing to ',trim(filename),' fh=',fh
!           call wrtout_gfsgrb(fh,filename)
!        endif
!     end if
     if (abs(fha-fhzer) .lt. 1.e-5) ta=0. ! reset accum time.
  enddo

end subroutine run

end module run_mod

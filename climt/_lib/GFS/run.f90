module run_mod
! time step loop for model run.
! Public subroutines:
! run: main time step loop (advances model state, writes out
! data at specified intervals).

use kinds, only: r_kind,r_double
use params, only: ndimspec, nlevs, ntmax, tstart, dt, nlons, nlats, nlevs,&
  fhzer,ntrac,ntrunc,ntout, explicit, adiabatic, ntrac, &
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
implicit none
private

real(c_double) :: t
bind(c) :: t
public :: take_one_step, t

contains

subroutine set_model_time(py_time) bind(c, name='gfs_set_model_time')
    real(c_double), intent(in):: py_time
    t = py_time
    !print *, 'model time is now: ', t
end subroutine


subroutine take_one_step() bind(c,name='gfs_take_one_step')

    !integer nt,ntstart
    !real(r_kind) fh,pstendmean,spdmax
    real(r_kind) fh
    real(8) t1,t2
    !real(r_kind), dimension(nlons,nlats) :: pstend
    integer(8) count, count_rate, count_max

    call system_clock(count, count_rate, count_max)
    t1 = count*1.d0/count_rate
    ! advance solution with RK
    call advance(t)
    fh = t/3600.
    t = t + dt ! absolute forecast time.
    call system_clock(count, count_rate, count_max)
    t2 = count*1.d0/count_rate
    !spdmax = maxval(sqrt(ug**2+vg**2)) ! max wind speed
    !pstend = (36.*psg*dlnpsdt)**2 ! ps tend variance (mb/hr)**2
    !pstendmean = sqrt(sum(pstend*areawts))
    !if (modulo(fh, 10.) .eq. 0.) then
    !    write(6,8998) fh,spdmax,minval(psg/100.),maxval(psg/100.),pstendmean,t2-t1
    !    8998 format('t = ',f0.3,' hrs, spdmax = ',f7.3,', min/max ps = ',f7.2,'/',f7.2,', pstend = ',f0.3,', cpu time = ',f0.3)
    !endif

end subroutine take_one_step

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
  integer nt, k, n
  !logical :: profile = .false. ! print out timing stats
  integer(8) count, count_rate, count_max
  real(8) t1,t2,t0,t4


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


  ! stage 1
  ! compute dynamics tendencies.
  call system_clock(count, count_rate, count_max)
  t1 = count*1.d0/count_rate
  call getdyntend(dvrtspecdt_orig,ddivspecdt_orig,dvirtempspecdt_orig,dtracerspecdt_orig,dlnpsspecdt_orig,1)
  call system_clock(count, count_rate, count_max)
  t2 = count*1.d0/count_rate

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

end subroutine advance

end module run_mod

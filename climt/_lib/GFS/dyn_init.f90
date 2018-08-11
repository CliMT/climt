module dyn_init
! initialize dynamics (including reading in initial conditions).
! also does IO of spectral data.
! public subroutines:
! init_dyn:  allocate spectral arrays (init_specdata),
!  read in initial conditions, initalize pressure variables (init_pressdata,
!  get ak,bk from IC file header), compute gradient of orography,
!  set up linear damping operators (disspec,diff_prof,damp_prof), initialize arrays
!  for semi-implicit time stepping.
! wrtout_sig: write out spectral data.
! readin_sig: read in spectral data.
 use kinds, only: r_kind, r_single, r_double
! JOY removing sigio dependency
! use sigio_module, only: sigio_sclose,sigio_swohdc,&
!  sigio_srohdc,sigio_aldata,sigio_data,sigio_head,sigio_sropen,sigio_srdata,sigio_axdata
 use params, only: &
 nlons,nlats,nlevs,ndimspec,ntrunc,initfile,dry,ndiss,efold,dcmip,polar_opt,&
 heldsuarez,explicit,tstart,idate_start,hdif_fac,hdif_fac2,fshk,ntrac,taustratdamp,&
 ntoz,ntclw,pdryini,massfix,toa_pressure
 !sighead
 !JOY removed sighead to remove sigio dependency
 use shtns, only: shtns_init, spectogrd, grdtospec, getgrad, getvrtdivspec, lap, lats, lons
 use spectral_data, only: vrtspec,divspec,virtempspec,tracerspec,topospec,lnpsspec,&
                          disspec,diff_prof,dmp_prof,init_specdata
 use pressure_data, only: ak,bk,ck,dbk,bkl,sl,si,psg,prs,pk,dpk,calc_pressdata
 use grid_data, only: init_griddata, dphisdx, dphisdy, phis, ug, vg, virtempg, &
 tracerg,lnpsg,dlnpdtg
! JOY removing stochastic routines
! use stoch_data, only: init_stochdata
 use physcons, only: rerth => con_rerth, rd => con_rd, cp => con_cp, &
                     omega => con_omega, grav => con_g, pi => con_pi, &
                     fv => con_fvirt
 use semimp_data, only: init_semimpdata
 use iso_c_binding


 implicit none
!namelist/dynamics_list/ak,bk
 private
 !public :: init_dyn, wrtout_sig, readin_sig, wrtout_gfsgrb
 !! removing all output calls
 ! JOY remove routines which read initial conditions from file
 !public :: init_dyn, wrtout_sig, readin_sig
 public :: init_dyn

 real(r_kind) :: slrd0=0
 contains
 !JOY subroutine to set topography from arbitrary array
 subroutine set_topography(surf_geop) bind(c, name='gfs_set_topography')

    real(r_kind), intent(in), dimension(nlons,nlats) :: surf_geop

    call grdtospec(surf_geop/grav, topospec)
    call spectogrd(grav*topospec, phis)
    call getgrad(grav*topospec, dphisdx, dphisdy, rerth)

 end subroutine set_topography

 subroutine init_dyn(py_damping_level, py_taustratdamp) bind(c, name='gfs_init_dynamics')
    integer k
    ! allocate arrays
    ! JOY data array inits will be done from python

    !call init_specdata()
    !call init_griddata()
    integer(c_int), intent(in):: py_damping_level
    real(c_double), intent(in):: py_taustratdamp

    ! initialize spherical harmonic lib
    !print *,'Initialising shtns'
    call shtns_init(nlons,nlats,ntrunc,nthreads=1,polar_opt=polar_opt)
        ! convert to ln(ps) in Pa.
    !print *,'calling spectogrd'
    !print *, 'Grid size', nlats, nlons, ntrunc
    !print *, shape(lnpsspec)
    !call spectogrd(lnpsspec, psg)
    !psg = 1000.*exp(psg) ! convert to Pa
    !call grdtospec(log(psg), lnpsspec) ! back to spectral.


    do k=1,nlevs
       dbk(k) = bk(k+1)-bk(k)
       bkl(k) = 0.5*(bk(k+1)+bk(k))
       ck(k)  = ak(k+1)*bk(k)-ak(k)*bk(k+1)
    enddo
         ! compute sigma coordinate quantities (bottom to top).
    do k=1,nlevs+1
        si(nlevs+2-k)= (ak(k) - toa_pressure)/(pdryini - toa_pressure)+bk(k) ! si are now sigmas
    enddo
    do k=1,nlevs
        sl(k) = 0.5*(si(k)+si(k+1))
         !sl(k) = ((si(k)**(rk+1.) - si(k+1)**(rk+1.))/&
         !        ((rk+1.)*(si(k)-si(k+1))))**(1./rk)
    enddo

    slrd0 = si(nlevs+1 - py_damping_level)
    taustratdamp = py_taustratdamp

    !call spectogrd(grav*topospec, phis)
    !print *,'min/max surface geopotential',minval(phis),maxval(phis)
    !call spectogrd(lnpsspec, psg)
    !psg = exp(psg)
    !print *,'min/max sfc pressure (hPa)',minval(psg/100.),maxval(psg/100.)
    ! initialize model interface and level pressures, related variables.
    !call calc_pressdata(lnpsg)
    ! compute gradient of surface orography
    !call getgrad(grav*topospec, dphisdx, dphisdy, rerth)
    ! hyper-diffusion operator (plus rayleigh damping in upper stratosphere)

    call setdampspec(ndiss,efold,hdif_fac,hdif_fac2,fshk,disspec,diff_prof,dmp_prof)

    ! initialize arrays for semi-implicit adjustments.
    if (.not. explicit) call init_semimpdata()
    ! initialize stochastic data.
    ! JOY removing stochastic routines
    !call init_stochdata()
 end subroutine init_dyn

 subroutine copyspecin(rspecdata,cspecdata)
    real(r_single), intent(in) :: rspecdata(2*ndimspec)
    complex(r_kind), intent(out) :: cspecdata(ndimspec)
    integer n,nn
    nn = 1
    ! factor of sqrt(2.*pi) accounts for difference in normalization
    ! between ncep libs and shtns (which uses orthonormalized norm)
    do n=1,ndimspec
       cspecdata(n) = sqrt(2.*pi)*cmplx(rspecdata(nn),rspecdata(nn+1))
       nn = nn + 2
    enddo
 end subroutine copyspecin

 subroutine copyspecout(cspecdata, rspecdata)
    real(r_single), intent(out) :: rspecdata(2*ndimspec)
    complex(r_kind), intent(in) :: cspecdata(ndimspec)
    integer n,nn
    nn = 1
    ! factor of sqrt(2.*pi) accounts for difference in normalization
    ! between ncep libs and shtns (which uses orthonormalized norm)
    do n=1,ndimspec
       rspecdata(nn) = real(cspecdata(n))/sqrt(2.*pi)
       rspecdata(nn+1) = imag(cspecdata(n))/sqrt(2.*pi)
       nn = nn + 2
    enddo
 end subroutine copyspecout


 subroutine twodtooned(data2,data1)
   real(r_kind), intent(in) :: data2(nlons,nlats)
   real(r_kind), intent(out) :: data1(nlons*nlats)
   integer i,j,n
   do n=1,nlons*nlats
      j = 1+(n-1)/nlons
      i = n-(j-1)*nlons
      data1(n) = data2(i,j)
   enddo
 end subroutine twodtooned

! subroutine dcmip_ics(testcase,dry)
!   ! DCMIP 2012 test case initial conditions (test cases 4 and 5)
!   ! (idealized baroclinic instability or trop cyclone test case initial conditions)
!   use dcmip_initial_conditions_test_4, only: test4_baroclinic_wave
!   use dcmip_initial_conditions_test_5, only: test5_tropical_cyclone
!   integer, intent(in) :: testcase
!   logical, intent(in) :: dry
!   real(8) lon,lat,z,X,p,u,v,w,t,phi_s,ps,rho,q,q1,q2
!   integer moist,zcoords,i,j,k,nt
!
!   if (testcase/10 .eq. 5) then
!      if (dry) then
!         print *,'cannot run dcmip trop cyclone test case with dry=.true.'
!         stop
!      endif
!      print *,'replacing initial conds with dcmip tropical cyclone case..'
!   else if (testcase/10 .eq. 4) then
!      print *,'replacing initial conds with dcmip baroclinic wave test case..'
!   else
!      print *,'dcmip must be 40 something or 50 something'
!      stop
!   endif
!   zcoords = 0 ! p, not z, is specified.
!   if (dry) then
!     moist = 0
!   else
!     moist = 1
!   endif
!   if (testcase/10 .eq. 4) psg = 1.e5
!   if (testcase/10 .eq. 5) then
!      psg = 1.015e5
!      lnpsg = log(psg)
!      call grdtospec(lnpsg,lnpsspec)
!      call calc_pressdata(lnpsg)
!      ! run once, just to calculate ps
!      ! ps then used to compute pressure on model levels for subsequent call
!      do k=1,nlevs
!      do j=1,nlats
!      do i=1,nlons
!         lon = lons(i,j); lat = lats(i,j); p = prs(i,j,k)
!         call test5_tropical_cyclone (lon,lat,p,z,zcoords,u,v,w,t,phi_s,ps,rho,q)
!         psg(i,j) = ps
!      enddo
!      enddo
!      enddo
!   endif
!   ! calculate pressure on model levels given ps.
!   lnpsg = log(psg)
!   call grdtospec(lnpsg,lnpsspec)
!   call calc_pressdata(lnpsg)
!   ! calculate u,v,t,q,phis and other tracers.
!   do k=1,nlevs
!      do j=1,nlats
!      do i=1,nlons
!         lon = lons(i,j); lat = lats(i,j); p = prs(i,j,k)
!         if (testcase/10 .eq. 4) then
!            call test4_baroclinic_wave (moist,X,lon,lat,p,z,zcoords,u,v,w,t,phi_s,ps,rho,q,q1,q2)
!         else if (testcase/10 .eq. 5) then
!            call test5_tropical_cyclone (lon,lat,p,z,zcoords,u,v,w,t,phi_s,ps,rho,q)
!         endif
!         ug(i,j,k) = u; vg(i,j,k) = v; phis(i,j) = phi_s
!         if (dry) then
!            virtempg(i,j,k) = t
!            tracerg(i,j,k,1) = 0.
!            tracerg(i,j,k,2) = q1
!            tracerg(i,j,k,3) = q2
!         else
!            virtempg(i,j,k) = t*(1. + fv*q)
!            tracerg(i,j,k,1) = q
!            if (testcase/10 .eq. 4) then
!               tracerg(i,j,k,2) = q1
!               tracerg(i,j,k,3) = q2
!            else if (testcase .eq. 52) then
!               tracerg(i,j,k,2:ntrac) = 0.
!            endif
!         endif
!       enddo
!       enddo
!       call getvrtdivspec(ug(:,:,k),vg(:,:,k),vrtspec(:,k),divspec(:,k),rerth)
!       call grdtospec(virtempg(:,:,k),virtempspec(:,k))
!       do nt=1,ntrac
!          call grdtospec(tracerg(:,:,k,nt),tracerspec(:,k,nt))
!       enddo
!   enddo
!   call grdtospec(phis/grav, topospec)
!
! end subroutine dcmip_ics


 subroutine heldsuarez_ics(toa_pressure)
   ! jablonowski and williamson (2006, QJR, p. 2943, doi: 10.1256/qj.06.12)
   ! initial perturbation on an isothermal state.
   real(r_kind), dimension(nlons,nlats) :: rnh,xnh,rsh,xsh, rand_vel
   real(r_kind) :: lonc,latc,up,pertrad
   real(c_double), intent(in):: toa_pressure
   integer k
   print *,'replacing initial conds with held and suarez test case..'
   lonc = pi/9.
   up = 1.
   pertrad = rerth/10.
   !JOY changing ICs
   latc = 0.
   !latc = 2.*pi/9.
   !xnh = sin(latc)*sin(lats) + cos(latc)*cos(lats)*cos(lons-lonc)
   xnh = sin(latc)*sin(lats) + cos(latc)*cos(lats)
   latc = -2.*pi/9.
   !xsh = sin(latc)*sin(lats) + cos(latc)*cos(lats)*cos(lons-lonc)
   xsh = sin(latc)*sin(lats) + cos(latc)*cos(lats)
   rnh = rerth*acos(xnh)
   rsh = 0.*rerth*acos(xsh)
   virtempg = 300.  ! isothermal state.
   vg = 0.
   psg = 1.e5
   lnpsg = log(psg)
   call grdtospec(lnpsg,lnpsspec)
   call calc_pressdata(lnpsg)
   do k=1,nlevs
      ! add a barotropic zonal wind perturbation (opp sign in each hemisphere)
      !JOY changed to single tropical jet
      call random_number(rand_vel)
      ug(:,:,k) = 0.8*rand_vel
      !ug(:,:,k) = up*(exp(-(rnh/pertrad)**2)+exp(-(rsh/pertrad)**2))
      call getvrtdivspec(ug(:,:,k),vg(:,:,k),vrtspec(:,k),divspec(:,k),rerth)
      call grdtospec(virtempg(:,:,k),virtempspec(:,k))
   enddo
   topospec = 0.
 end subroutine heldsuarez_ics

 subroutine setdampspec(ndiss,efold,hdif_fac,hdif_fac2,fshk,disspec,diff_prof,dmp_prof)
   ! set hyper-diffusion and upper level rayleigh damping parameters/structure.
   ! if efold and/or ndiss are zero, default GFS parameters are used.
   integer, intent(inout) :: ndiss
   real(r_kind), intent(inout) :: efold,fshk
   real(r_kind), intent(in) :: hdif_fac,hdif_fac2
   real(r_kind), intent(out),dimension(ndimspec) :: disspec
   real(r_kind), dimension(nlevs),intent(out) :: diff_prof,dmp_prof
   integer k
   real(r_kind) dmp_prof1
   ! if ndiss=0, use default value of 8
   if (ndiss == 0) ndiss = 8
   ! if efold <= 0, use GFS defaults.
   if (efold .le. 0) then
      if (ntrunc > 170) then
         efold = 3600./(hdif_fac2*(ntrunc/170.)**4*1.1)
      else if (ntrunc == 126 .or. ntrunc == 170) then
         efold = 1./(hdif_fac2*12.E15/(RERTH**4)*(80.*81.)**2)
      else
         efold = 1./(hdif_fac2*3.E15/(RERTH**4)*(80.*81.)**2)
      end if
      efold = 2.*efold ! mysterious factor 2 in deldifs.f
      !print *,'ndiss,efold =',ndiss,efold
   end if
   if (fshk .le. 0) then
      ! factor to multiply diffusion per scale height.
      fshk = 1.0*hdif_fac
      if (ntrunc > 170) fshk = 2.2*hdif_fac
      if (ntrunc == 126) fshk = 1.5*hdif_fac
   end if
   !slrd0=0.00        ! SIGMA LEVEL AT WHICH TO BEGIN RAYLEIGH MOMTUM DAMPING
   dmp_prof1=1./taustratdamp ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
                      ! ABOVE BEGINNING SIGMA LEVEL FOR RAYLEIGH DAMPING
   dmp_prof = 0.
   diff_prof = sl**(log(1./fshk))
   !print *,'profiles for diffusion and linear damping:'
   !print *,'(level, sigma, diffusion enhancment, linear drag coeff)'
   do k=1,nlevs
      if (sl(k) .lt. slrd0) dmp_prof(k)=dmp_prof1*log(slrd0/sl(k))
   !   print *,k,sl(k),diff_prof(k),dmp_prof(k)
   enddo
   disspec(:) = -(1./efold)*(lap(:)/minval(lap))**(ndiss/2)
   return
 end subroutine setdampspec

end module dyn_init

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
 ntoz,ntclw,pdryini,massfix
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

 contains
 !JOY subroutine to set topography from arbitrary array
 subroutine set_topography(surf_geop) bind(c, name='gfs_set_topography')

    real(r_kind), intent(in), dimension(nlons,nlats) :: surf_geop

    call grdtospec(surf_geop/grav, topospec)
    call spectogrd(grav*topospec, phis)
    call getgrad(grav*topospec, dphisdx, dphisdy, rerth)

 end subroutine set_topography

 subroutine init_dyn() bind(c, name='gfs_init_dynamics')
    integer k
    ! allocate arrays
    ! JOY data array inits will be done from python

    !call init_specdata()
    !call init_griddata()
    
    
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
    !JOY adding read from namelist to eliminate using sigio to read ak,bk
    ak = (/ 200.00000000,566.898010,1290.53296,2210.97900&
    ,3376.51611,4844.03613,6678.60791,8913.76660,&
    11343.6543,13671.4268,15613.5645,16912.1309,&
    17364.6582,16959.9941,15812.9258,14089.5352,&
    11991.4277,9731.80664,7508.53223,5481.14404,&
    3757.14209,2389.20508,1381.52600,701.453003,&
    292.576996,86.4570007,11.6350002,8.99999961E-03,0.00000000 /)

    bk = (/ 0.00000000,0.00000000,0.00000000,0.00000000,&
            0.00000000,0.00000000,0.00000000,3.78680008E-04,&
            3.93276010E-03,1.43262697E-02,3.49500999E-02,&
            6.86749965E-02,0.117417611,0.180667043,0.256084293,&
            0.340293199,0.429195374,0.518456697,0.604054928,&
            0.682746828,0.752347112,0.811784863,0.860974848,&
            0.900580883,0.931749582,0.955872416,0.974401832,&
            0.988725841,1.00000000 /)


    do k=1,nlevs
       dbk(k) = bk(k+1)-bk(k)
       bkl(k) = 0.5*(bk(k+1)+bk(k))
       ck(k)  = ak(k+1)*bk(k)-ak(k)*bk(k+1)
    enddo
         ! compute sigma coordinate quantities (bottom to top).
    do k=1,nlevs+1
        si(nlevs+2-k)= ak(k)/101300.0+bk(k) ! si are now sigmas
    enddo
    do k=1,nlevs
        sl(k) = 0.5*(si(k)+si(k+1))
         !sl(k) = ((si(k)**(rk+1.) - si(k+1)**(rk+1.))/&
         !        ((rk+1.)*(si(k)-si(k+1))))**(1./rk)
    enddo

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
 
 ! JOY remove routines which read initial conditions from file
! subroutine readin_sig(filename,vrtspecin,divspecin,virtempspecin,tracerspecin,lnpsspecin,topospecin,sigheadin)
!    complex(r_kind), intent(out), dimension(ndimspec,nlevs) :: vrtspecin,divspecin,virtempspecin
!    complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: tracerspecin
!    complex(r_kind), intent(out), dimension(ndimspec) :: lnpsspecin,topospecin
!    type(sigio_head), intent(out) :: sigheadin
!    type(sigio_data) sigdata
!    character(*), intent(in) :: filename
!    integer lu,iret,k,nt
!    ! read initial conditions
!    lu = 7
!    call sigio_srohdc(lu,trim(filename),sigheadin,sigdata,iret)
!    if (iret .ne. 0) then
!      print *,'error reading ',trim(filename),iret
!      stop
!    endif
!    ! convert spectral arrays to double precision complex,
!    ! re-normalize coefficients.
!    call copyspecin(sigdata%ps, lnpsspecin)
!    call copyspecin(sigdata%hs, topospecin)
!    do k=1,nlevs
!       call copyspecin(sigdata%z(:,k),vrtspecin(:,k))
!       call copyspecin(sigdata%d(:,k),divspecin(:,k))
!       call copyspecin(sigdata%t(:,k),virtempspecin(:,k))
!       do nt=1,ntrac
!          call copyspecin(sigdata%q(:,k,nt),tracerspecin(:,k,nt))
!       enddo
!    enddo
!    call sigio_axdata(sigdata,iret)
!    call sigio_sclose(lu,iret)
! end subroutine readin_sig
!
! subroutine wrtout_sig(fh,filename)
!    ! write out spectral data
!    ! (this probably belongs in a separate module)
!    real(r_kind), dimension(nlons,nlats) :: psg 
!    complex(r_kind), dimension(ndimspec) :: lnpsspec_tmp
!    type(sigio_data) sigdata
!    real(r_kind), intent(in) :: fh
!    character(len=500), intent(in) :: filename
!    integer k,lu,iret,lu2,nt
!    ! convert to ln(ps) in Pa.
!    call spectogrd(lnpsspec, psg)
!    psg = exp(psg)/1000. ! convert to cb
!    call grdtospec(log(psg), lnpsspec_tmp) ! back to spectral.
!    !print *,'fhour =',int(fh),', min/max psg = ',minval(10.*psg),maxval(10.*psg)
!    lu = 7; lu2 = 8
!    call sigio_aldata(sighead,sigdata,iret)
!    if (iret .ne. 0) then
!      print *,'error allocating sigdata',iret
!      stop
!    endif
!    sighead%fhour = fh
!    sighead%pdryini = pdryini/1000. ! convert to cb from Pa.
!    call copyspecout(lnpsspec_tmp, sigdata%ps)
!    call copyspecout(topospec, sigdata%hs)
!    do k=1,nlevs
!       call copyspecout(vrtspec(:,k), sigdata%z(:,k))
!       call copyspecout(divspec(:,k), sigdata%d(:,k))
!       call copyspecout(virtempspec(:,k), sigdata%t(:,k))
!       do nt=1,ntrac
!         call copyspecout(tracerspec(:,k,nt), sigdata%q(:,k,nt))
!       enddo
!    enddo
!    call sigio_swohdc(lu2,filename,sighead,sigdata,iret)
!    if (iret .ne. 0) then
!      print *,'error writing ',trim(filename),iret
!      stop
!    endif
!    call sigio_axdata(sigdata,iret)
!    call sigio_sclose(lu,iret)
!    call sigio_sclose(lu2,iret)
! end subroutine wrtout_sig

! subroutine wrtout_gfsgrb(fh,filename)
!    use gfsio_module, only: gfsio_gfile, gfsio_open, gfsio_writerecvw34, &
!                            gfsio_init, gfsio_close
!    ! write out gfsio grib data
!    ! (this probably belongs in a separate module)
!    type(gfsio_gfile)    :: gfile
!    real(r_kind), intent(in) :: fh
!    character(len=500), intent(in) :: filename
!    integer k,iret,nt,reclev(2+nlevs*8)
!    character(len=8) recname(2+nlevs*8)
!    character(len=16) reclevtyp(2+nlevs*9)
!    real(4) pdryini4,fhour4,tmpg(nlons,nlats)
!
!    pdryini4 = pdryini
!    fhour4    = fh
!    recname(1)   = 'hgt'
!    reclevtyp(1) = 'sfc'
!    reclev(1)    = 1
!    recname(2)   = 'pres'
!    reclevtyp(2) = 'sfc'
!    reclev(2)    = 1
!    do k=1,nlevs
!      recname(k+2)            = 'pres'
!      reclevtyp(k+2)          = 'layer'
!      reclev(k+2)             = k
!      recname(k+2+nlevs)      = 'dpres'
!      reclevtyp(k+2+nlevs)    = 'layer'
!      reclev(k+2+nlevs)       = k
!      recname(k+2+nlevs*2)    = 'tmp'
!      reclevtyp(k+2+nlevs*2)  = 'layer'
!      reclev(k+2+nlevs*2)     = k
!      recname(k+2+nlevs*3)    = 'ugrd'
!      reclevtyp(k+2+nlevs*3)  = 'layer'
!      reclev(k+2+nlevs*3)     = k
!      recname(k+2+nlevs*4)    = 'vgrd'
!      reclevtyp(k+2+nlevs*4)  = 'layer'
!      reclev(k+2+nlevs*4)     = k
!      recname(k+2+nlevs*5)    = 'spfh'
!      reclevtyp(k+2+nlevs*5)  = 'layer'
!      reclev(k+2+nlevs*5)     = k
!      recname(k+2+nlevs*6)    = 'o3mr'
!      reclevtyp(k+2+nlevs*6)  = 'layer'
!      reclev(k+2+nlevs*6)     = k
!      recname(k+2+nlevs*7)    = 'clwmr'
!      reclevtyp(k+2+nlevs*7)  = 'layer'
!      reclev(k+2+nlevs*7)     = k
!      recname(k+2+nlevs*8)    = 'vvel'
!      reclevtyp(k+2+nlevs*8)  = 'layer'
!      reclev(k+2+nlevs*8)     = k
!    enddo
!
!    !print *,' calling gfsio_open idate=',idate_start,' fhour=',fhour4
! 
!    call gfsio_init(iret)
!    call gfsio_open(gfile,trim(filename),'write',iret,&
!         version=sighead%ivs,fhour=fhour4,idate=idate_start,nrec=2+nlevs*9,&
!         latb=nlats,lonb=nlons,levs=nlevs,jcap=ntrunc,itrun=sighead%itrun,&
!         iorder=sighead%iorder,irealf=sighead%irealf,igen=sighead%igen,latf=nlats,lonf=nlons,&
!         latr=nlats,lonr=nlons,ntrac=ntrac,icen2=sighead%icen2,iens=sighead%iens,&
!         idpp=sighead%idpp,idsl=sighead%idsl,idvc=sighead%idvc,idvm=sighead%idvm,&
!         idvt=sighead%idvt,idrun=sighead%idrun,&
!         idusr=sighead%idusr,pdryini=pdryini4,ncldt=sighead%ncldt,nvcoord=sighead%nvcoord,&
!         vcoord=sighead%vcoord,recname=recname,reclevtyp=reclevtyp,&
!         reclev=reclev,Cpi=sighead%cpi,Ri=sighead%ri)
!
!    call twodtooned(phis/grav,tmpg)
!    call gfsio_writerecvw34(gfile,'hgt','sfc',1,tmpg,iret)
!    call twodtooned(psg,tmpg)
!    call gfsio_writerecvw34(gfile,'pres','sfc',1,tmpg,iret)
!    do k=1,nlevs
!       call twodtooned(prs(:,:,k),tmpg)
!       call gfsio_writerecvw34(gfile,'pres','layer',k,&
!                               tmpg, iret)
!    enddo
!    do k=1,nlevs
!       call twodtooned(dpk(:,:,nlevs-k+1),tmpg)
!       call gfsio_writerecvw34(gfile,'dpres','layer',k, &
!                               tmpg, iret)
!    enddo
!    do k=1,nlevs
!       call twodtooned(virtempg(:,:,k),tmpg)
!       tmpg = tmpg/(1.+fv*tracerg(:,:,k,1))
!       call gfsio_writerecvw34(gfile,'tmp','layer',k,&
!                               tmpg, iret)
!    enddo
!    do k=1,nlevs
!       call twodtooned(ug(:,:,k),tmpg)
!       call gfsio_writerecvw34(gfile,'ugrd','layer',k, &
!                               tmpg, iret)
!    enddo
!    do k=1,nlevs
!       call twodtooned(vg(:,:,k),tmpg)
!       call gfsio_writerecvw34(gfile,'vgrd','layer',k,&
!                               tmpg, iret)
!    enddo
!    do k=1,nlevs
!       call twodtooned(tracerg(:,:,k,1),tmpg)
!       call gfsio_writerecvw34(gfile,'spfh','layer',k,&
!                               tmpg, iret)
!    enddo
!    if (ntoz .gt. 0) then
!       do k=1,nlevs
!          call twodtooned(tracerg(:,:,k,ntoz),tmpg)
!          call gfsio_writerecvw34(gfile,'o3mr','layer',k,&
!                               tmpg, iret)
!       enddo
!    endif
!    if (ntclw .gt. 0) then
!       do k=1,nlevs
!          call twodtooned(tracerg(:,:,k,ntclw),tmpg)
!          call gfsio_writerecvw34(gfile,'clwmr','layer',k,&
!                                  tmpg, iret)
!       enddo
!    endif
!    do k=1,nlevs
!       call twodtooned(dlnpdtg(:,:,k),tmpg)
!       !tmpg = tmpg*prs(:,:,k)
!       tmpg = tmpg*0.5*(pk(:,:,nlevs-k+2)+pk(:,:,nlevs-k+1))
!       call gfsio_writerecvw34(gfile,'vvel','layer',k,&
!                               tmpg, iret, precision=6)
!    enddo

!    call gfsio_close(gfile,iret)

! end subroutine wrtout_gfsgrb

 subroutine twodtooned(data2,data1)
   real(r_kind), intent(in) :: data2(nlons,nlats)
   real(4), intent(out) :: data1(nlons*nlats)
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


 subroutine heldsuarez_ics()
   ! jablonowski and williamson (2006, QJR, p. 2943, doi: 10.1256/qj.06.12)
   ! initial perturbation on an isothermal state.
   real(r_kind), dimension(nlons,nlats) :: rnh,xnh,rsh,xsh, rand_vel
   real(r_kind) :: lonc,latc,up,pertrad
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
   real(r_kind) slrd0,dmp_prof1
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
      print *,'ndiss,efold =',ndiss,efold
   end if
   if (fshk .le. 0) then
      ! factor to multiply diffusion per scale height.
      fshk = 1.0*hdif_fac
      if (ntrunc > 170) fshk = 2.2*hdif_fac
      if (ntrunc == 126) fshk = 1.5*hdif_fac
   end if
   slrd0=0.01        ! SIGMA LEVEL AT WHICH TO BEGIN RAYLEIGH MOMTUM DAMPING
   dmp_prof1=1./taustratdamp ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
                      ! ABOVE BEGINNING SIGMA LEVEL FOR RAYLEIGH DAMPING
   dmp_prof = 0.
   diff_prof = sl**(log(1./fshk))
   print *,'profiles for diffusion and linear damping:'
   print *,'(level, sigma, diffusion enhancment, linear drag coeff)'
   do k=1,nlevs
      if (sl(k) .lt. slrd0) dmp_prof(k)=dmp_prof1*log(slrd0/sl(k))
      print *,k,sl(k),diff_prof(k),dmp_prof(k)
   enddo
   disspec(:) = -(1./efold)*(lap(:)/minval(lap))**(ndiss/2)
   return
 end subroutine setdampspec

end module dyn_init

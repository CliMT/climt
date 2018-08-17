module params
! holds model parameters
! Public subroutines:
! read_namelist: read namelist
 use kinds, only: r_kind,r_double
! JOY removing sigio dependency
! use sigio_module, only: sigio_head, sigio_srhead, sigio_sropen, sigio_sclose

!JOY adding iso_c_binding required for relevant variables
use iso_c_binding, only: c_double,c_int,c_bool
 implicit none
 private

 public :: initfile,sfcinitfile,dt,ntmax,ndimspec,nlons,nlats,&
 tstart,ndiss,efold,nlevs,ntrunc,dry,explicit,heldsuarez,dcmip,&
 ntout,fhdfi,fhzer,idate_start,adiabatic,hdif_fac,hdif_fac2,fshk,ntrac,ntoz,ntclw,&
 pdryini, toa_pressure, massfix,timestepsperhr,ncw,taustratdamp,polar_opt,ntdfi,gfsio_out,sigio_out,&
! gfs phys parameters.
 nmtvr,fhlwr,fhswr,ictm,isol,ico2,iaer,ialb,iems,isubc_sw,isubc_lw,&
 iovr_sw,iovr_lw,newsas,ras,sashal,num_p3d,num_p2d,crick_proof,ccnorm,&
 norad_precip,crtrh,cdmbgwd,ccwf,dlqf,ctei_rm,psautco,prautco,evpco,wminco,flgmin,&
 old_monin,cnvgwd,mom4ice,shal_cnv,cal_pre,trans_trac,nst_fcst,moist_adj,mstrat,&
 ngptc,pre_rad,bkgd_vdif_m,bkgd_vdif_s,bkgd_vdif_h,gloopb_filter,&
! iau parameters
 iau,iaufiles_fg,iaufiles_anl,iaufhrs,iau_delthrs,&
! vorticity confinement parameters
 vcamp,svc,svc_tau,svc_lscale,iseed_svc,&
! stochastic physics tendency parameters
 sppt,sppt_tau,sppt_lscale,iseed_sppt, &
! additive stochastic humidity perturbations
 shum,shum_tau,shum_lscale,iseed_shum,clipsupersat,&
! additive noise perts
 addnoise,addnoise_tau,addnoise_lscale,addnoise_vfilt,iseed_addnoise,&
 addnoise_kenorm, addnoise_dissfact,addnoise_vrtonly


 !! JOY adding bind(c) arguments wherever needed

 character(len=500) :: initfile ! init cond filename
 character(len=500) :: sfcinitfile ! surface init cond filename
! integer(c_int)     :: fhmax ! hours to run
! integer(c_int)     :: fhout ! interval for IO
 integer            :: fhzer ! interval to zero accumulated arrays
! half window length (hrs) for digital filter launch (=0 mean no dfi)
 integer            :: fhdfi=0
 real(c_double)     :: deltim=0    ! namelist input time step (secs)

 real(c_double)     :: dt    ! time step (secs) (=deltim or 3600/timestepsperhr)
! bind(c) :: fhmax,fhout
 bind(c) :: deltim,dt

 real(c_double) :: pdryini ! initial dry ps
 bind(c) :: pdryini

 real(c_double) :: toa_pressure !top of model pressure
 bind(c) :: toa_pressure

 logical(c_bool)  :: massfix=.true. ! apply dry mass 'fixer'
 bind(c) :: massfix

 logical    :: gfsio_out=.false. ! write out 'gfsio' grib 1 files.
 logical    :: sigio_out=.true. ! write out 'sigma' spectral binary files.
 integer    :: ntmax ! time steps to run
 integer    :: ntdfi ! number of time steps in dfi window is 2*ntdfi+1

 integer(c_int)    :: nlons=1 ! number of longitudes on grid
 integer(c_int)    :: nlats=1 ! number of latitudes on grid
 integer(c_int)    :: nlevs=1 ! number of levels on grid
 integer(c_int)    :: ntrunc=1 ! spectral truncation
 integer(c_int)    :: ndimspec=1 ! spectral array dimension
 

! JOY removing sigio dependency
! type(sigio_head),save  :: sighead ! header struct from initfile
 logical(c_bool)    :: dry = .false. ! no moisture, cloud condensate or ozone.
 logical(c_bool)    :: adiabatic = .false. ! don't call physics
 ! held-suarez forcing
 logical(c_bool)    :: heldsuarez = .false.
 bind(c) :: dry,adiabatic,heldsuarez

 logical    :: iau = .false. ! iau forcing included
 integer    :: iau_delthrs = 6 ! iau time interval (to scale increments)
 character(len=120), dimension(7) ::  iaufiles_fg,iaufiles_anl
 real(r_kind), dimension(7) :: iaufhrs
 ! dcmip test cases
 integer    :: dcmip = -1 ! 4x for baroclinic wave, 5x for tropical cyclone.
 ! 3-stage, 2nd order RK time differencing is used.
 ! explicit=.true. means use explicit time differencing
 ! if .false., semi-implicit or additive RK (Ascher, Spiteri and Ruth 1997
 ! DOI:10.1016/S0168-9274(97)00056-1) is used.
 logical    :: explicit = .false. ! use semi-implicit RK
 ! starting forecast time in seconds (read in from initfile)
 real(r_kind) :: tstart
 integer    :: idate_start(4) ! starting date (hr,month,day,year)
 integer    :: ntout ! time step interval for IO
 integer    :: ndiss=0 ! hyperdiffusion order (0 for GFS defaults)
 real(r_double) :: polar_opt=1.e-10 ! polar optimization threshold for transforms
 ! efolding scale for smallest resolvable wave (0 for GFS defaults)
 real(r_kind) :: efold=0. 
 real(r_kind) :: hdif_fac=1.0 ! multiplier for height-dep part of hyper-diff
 real(r_kind) :: hdif_fac2=1.0 ! multiplier to increase hyper-diffusion
 ! amplitude of vertically varying part of hyper-diff (1 means no variation,
 ! zero gives GFS resolution dependent defaults)
 real(r_kind) :: fshk=0 
 integer(c_int) :: ntrac=3 ! number of tracers (including specific humidity)

   
!JOY bind the grid size variables   
 bind(c) :: nlons, nlats, nlevs, ntrunc, ndimspec, ntrac

 integer :: ntoz=2 ! ozone tracer number
 integer :: ntclw=3 ! cloud condensate tracer number
 integer :: nmtvr=14 ! number of fields in mtnvar file.
 integer :: ngptc=24 ! pass this many profiles to gbphys in one call.
 real(r_kind) :: taustratdamp=2.*86400. ! extra linear drag near top of model
! parameters relevant for GFS physics
 ! interval in hours to call long-wave radiation (0 means every time step)
 real(r_kind) :: fhlwr=0 
 ! interval in hours to call short-wave radiation (0 means every time step)
 real(r_kind) :: fhswr=0 
 ! ictm controls source for controls source for time sensitive external data (e.g. CO2,
 ! solcon, aerosols, etc)
 integer :: ictm=0 ! use data at initial cond time, or latest available.
 ! isol controls solar constant data source.
 integer :: isol=0 ! use prescribed solar constant
 ! ico2 controls co2 data source for radiation
 integer :: ico2=0 ! prescribed global mean co2
 ! iaer controls aerosols scheme selection
 integer :: iaer=1 ! default aerosol
 ! ialb controls surface albedo for sw radiation
 integer :: ialb=0 ! use climo albedo based on sfc type
 ! iems controls surface emissivity and sfc air/ground temp for lw radiation
 integer :: iems=0 ! used fixed value of 1.0
 ! isubc_sw,lw controls sub-column cloud approximation in radiation
 integer :: isubc_sw=0 ! sw clouds without sub-column approximation
 integer :: isubc_lw=0 ! lw clouds without sub-column approximation
 ! iovr_sw,lw controls cloud overlapping method in radiation
 integer :: iovr_sw=1 ! max-random overlap clouds
 integer :: iovr_lw=1 ! max-random overlap clouds
 logical :: newsas=.true. ! 'new. SAS convection scheme
 logical :: ras=.false. ! RAS convection scheme.
 logical :: sashal=.true. ! 'new' mass-sflux based shallow convection
 logical :: mstrat=.false. ! flag for moorthi approach for stratus
 logical :: pre_rad=.false. ! flag for debugging in gbphys
 integer :: num_p3d=4 ! # of 3d microphysics fields (value for Zhao microphysics)
 integer :: num_p2d=3 ! # of 2d microphysics fields (value for Zhao)
 logical :: crick_proof = .false.
 logical :: ccnorm = .false.
 logical :: norad_precip = .false. ! only used for Ferrier microphysics
 ! range of droplet number concentrations for Ferrier scheme.
 integer :: ncw(2)=(/50,150/)
 ! critical relative humidity at the surface, PBL and atmosphere top   
 real(r_kind) :: crtrh(3)=(/0.85,0.85,0.85/)
 ! multiplication factors for cdmb and gwd  (Mtn Blking and GWD tuning factors)
 real(r_kind) :: cdmbgwd(2)=(/1.0,1.0/)
 ! RAS convection mult factor for critical cloud work function.
 real(r_kind) :: ccwf(2)=(/1.0,1.0/)
 ! factor for cloud condensate detrainment from cloud edges in RAS.
 real(r_kind) :: dlqf(2)=(/0.0,0.0/)
 ! critical cloud top entrainment instability criteria (for mstrat=.true.)
 real(r_kind) :: ctei_rm(2)=(/10.0,10.0/)
 ! auto conversion coeff from ice to snow for Zhao microphysics
 real(r_kind) :: psautco(2)=(/4.0e-4,4.0e-4/)
 ! auto conversion coeff from cloud to rain for Zhao microphysics
 real(r_kind) :: prautco(2)=(/1.0e-4,1.0e-4/)
 real(r_kind) :: evpco=2.e-5 ! Zhao scheme evap coefficient for lg-scale rain
 real(r_kind) :: wminco(2)=(/1.e-5,1.e-5/) !  water and ice minimum threshold for Zhao 
 real(r_kind) :: flgmin(2)=(/0.2,0.2/) ! (Ferrier only) range of minimum large ice fraction
 real(r_kind) :: vcamp=0. ! vorticity confinement amplitude
 real(r_kind) :: svc=0.   ! stochastic vorticity confinement amplitude
 real(r_kind) :: svc_tau=0.      ! stochastic vorticity confinement time scale
 real(r_kind) :: svc_lscale=0.   ! stochastic vorticity confinement length scale
 integer :: iseed_svc=0 ! random seed for stochastic vc (zero means use clock)
 real(r_kind) :: sppt=0.  ! stochastic physics tendency amplitude
 real(r_kind) :: sppt_tau=0.  ! stochastic physics tendency time scale
 real(r_kind) :: sppt_lscale=0.  ! stochastic dynamics tendency length scale
 integer :: iseed_sppt=0 ! random seed for sppt (0 means use system clock)
 real(r_kind) :: shum=0.  ! stochastic humidity pert amplitude
 real(r_kind) :: shum_tau=0.  ! stochastic humidity pert time scale
 real(r_kind) :: shum_lscale=0.  ! stochastic humidity pert length scale
 integer :: iseed_shum=0 ! random seed for stochastic humid pert (0 means use system clock)
 logical :: addnoise_kenorm=.false. ! add noise in ke norm (instead of psi norm)
 logical :: addnoise_dissfact=.false. ! modulate noise amplitude with ke diss
 logical :: addnoise_vrtonly=.true. ! include balanced t noise (if false)
 real(r_kind) :: addnoise=0.  ! add noise pert amplitude
 real(r_kind) :: addnoise_tau=0.  ! add noise pert time scale
 real(r_kind) :: addnoise_lscale=0.  ! add noise pert length scale
 integer :: iseed_addnoise=0 ! random seed for add noise pert (0 means use system clock)
 integer :: addnoise_vfilt=0 ! add noise 1-2-1 vertical filter parameter
 logical :: old_monin = .false. ! flag for old Monin-Obhukov surface layer
 logical :: cnvgwd = .false. ! flag for convective gravity wave drag
 logical :: mom4ice = .false. ! flag for MOM4 sea-ice scheme
 logical :: shal_cnv = .true. ! use shallow convection?
 logical :: cal_pre = .false. ! true for huiya's precip type algorithm
 logical :: trans_trac = .true. ! convective transport of tracers? (RAS only)
 integer :: nst_fcst=0 ! 0 - AM only, 1 - uncoupled, 2 - coupled
 logical :: moist_adj = .false. 
 logical :: gloopb_filter = .true. ! apply spectral filter to physics tendencies
 ! make sure stochastic perts don't create neg or supersat humidities.
 logical :: clipsupersat=.false. 

 real(r_kind) :: bkgd_vdif_m = 3.0 ! background vertical diffusion for momentum
 real(r_kind) :: bkgd_vdif_h = 1.0 ! background vertical diffusion for heat, q
 real(r_kind) :: bkgd_vdif_s = 0.2 ! sigma threshold for background mom. diffusn 
 ! if dt not given, but timestepsperhr is, dt=3600/timestepsperhr
 real(r_double) :: timestepsperhr = -1
 !JOY added to remove sigio dependency
 integer(c_int) :: ntracin = 0
 bind(c) :: ntracin

contains

    subroutine set_dry_pressure(py_pdryini)bind(c, name='gfs_set_dry_pressure')

        real(c_double), intent(in):: py_pdryini

        pdryini = py_pdryini

    end subroutine

    subroutine set_model_top_pressure(py_ptoa)bind(c, name='gfs_set_model_top_pressure')

        real(c_double), intent(in):: py_ptoa

        toa_pressure = py_ptoa

    end subroutine

    subroutine gfs_set_time_step(py_dt)bind(c, name='gfs_set_time_step')

        real(c_double), intent(in):: py_dt

        dt = py_dt
        deltim = py_dt

    end subroutine

    subroutine gfs_set_model_dimensions(py_nlats, py_nlons, py_nlevs, py_ntrunc,&
                    py_ndimspec, py_ntrac)bind(c,name='gfs_set_model_dimensions')

        integer(c_int), intent(in):: py_nlats, py_nlons, py_nlevs, py_ntrunc, py_ndimspec, py_ntrac

        nlats = py_nlats
        nlons = py_nlons
        nlevs = py_nlevs
        ntrunc = py_ntrunc
        ndimspec = py_ndimspec
        ntrac = py_ntrac

    end subroutine
end module params

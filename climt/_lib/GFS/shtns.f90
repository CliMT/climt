      module shtns

! fortran95 interface to the SHTNS library
! (http://users.isterre.fr/nschaeff/SHTns/index.html).
! this version hardwired for gaussian grids.

! The available routines are:

! SUBROUTINE GRDTOSPEC(datagrid,dataspec):
! converts input gridded data array (datagrid) to complex spectral
! coefficients (dataspec).

! SUBROUTINE SPECTOGRD(dataspec,datagrid):
! converts input spectral coefficient array (dataspec) to a grid (datagrid).

! SUBROUTINE GETUV(vrtspec,divspec,ugrid,vgrid,rsphere):
! given input spectral coefficients of vorticity and divergence
! (vrtspec,divspec) calculates gridded winds (ugrid,vgrid).
! rsphere is the radius of the sphere.

! SUBROUTINE GETVRTDIVSPEC(ugrid,vgrid,vrtspec,divspec,rsphere):
! given input gridded winds (ugrid,vgrid) calculates spectral coefficients
! of vorticity and divergence (vrtspec,divspec). rsphere is the radius of
! the sphere.

! SUBROUTINE GETGRAD(dataspec,gradx,grady,rsphere):
! calculates gridded vector gradient (gradx,grady) given input
! spectral coefficients. rsphere is the radius of the sphere.

! SUBROUTINE SHTNS_CLEANUP():
! Garbage collection. Deallocates all memory.
! Call this when you are done with the module to free up memory.

! SUBROUTINE SHTNS_INIT(NLONS,NLATS,NTRUNC):
! initialize transform machinery for a given grid and spectral
! truncation.  Called automatically by other routines if
! grid (as described by input arrays) changes.

! module data:
! lats,lons: 2-d arrays with lat/lon values (radians) of transform grid.
! nlm: number of complex spectral coeffs.
! current_nlon, current_nlat: current values of nlons,nlats
! gauwts: 1-d array (length current_nlat) of gaussian quadrature weights.
! current_ntrunc: current spectral truncation
! degree: array with degree (n) for each spectral coefficient.
! order: array with order (m) for each spectral coefficient.
! lap: lapacian operator (-n*(n+1))
! invlap: inverse of laplacian operator.

! All of these routines use triangular truncation.
! The spherical harmonic coefficients are normalized using an
! 'orthonormalized' norm
! (http://en.wikipedia.org/wiki/Spherical_harmonics#Conventions).

      use kinds, only: r_kind, default_real, r_double
      use iso_c_binding
      implicit none
      private
      public :: shtns_init,grdtospec,spectogrd,getuv,getvrtdivspec,&
                getgrad,shtns_destroy
      public :: gauwts, lats, lons, nlm, degree, order, lap, invlap,&
                current_nlon, current_nlat, current_ntrunc, areawts,&
            print_bool
      INTEGER, PARAMETER :: SHT_NATIVE_LAYOUT=0
      INTEGER, PARAMETER :: SHT_THETA_CONTIGUOUS=256
      INTEGER, PARAMETER :: SHT_PHI_CONTIGUOUS=512
      INTEGER, PARAMETER :: SHT_NO_CS_PHASE=1024
      INTEGER, PARAMETER :: SHT_REAL_NORM=2048
      INTEGER, PARAMETER :: SHT_ORTHONORMAL=0
      INTEGER, PARAMETER :: SHT_FOURPI=1
      INTEGER, PARAMETER :: SHT_SCHMIDT=2
      INTEGER, PARAMETER :: SHT_GAUSS=0
      INTEGER, PARAMETER :: SHT_AUTO=1
      INTEGER, PARAMETER :: SHT_REG_FAST=2
      INTEGER, PARAMETER :: SHT_REG_DCT=3
      INTEGER, PARAMETER :: SHT_QUICK_INIT=4
      INTEGER, PARAMETER :: SHT_REG_POLES=5
      INTEGER, PARAMETER :: SHT_GAUSS_FLY=6
      REAL(r_double), PARAMETER :: SHT_DEFAULT_POLAR_OPT=1.d-10
      INTEGER, SAVE      :: current_nlon = -1
      INTEGER, SAVE      :: current_nlat = -1
      INTEGER, SAVE      :: current_ntrunc = -1
      INTEGER            :: nlm
      INTEGER :: print_bool
! arrays allocated when nlon or nlat or ntrunc change.
      REAL(r_kind), DIMENSION(:), pointer :: lap, invlap, gauwts
      REAL(r_kind), DIMENSION(:,:), pointer :: lats
      REAL(r_kind), DIMENSION(:,:), pointer :: lons, areawts
      INTEGER, DIMENSION(:), pointer :: degree, order
      real(r_double) :: popt ! polar optimization thresh
      integer :: nth ! number of threads to use

      contains

      subroutine sht_assign_arrays(py_lap, py_invlap, py_gauwts,&
              py_lons, py_lats, py_areawts, py_degree,&
              py_order, py_nlons, py_nlats, py_nlm)&
              bind(c, name='gfs_assign_sht_arrays')

          REAL(c_double), DIMENSION(py_nlm), target :: py_lap, py_invlap
          REAL(c_double), DIMENSION(py_nlats), target :: py_gauwts
          REAL(c_double), DIMENSION(py_nlons,py_nlats), target :: py_lats
          REAL(c_double), DIMENSION(py_nlons,py_nlats), target :: py_lons
          REAL(c_double), DIMENSION(py_nlons,py_nlats), target :: py_areawts
          INTEGER(c_int), DIMENSION(py_nlm), target :: py_degree, py_order
          integer(c_int), intent(in) :: py_nlons, py_nlats, py_nlm

          lap => py_lap
          invlap => py_invlap
          gauwts => py_gauwts
          lats => py_lats
          lons => py_lons
          areawts => py_areawts
          degree => py_degree
          order => py_order

      end subroutine sht_assign_arrays

      !subroutine get_lon_lat(longitudes, latitudes) &
      !        bind(c,name='gfs_get_lon_lat')

      !real(c_double), intent(out), dimension(current_nlon,current_nlat)&
      !          :: longitudes, latitudes

      !if(allocated(lats)) then
      !    latitudes(:,:) = lats(:,:)
      !else
      !    print *, 'latitudes not allocated'
      !endif

      !if(allocated(lons)) then
      !    longitudes(:,:) = lons(:,:)
      !else
      !    print *, 'longitudes not allocated'
      !endif

      !end subroutine get_lon_lat

      subroutine shtns_init(nlon,nlat,ntrunc,nthreads,polar_opt)
! initialize library, allocate arrays.
      integer, intent(in) :: nlon,nlat,ntrunc
      integer, intent(in), optional :: nthreads
      real(r_double), intent(in), optional :: polar_opt
      real(r_double), dimension(:), allocatable :: lats1,gauwts1
      real(r_double) pi
      integer m,n,i,j
      print_bool = 0
      if (present(nthreads)) then
         nth = nthreads
      else
         nth = 1
      endif
      if (present(polar_opt)) then
         popt = polar_opt
      else
         popt = SHT_DEFAULT_POLAR_OPT
      endif
      call shtns_destroy()
      call shtns_use_threads(0) ! JOY auto compute number of openmp threads.
      !print *, 'In set_size', ntrunc
      call shtns_set_size(ntrunc,ntrunc,1,SHT_ORTHONORMAL+SHT_NO_CS_PHASE)
      !print *, 'In precompute'
      !print *, nlat, nlon
      call shtns_precompute(SHT_GAUSS_FLY,SHT_PHI_CONTIGUOUS,popt,nlat,nlon)
      call shtns_calc_nlm(nlm,ntrunc,ntrunc,1)
      current_nlat = nlat
      current_nlon = nlon
      current_ntrunc = ntrunc
      if (nlm .ne. (ntrunc+1)*(ntrunc+2)/2) then
         print *,'error: nlm not what expected',nlm,ntrunc
         stop
      endif
      !if (.not. allocated(lats)) allocate(lats(nlon,nlat))
      !if (.not. allocated(gauwts)) allocate(gauwts(nlat))
      !if (.not. allocated(lons)) allocate(lons(nlon,nlat))
      allocate(lats1(nlat))
      allocate(gauwts1(nlat/2))
      call shtns_cos_array(lats1)
      call shtns_gauss_wts(gauwts1)
      ! lats is 2d-array (lon,lat), lats1 is 1d
      do i=1,nlon
         lats(i,:) = asin(lats1)
      enddo
      do j=1,nlat/2
         gauwts(j) = gauwts1(j)
         gauwts(nlat-j+1) = gauwts1(j)
      enddo
      !allocate(areawts(nlon,nlat))
      do i=1,nlon
         areawts(i,:) = gauwts(:)
      enddo
      areawts = areawts/sum(areawts)
      pi = 4.*atan(1.0)
      do j=1,nlat
      do i=1,nlon
         lons(i,j) = 2.*pi*real(i-1)/nlon
      enddo
      enddo
      deallocate(lats1,gauwts1)
      !if (.not. allocated(degree)) allocate(degree(nlm))
      !if (.not. allocated(order)) allocate(order(nlm))
      degree = (/((n,n=m,ntrunc),m=0,ntrunc)/)
      order = (/((m,n=m,ntrunc),m=0,ntrunc)/)
      !if (.not. allocated(lap)) allocate(lap(nlm))
      !if (.not. allocated(invlap)) allocate(invlap(nlm))
      lap = -degree*(degree+1.0)
      invlap = 0
      invlap(2:) = 1./lap(2:)
      end subroutine shtns_init

      subroutine shtns_destroy()
! deallocate arrays.
      call shtns_reset()
      !if (allocated(lats)) deallocate(lats)
      !if (allocated(gauwts)) deallocate(gauwts)
      !if (allocated(areawts)) deallocate(areawts)
      !if (allocated(lons)) deallocate(lons)
      !if (allocated(degree)) deallocate(degree)
      !if (allocated(order)) deallocate(order)
      !if (allocated(lap)) deallocate(lap)
      !if (allocated(invlap)) deallocate(invlap)
      !nullify(lats,gauwts,areawts,lons,degree,order,lap,invlap)
      end subroutine shtns_destroy

      subroutine grdtospec(datagrid,dataspec)
! converts gridded input array (datagrid) to complex spectral coefficients
! (dataspec).
      real(c_double), dimension(:,:), intent(in) :: datagrid
      complex(r_kind), dimension(:), intent(out) :: dataspec
      real(r_double), dimension(:,:), allocatable :: datagrid_tmp
      complex(r_double), dimension(:), allocatable :: dataspec_tmp
      integer nlon,nlat,ntrunc
      ! compute array dimensions and infer truncation limit
      ! from size of input arrays
      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1

      if (nlon .ne. current_nlon .or. &
          nlat .ne. current_nlat .or. &
          ntrunc .ne. current_ntrunc) then
          call shtns_init(nlon,nlat,ntrunc,nth,popt)
      end if
      if (default_real .ne. 2) then
         ! if inputs not double precision convert,
         ! then convert back.
         !if (print_bool == 1) then
         !   print *, datagrid(4,1:3)
         !   print *, 'grdtospec not default_real: ', loc(datagrid)
         !   print *, ''
         !end if

         allocate(datagrid_tmp(nlon,nlat))
         allocate(dataspec_tmp(nlm))
         datagrid_tmp = datagrid
         call shtns_spat_to_sh(datagrid_tmp, dataspec_tmp)
         dataspec = dataspec_tmp
         deallocate(datagrid_tmp,dataspec_tmp)
      else
         call shtns_spat_to_sh(datagrid, dataspec)
      endif
      end subroutine grdtospec

      subroutine spectogrd(dataspec,datagrid)
! converts complex spectral coefficients (dataspec) to
! gridded data array (datagrid).
      real(r_kind), dimension(:,:), intent(out) :: datagrid
      complex(r_kind), dimension(:), intent(in) :: dataspec
      real(r_double), dimension(:,:), allocatable :: datagrid_tmp
      complex(r_double), dimension(:), allocatable :: dataspec_tmp
      integer nlon,nlat,ntrunc
      ! compute array dimensions and infer truncation limit
      ! from size of input arrays
      nlon = size(datagrid,1)
      nlat = size(datagrid,2)
      ntrunc = nint((-1.+sqrt(1+8*float(size(dataspec))))/2.)-1
      if (nlon .ne. current_nlon .or. &
          nlat .ne. current_nlat .or. &
          ntrunc .ne. current_ntrunc) then
          call shtns_init(nlon,nlat,ntrunc,nth,popt)
      end if
      if (default_real .ne. 2) then
         allocate(datagrid_tmp(nlon,nlat))
         allocate(dataspec_tmp(nlm))
         dataspec_tmp = dataspec
         call shtns_sh_to_spat(dataspec_tmp, datagrid_tmp)
         datagrid = datagrid_tmp
         deallocate(datagrid_tmp,dataspec_tmp)
      else
         call shtns_sh_to_spat(dataspec, datagrid)
      endif
      end subroutine spectogrd

      subroutine getuv(vrtspec,divspec,ugrid,vgrid,rsphere)
! given spectral coefficients of vorticity and divergence
! (vrtspec,divspec) compute gridded winds (ugrid,vgrid).
! rsphere is the radius of the sphere in meters.
      real(r_kind), intent(in) :: rsphere
      real(r_kind), dimension(:,:), intent(out) :: ugrid,vgrid
      complex(r_kind), dimension(:), intent(in) :: vrtspec,divspec
      real(r_double), dimension(:,:), allocatable :: ugrid_tmp,vgrid_tmp
      complex(r_double), dimension(:), allocatable :: &
      vrtspec_tmp,divspec_tmp
      integer nlon,nlat,ntrunc
      ! compute array dimensions and infer truncation limit
      ! from size of input arrays
      nlon = size(ugrid,1)
      nlat = size(ugrid,2)
      ntrunc = nint((-1.+sqrt(1+8*float(size(vrtspec))))/2.)-1
      if (nlon .ne. current_nlon .or. &
          nlat .ne. current_nlat .or. &
          ntrunc .ne. current_ntrunc) then
          call shtns_init(nlon,nlat,ntrunc,nth,popt)
      end if
      if (default_real .ne. 2) then
         allocate(ugrid_tmp(nlon,nlat))
         allocate(vrtspec_tmp(nlm))
         allocate(vgrid_tmp(nlon,nlat))
         allocate(divspec_tmp(nlm))
         vrtspec_tmp = invlap*rsphere*vrtspec
         divspec_tmp = invlap*rsphere*divspec
         call shtns_sphtor_to_spat(vrtspec_tmp,divspec_tmp,ugrid_tmp,vgrid_tmp)
         ugrid = ugrid_tmp
         vgrid = vgrid_tmp
         deallocate(ugrid_tmp,vgrid_tmp,vrtspec_tmp,divspec_tmp)
      else
         allocate(vrtspec_tmp(nlm))
         allocate(divspec_tmp(nlm))
         vrtspec_tmp = invlap*rsphere*vrtspec
         divspec_tmp = invlap*rsphere*divspec
         call shtns_sphtor_to_spat(vrtspec_tmp,divspec_tmp,ugrid,vgrid)
         deallocate(vrtspec_tmp,divspec_tmp)
      endif
      end subroutine getuv

      subroutine getvrtdivspec(ugrid,vgrid,vrtspec,divspec,rsphere)
! calculate spectral coefficients of vorticity and divergence
! (vrtspec,divspec) given input gridded winds (ugrid,vgrid).
! rsphere is the radius of the sphere in meters.
      real(r_kind), intent(in) :: rsphere
      real(r_kind), dimension(:,:), intent(in) :: ugrid,vgrid
      complex(r_kind), dimension(:), intent(out) :: vrtspec,divspec
      real(r_double), dimension(:,:), allocatable :: ugrid_tmp,vgrid_tmp
      complex(r_double), dimension(:), allocatable :: &
      vrtspec_tmp,divspec_tmp
      integer nlon,nlat,ntrunc
      ! compute array dimensions and infer truncation limit
      ! from size of input arrays.
      nlon = size(ugrid,1)
      nlat = size(ugrid,2)
      ntrunc = nint((-1.+sqrt(1+8*float(size(vrtspec))))/2.)-1
      if (nlon .ne. current_nlon .or. &
          nlat .ne. current_nlat .or. &
          ntrunc .ne. current_ntrunc) then
          call shtns_init(nlon,nlat,ntrunc,nth,popt)
      end if
      if (default_real .ne. 2) then
         allocate(ugrid_tmp(nlon,nlat))
         allocate(vrtspec_tmp(nlm))
         allocate(vgrid_tmp(nlon,nlat))
         allocate(divspec_tmp(nlm))
         ugrid_tmp = ugrid
         vgrid_tmp = vgrid
         call shtns_spat_to_sphtor(ugrid_tmp, vgrid_tmp, vrtspec_tmp, divspec_tmp)
         vrtspec = vrtspec_tmp
         divspec = divspec_tmp
         deallocate(ugrid_tmp,vgrid_tmp,vrtspec_tmp,divspec_tmp)
      else
         call shtns_spat_to_sphtor(ugrid, vgrid, vrtspec, divspec)
      endif
      vrtspec = (lap/rsphere)*vrtspec; divspec = (lap/rsphere)*divspec
      end subroutine getvrtdivspec

      subroutine getgrad(divspec,ugrid,vgrid,rsphere)
! given spectral coefficients compute gridded vector gradient.
! rsphere is the radius of the sphere in meters.
      real(r_kind), intent(in) :: rsphere
      real(r_kind), dimension(:,:), intent(out) :: ugrid,vgrid
      complex(r_kind), dimension(:), intent(in) :: divspec
      complex(r_kind), dimension(:), allocatable :: vrtspec
      real(r_double), dimension(:,:), allocatable :: ugrid_tmp,vgrid_tmp
      complex(r_double), dimension(:), allocatable :: &
      vrtspec_tmp,divspec_tmp
      integer nlon,nlat,ntrunc
      ! compute array dimensions and infer truncation limit
      ! from size of input arrays
      nlon = size(ugrid,1)
      nlat = size(ugrid,2)
      ntrunc = nint((-1.+sqrt(1+8*float(size(divspec))))/2.)-1
      if (nlon .ne. current_nlon .or. &
          nlat .ne. current_nlat .or. &
          ntrunc .ne. current_ntrunc) then
          call shtns_init(nlon,nlat,ntrunc,nth,popt)
      end if
      if (default_real .ne. 2) then
         allocate(ugrid_tmp(nlon,nlat))
         allocate(vrtspec_tmp(nlm))
         allocate(vgrid_tmp(nlon,nlat))
         allocate(divspec_tmp(nlm))
         vrtspec_tmp = 0.
         divspec_tmp = divspec
         call shtns_sphtor_to_spat(vrtspec_tmp,divspec_tmp,ugrid_tmp,vgrid_tmp)
         ugrid = ugrid_tmp
         vgrid = vgrid_tmp
         deallocate(ugrid_tmp,vgrid_tmp,vrtspec_tmp,divspec_tmp)
      else
         allocate(vrtspec(nlm))
         vrtspec = 0.
         call shtns_sphtor_to_spat(vrtspec,divspec,ugrid,vgrid)
         deallocate(vrtspec)
      endif
      ugrid = ugrid/rsphere; vgrid = vgrid/rsphere
      end subroutine getgrad

! that's it!

      end module shtns

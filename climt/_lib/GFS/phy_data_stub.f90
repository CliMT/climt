module phy_data
! stub
! public subroutines:
! init_phydata: allocate and populate arrays.
! destroy_phydata: deallocate arrays.
 use kinds, only: r_kind,r_double
 use params, only: nlats,nlons,nlevs,ntrac,ndimspec

 implicit none
 private
 
 !Temporary variable to store tendencies. This will eventually
 ! be replaced by a tendency function outside the code.

 complex(r_kind), dimension(:,:), allocatable :: spectral_vort_tend, spectral_div_tend, spectral_virtemp_tend

 complex(r_kind), dimension(:,:,:), allocatable :: spectral_tracer_tend
 complex(r_kind), dimension(:), allocatable :: spectral_lnps_tend

 public :: init_phydata, destroy_phydata, wrtout_flx, wrtout_sfc, pwat
 public :: spectral_vort_tend, spectral_div_tend, spectral_virtemp_tend, spectral_tracer_tend, spectral_lnps_tend

 real(r_kind), allocatable, dimension(:,:) :: pwat

 contains

 subroutine init_phydata()
    allocate(pwat(nlons,nlats))
    pwat = 0.
    allocate(spectral_vort_tend(ndimspec,nlevs), spectral_div_tend(ndimspec,nlevs), spectral_virtemp_tend(ndimspec,nlevs))
    allocate(spectral_tracer_tend(ndimspec,nlevs,ntrac))
    allocate(spectral_lnps_tend(ndimspec))
 end subroutine init_phydata

 subroutine destroy_phydata()
    deallocate(pwat)
    deallocate(spectral_lnps_tend)
    deallocate(spectral_vort_tend,spectral_div_tend,spectral_virtemp_tend)
    deallocate(spectral_tracer_tend)
 end subroutine destroy_phydata

 subroutine wrtout_sfc(fhour,filename)
   implicit none
   real(r_kind), intent(in) :: fhour
   character(len=120) filename
 end subroutine wrtout_sfc

 subroutine wrtout_flx(fhour,ta,filename)
   implicit none
   real(r_kind),intent(in) :: fhour
   real(r_double),intent(in) :: ta
   character(len=120),intent(in) :: filename
 end subroutine wrtout_flx

end module phy_data

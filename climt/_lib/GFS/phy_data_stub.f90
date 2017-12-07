module phy_data
 use kinds, only: r_kind,r_double
 use params, only: nlats,nlons,nlevs,ntrac,ndimspec
 use iso_c_binding
 ! stub
 ! public subroutines:
 ! init_phydata: allocate and populate arrays.
 ! destroy_phydata: deallocate arrays.

 implicit none
 private
 
 !Temporary variable to store tendencies. This will eventually
 ! be replaced by a tendency function outside the code.

 complex(r_kind), dimension(:,:), pointer :: spectral_vort_tend, spectral_div_tend, spectral_virtemp_tend

 complex(r_kind), dimension(:,:,:), pointer :: spectral_tracer_tend
 complex(r_kind), dimension(:), pointer :: spectral_lnps_tend

 public :: init_phydata, destroy_phydata, wrtout_flx, wrtout_sfc, pwat
 public :: spectral_vort_tend, spectral_div_tend, spectral_virtemp_tend, spectral_tracer_tend, spectral_lnps_tend

 real(r_kind), pointer, dimension(:,:) :: pwat

 contains

     subroutine init_phydata(py_spec_vort_tend, py_spec_div_tend, py_spec_theta_tend,&
             py_spec_tracer_tend, py_spec_lnps_tend)bind(c, name='gfs_initialise_spectral_physics_arrays')

         complex(c_double), dimension(:,:), target, intent(inout):: py_spec_vort_tend(ndimspec, nlevs),&
            py_spec_div_tend(ndimspec, nlevs), py_spec_theta_tend(ndimspec, nlevs)
         complex(c_double), target, intent(inout):: py_spec_tracer_tend(ndimspec, nlevs, ntrac)
         complex(c_double), target, intent(inout):: py_spec_lnps_tend(ndimspec)

         !print *, 'Initialising physics'
         spectral_vort_tend => py_spec_vort_tend
         spectral_div_tend => py_spec_div_tend
         spectral_virtemp_tend => py_spec_theta_tend
         spectral_lnps_tend => py_spec_lnps_tend
         spectral_tracer_tend => py_spec_tracer_tend

     end subroutine init_phydata

     subroutine init_phydata_pwat(py_pwat)bind(c, name='gfs_initialise_grid_physics_arrays')
         real(c_double), target, intent(inout):: py_pwat(nlons, nlats)

         pwat => py_pwat
         pwat = 0.
     end subroutine

     subroutine destroy_phydata()
         nullify(pwat)
         nullify(spectral_lnps_tend)
         nullify(spectral_vort_tend,spectral_div_tend,spectral_virtemp_tend)
         nullify(spectral_tracer_tend)
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

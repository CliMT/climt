module spectral_data
! spectral model state variable data
! Public subroutines:
! init_specdata: allocate arrays.
! destroy_specdata: deallocate arrays.
 use kinds, only: r_kind
 use params, only: ndimspec, nlevs, ntrac
 use iso_c_binding
 implicit none
 private
 public :: init_specdata,destroy_specdata
 complex(r_kind), public, dimension(:,:), pointer :: &
 vrtspec,divspec,virtempspec
 complex(r_kind), public, dimension(:,:,:), pointer :: tracerspec
! ndimspec by nlevs complex arrays (ndimspec = (ntrunc+1)*(ntrunc)/2):
! these arrays are update each timestep in subroutine advance in run.f90.
! vrtspec: vorticity
! divspec: divergence
! virtempspec: virtural temp
! tracerspec: tracers (ntrac=1 is  specific humidity, ntrac=ntoz ozone, ntrac=ntclw cloud
! condensate)
 complex(r_kind), public, dimension(:), pointer :: lnpsspec,topospec
 real(r_kind), public, dimension(:), pointer :: disspec,diff_prof,dmp_prof
! ndimspec 1-d arrays:
! lnpsspec: log(ps) in Pa
! topospec: orography in meters (static, created in init_dyn)
! disspec: hyperdiffusion operator (static, created in init_dyn)
! diff_prof: vertical profile of hyperdiffusion coefficient.
! dmp_prof: vertical profile of linear momentum drag (active
! at top of model).

 contains
!
! Modifying the specdata init routines to accept arrays from python   
 subroutine init_specdata(pyVrtSpec, pyDivSpec, pyVirtTempSpec, pyTracerSpec, &
         pyTopoSpec, pyLnPsSpec, pyDissSpec, pyDmpProf, pyDiffProf) & 
            bind(c, name='gfs_initialise_spectral_arrays')
    
    complex(c_double_complex), target, intent(inout) :: pyVrtSpec(ndimspec,nlevs), &
    pyDivSpec(ndimspec, nlevs), pyVirtTempSpec(ndimspec, nlevs), &
    pyTracerSpec(ndimspec, nlevs, ntrac), pyTopoSpec(ndimspec), &
    pyLnPsSpec(ndimspec)

    real(c_double), target, intent(inout) :: pyDmpProf(nlevs), &
            pyDiffProf(nlevs), pyDissSpec(ndimspec)

    !print *, "initialising Spectral arrays"
    !print *, shape(pyVrtSpec), shape(pyTracerSpec)

    vrtspec => pyVrtSpec
    divspec => pyDivSpec
    tracerspec => pyTracerSpec
    virtempspec => pyVirtTempSpec

    topospec => pyTopoSpec
    lnpsspec => pyLnPsSpec

    disspec => pyDissSpec
    dmp_prof => pyDmpProf
    diff_prof => pyDiffProf

 end subroutine init_specdata
 subroutine destroy_specdata()
    !print *, 'Deallocating spectral data'
    nullify(vrtspec,divspec,virtempspec)
    nullify(tracerspec)
    nullify(topospec,lnpsspec,disspec,dmp_prof,diff_prof)
 end subroutine destroy_specdata


end module spectral_data

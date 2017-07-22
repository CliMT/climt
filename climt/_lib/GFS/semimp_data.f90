module semimp_data
    ! static data arrays used in semi-implicit time integration scheme.
    ! Public subroutines:
    ! init_semimpdata: allocate and populate arrays.
    ! destroy_semimpdata: deallocate arrays.
    use kinds, only: r_kind, default_real, r_double
    use params, only: dt,ntrunc,nlons,nlats,nlevs
    use pressure_data, only: ak,bk
    use physcons, only: rd => con_rd, cp => con_cp, rerth => con_rerth,&
        kappa => con_rocp
    use iso_c_binding
    implicit none
    private
    public :: init_semimpdata, destroy_semimpdata
    real(r_kind), pointer, public, dimension(:,:,:,:) :: d_hyb_m
    real(r_kind), pointer, public, dimension(:,:) :: amhyb,bmhyb
    real(r_kind), pointer, public, dimension(:) ::  &
        tref,pkref,dpkref,alfaref,svhyb,tor_hyb

    ! scheme derived from Ascher, Spiteri and Ruuth 1997.
    ! Applied Numerical Mathematics DOI:10.1016/S0168-9274(97)00056-1
    ! it's a combination of schemes in sxns 2.5 and 2.6
    ! from that paper (tableaus from 2.5, value of delta from 2.6)
    !  
    ! Double tableau:  c = sqrt(2)/2 = 1/sqrt(2)
    ! *explicit part*  *implicit part* 
    !
    !   0 | 0             0 | 0
    ! 1-c | 1-c 0       1-c | 0 1-c
    !   1 | -c 1+c 0      1 | 0 c 1-c
    ! --------------   ---------------
    !       0  c 1-c          0 c 1-c
    !
    ! when used to discretize 'fast'/'slow' double oscillator eqn
    ! this scheme is:
    ! stable for 'slow' courant numbers < 2*c = sqrt(2) ~ 1.414
    ! stable for all 'fast' courant numbers


    !real(r_kind),parameter :: halfroot2    = sqrt(2.)/2.
    !real(r_kind),parameter, public :: a21  = 1.-halfroot2
    !real(r_kind),parameter, public :: a31  = -halfroot2
    !real(r_kind),parameter, public :: a32  = 1.+halfroot2
    !real(r_kind),parameter, public :: aa21 = 0.
    !real(r_kind),parameter, public :: aa22 = 1.-halfroot2 
    !real(r_kind),parameter, public :: aa31 = 0.
    !real(r_kind),parameter, public :: aa32 = halfroot2
    !real(r_kind),parameter, public :: aa33 = 1.-halfroot2 
    !real(r_kind),parameter, public :: b1   = 0.  
    !real(r_kind),parameter, public :: b2   = halfroot2
    !real(r_kind),parameter, public :: b3   = 1.-halfroot2
    !real(r_kind),parameter, public :: bb1  = 0. 
    !real(r_kind),parameter, public :: bb2  = halfroot2
    !real(r_kind),parameter, public :: bb3  = 1.-halfroot2
    !real(r_kind),parameter, public :: bb4  = 0.

    ! Lorenz 3-cycle + trapezoidal
    !
    ! has very weak instability for purely imaginary double
    ! oscillator.
    ! Largest amp factor for 'slow' courant # < sqrt(3) ~ 1.73,
    ! which is not a problem for models with weak damping.
    !
    ! Double tableau:
    !
    !  *explicit part*                *implicit part* 
    !
    !   0  | 0                       0  | 0
    !  1/3 | 1/3  0                 1/3 | 1/6 1/6 0
    !  2/3 | 1/6  1/2 0             2/3 | 1/3  0 1/3  0
    !   1  | 1/2 -1/2 1  0           1  | 3/8  0 3/8 1/4 
    ! ---------------------         ----------------------
    !        1/2 -1/2 1  0                3/8  0 3/8 1/4

    ! real(r_kind),parameter, public :: a21  = 1./3.
    ! real(r_kind),parameter, public :: a31  = 1./6.
    ! real(r_kind),parameter, public :: a32  = 1./2.
    ! real(r_kind),parameter, public :: aa21 = 1./6.
    ! real(r_kind),parameter, public :: aa22 = 1./6.
    ! real(r_kind),parameter, public :: aa31 = 1./3.
    ! real(r_kind),parameter, public :: aa32 = 0.
    ! real(r_kind),parameter, public :: aa33 = 1./3.
    ! real(r_kind),parameter, public :: b1   = 1./2.
    ! real(r_kind),parameter, public :: b2   = -1./2.
    ! real(r_kind),parameter, public :: b3   = 1.0 
    ! real(r_kind),parameter, public :: bb1  = 3./8.
    ! real(r_kind),parameter, public :: bb2  = 0.
    ! real(r_kind),parameter, public :: bb3  = 3./8.  
    ! real(r_kind),parameter, public :: bb4  = 1./4.

    ! ssp/tvd rk3 (shu 1988, http://epubs.siam.org/doi/abs/10.1137/0909073) IMEX
    ! minimizes spurious oscillations - good for advection of tracers
    ! (no negative entries in Butcher tableau)
    ! stability similar to Lorenz 3-cycle IMEX
    !
    ! Double tableau:
    !
    !  *explicit part*                *implicit part* 
    !
    !   0  | 0                       0  | 0
    !   1  | 1   0                   1  | 0.635  0.365 
    !  1/2 | 1/4 1/4 0              1/2 | 0.3175 0.0   0.1825  
    !   1  | 1/6 1/6 2/3 0           1  | 0.35   0.0   0.3    0.35
    ! ---------------------         -------------------------------
    !        1/6 1/6 2/3 0                0.35   0.0   0.3    0.35
    !
    real(r_kind),parameter, public :: a21  = 1.0
    real(r_kind),parameter, public :: a31  = 0.25  
    real(r_kind),parameter, public :: a32  = 0.25  
    real(r_kind),parameter, public :: aa21 = 0.635
    real(r_kind),parameter, public :: aa22 = 0.365
    real(r_kind),parameter, public :: aa31 = 0.3175 
    real(r_kind),parameter, public :: aa32 = 0.
    real(r_kind),parameter, public :: aa33 = 0.1825  
    real(r_kind),parameter, public :: b1   = 1./6.
    real(r_kind),parameter, public :: b2   = 1./6.
    real(r_kind),parameter, public :: b3   = 2./3.
    real(r_kind),parameter, public :: bb1  = 0.35  
    real(r_kind),parameter, public :: bb2  = 0.0
    real(r_kind),parameter, public :: bb3  = 0.3    
    real(r_kind),parameter, public :: bb4  = 0.35  

    real(r_kind),parameter, public :: ref_temp = 300.
    real(r_kind),parameter, public :: ref_press = 800.e2

contains

    subroutine init_semimp_arrays(py_d_hyb_m, py_amhyb, py_bmhyb,&
        py_tref, py_pkref, py_dpkref, py_alfaref, py_svhyb, py_tor_hyb)bind(c, name='gfs_initialise_timestepper_arrays')

    real(c_double), target :: py_d_hyb_m(nlevs, nlevs, ntrunc+1, 3)
    real(c_double), target:: py_amhyb(nlevs, nlevs), py_bmhyb(nlevs, nlevs)
    real(c_double), target:: py_tref(nlevs), py_pkref(nlevs+1), py_dpkref(nlevs),&
        py_alfaref(nlevs), py_svhyb(nlevs), py_tor_hyb(nlevs)

    d_hyb_m => py_d_hyb_m
    amhyb => py_amhyb
    bmhyb => py_bmhyb
    tref => py_tref
    pkref => py_pkref
    dpkref => py_dpkref
    alfaref => py_alfaref
    svhyb => py_svhyb
    tor_hyb => py_tor_hyb

    end subroutine

    subroutine init_semimpdata()
        real(r_double) rnn1
        integer irow,icol,icolbeg,i,j,k,n,icolend,nn,iret
        real(r_double), allocatable, dimension(:,:) :: yecm,tecm,ym,rim
        real(r_double), allocatable, dimension(:) :: vecm
        integer, allocatable, dimension(:) :: ipiv
        real(r_kind) epstiny,const

        ! make sure progression in time the same in explicit and implicit
        ! parts of additive runge-kutta scheme.
        epstiny = 1.e-5
        if (abs(a21-(aa21+aa22)) > epstiny .or. abs(a31+a32-(aa31+aa32+aa33)) > epstiny) then
            print *,a21,a31+a32,abs(a21-(aa21+aa22)) 
            print *,aa22,aa32+aa33,abs(a31+a32-(aa31+aa32+aa33))
            print *,'time progression different in explicit and implicit parts of RK'
            print *,'fix in semimp_data.f90!'
            stop
        endif
        ! check to make sure b1+b2+b3=bb1+bb2+bb3+bb4=1
        if (abs(b1+b2+b3-1.) > epstiny .or. abs(bb1+bb2+bb3+bb4-1) > epstiny) then
            print *,'b1+b2+b3 or bb1+bb2+bb3+bb4 != 1'
            print *,'fix in semimp_data.f90!'
            stop
        endif
        ! check to make sure scheme is 2nd order
        if (abs(a21*b2+(a31+a32)*b3-0.5) > epstiny .or.&
            abs((aa21+aa22)*bb2+(aa31+aa32+aa33)*bb3+bb4-0.5) > epstiny) then
            print *,a21*b2 + (a31+a32)*b3,abs(a21*b2+(a31+a32)*b3+bb4-0.5)
            print *,(aa21+aa22)*bb2+(aa31+aa32+aa33)*bb3 ,&
                abs((aa21+aa22)*bb2+(aa31+aa32+aa33)*bb3+bb4-0.5)
            print *,'scheme not 2nd order'
            print *,'fix in semimp_data.f90!'
            stop
        endif

        ! module vars
        allocate(ipiv(nlevs))
        !allocate(amhyb(nlevs,nlevs),bmhyb(nlevs,nlevs))
        !allocate(tref(nlevs),pkref(nlevs+1),dpkref(nlevs),alfaref(nlevs))
        !allocate(svhyb(nlevs),tor_hyb(nlevs))
        !allocate(d_hyb_m(nlevs,nlevs,ntrunc+1,3))
        ! temp storage.
        allocate(yecm(nlevs,nlevs),tecm(nlevs,nlevs))
        allocate(ym(nlevs,nlevs))
        allocate(vecm(nlevs))
        allocate(rim(nlevs,nlevs)) ! identity matrix

        tref = ref_temp
        do k=1,nlevs+1
        pkref(k) = ak(k) + bk(k)*ref_press
        enddo
        do k=1,nlevs
        dpkref(k) = pkref(k+1)-pkref(k)
        enddo
        alfaref(1) = log(2.)
        do k=2,nlevs
        alfaref(k) = 1.-(pkref(k)/dpkref(k))*log(pkref(k+1)/pkref(k))
        enddo
        yecm=0.
        do irow=1,nlevs
        yecm(irow,irow)=alfaref(irow)*rd
        icolbeg=irow+1
        if(icolbeg<=nlevs)then
            do icol=icolbeg,nlevs
            yecm(irow,icol)=rd*log( pkref(icol+1)/pkref(icol) )
            enddo
        endif
        enddo
        tecm=0.
        do irow=1,nlevs
        tecm(irow,irow)=kappa*tref(irow)*alfaref(irow)
        icolend=irow-1
        do icol=1,icolend
        tecm(irow,icol)=(kappa*tref(irow)*dpkref(icol)/ &
            dpkref(irow))*log(pkref(irow+1)/pkref(irow))
        enddo
        enddo
        vecm=dpkref/ref_press
        ! amhyb is operator for linearized pres grad term (geopot part)
        ! (operatores on virt temp)
        ! bmhyb is operator for linearized energy conv term
        ! (operates on div)
        do j=1,nlevs
        svhyb(j)=vecm(nlevs+1-j) ! for linearized lnps tend term (div)
        do k=1,nlevs
        amhyb(k,j)=yecm(nlevs+1-k,nlevs+1-j)
        bmhyb(k,j)=tecm(nlevs+1-k,nlevs+1-j)
        enddo
        enddo
        amhyb=amhyb/rerth**2 
        ! times lnps in linearized pgf term (lapacian in div eqn)
        tor_hyb=rd*tref/rerth**2

        rim = 0.
        do k=1,nlevs
        rim(k,k) = 1.
        enddo
        ! computations that do not depend on wavenumber
        do i=1,nlevs
        do j=1,nlevs
        ym(i,j) = tor_hyb(i)*svhyb(j)
        enddo
        do k=1,nlevs
        do j=1,nlevs
        ym(i,j) = ym(i,j) + amhyb(i,k)*bmhyb(k,j)
        enddo
        enddo
        enddo

        ! computations that do depend on wavenumber
        ! enabling openmp for this loop doesn't work with intel MKL
        do k=1,3
        if (k .eq. 1)  const = aa22
        if (k .eq. 2)  const = aa33
        if (k .eq. 3)  const = bb4
        !!$omp parallel do private(nn,n,rnn1,yecm,ipiv,iret,vecm)
        do nn=1,ntrunc+1
        n = nn-1
        rnn1 = n*(n+1)
        yecm = rim + (const*dt)**2*rnn1*ym
        ! invert matrix using LAPACK, save in d_hyb_m
        call dgetrf(nlevs,nlevs,yecm,nlevs,ipiv,iret)
        call dgetri(nlevs,yecm,nlevs,ipiv,vecm,nlevs,iret)
        d_hyb_m(:,:,nn,k) = yecm
        enddo
        !!$omp end parallel do 
        enddo
        deallocate(rim,yecm,tecm,ym,vecm,ipiv)

    end subroutine init_semimpdata

    subroutine destroy_semimpdata()
        nullify(amhyb,bmhyb,d_hyb_m)
        nullify(tref,pkref,dpkref,alfaref,svhyb,tor_hyb)
    end subroutine destroy_semimpdata

end module semimp_data

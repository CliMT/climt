 module phy_run
! compute physics tendencies for simplifed held-suarez forcing
! getphytend: compute physics tendencies.
! Public subroutines:
! getphytend: compute newtonian damping of temperature, drag terms
! in vort and div eqns.

 use params, only: nlevs,nlons,nlats,ndimspec,ntrac,pdryini
 use kinds, only: r_kind, r_double
 use shtns, only: grdtospec, lats, areawts, print_bool
 !use grid_data, only: vrtg,divg,virtempg
 !use pressure_data, only:  prs,psg
 !use spectral_data, only: lnpsspec

 use phy_data, only: spectral_vort_tend, spectral_div_tend, spectral_virtemp_tend, spectral_tracer_tend, spectral_lnps_tend
 use,intrinsic  :: iso_c_binding


 implicit none
 private


 ! Define interface of call-back routine.
   !    ABSTRACT INTERFACE
   !      SUBROUTINE callback (x)
 !         USE, INTRINSIC :: ISO_C_BINDING
 !          REAL(KIND=C_DOUBLE), INTENT(IN), VALUE :: x
!         END SUBROUTINE callback
!       END INTERFACE
 ! Define interface of call-back routine.
       ABSTRACT INTERFACE
         SUBROUTINE callback ()
           USE, INTRINSIC :: ISO_C_BINDING
         END SUBROUTINE callback
       END INTERFACE

 public :: getphytend, set_tendencies
 !public :: getphytend, gfsRegisterPhysicsCallback, set_tendencies

 !procedure(callback), pointer :: pyPhysicsCallback => null()

 integer(c_int) :: phy_trac
 contains

 !subroutine gfsRegisterPhysicsCallback(callback) bind(c,name='gfs_register_physics_callback')

 !   type(c_funptr), intent(in), value :: callback

 !   print *, 'Registering callback'
 !   print *, callback
 !   call c_f_procpointer(callback, pyPhysicsCallback)

 !   call pyPhysicsCallback()

 !end subroutine gfsRegisterPhysicsCallback

 subroutine getphytend(dvrtspecdt,ddivspecdt,dvirtempspecdt,dtracerspecdt,dlnpsspecdt,t,dt)
     complex(r_kind), intent(inout), dimension(ndimspec,nlevs) :: &
         dvrtspecdt,ddivspecdt,dvirtempspecdt
     complex(r_kind), intent(out), dimension(ndimspec,nlevs,ntrac) :: &
         dtracerspecdt
     real(r_double), intent(in) :: t,dt
     complex(r_kind), intent(inout), dimension(ndimspec) :: dlnpsspecdt
     integer k, trac


     phy_trac = ntrac
     !$omp parallel do private(k)
     do k=1,nlevs
     dvrtspecdt(:,k) = spectral_vort_tend(:,k)
     ddivspecdt(:,k) = spectral_div_tend(:,k)
     dvirtempspecdt(:,k) = spectral_virtemp_tend(:,k)
     do trac=1,phy_trac
         dtracerspecdt(:,k,trac) = spectral_tracer_tend(:,k,trac)
     enddo
     enddo
     !$omp end parallel do
     dlnpsspecdt(:) = spectral_lnps_tend(:)


 end subroutine getphytend


subroutine set_tendencies(dvrtgdt,ddivgdt,dvirtempgdt,dlnpsgdt,dtracergdt) bind(c, name='gfs_set_tendencies')

    real(c_double), intent(in), dimension(nlons,nlats,nlevs) :: &
        dvrtgdt,ddivgdt,dvirtempgdt
    real(c_double), intent(in), dimension(nlons,nlats,nlevs,ntrac) :: &
        dtracergdt
    real(c_double), intent(in), dimension(nlons,nlats) :: dlnpsgdt
    !complex(r_kind), dimension(:,:),allocatable :: forcingspec
    real(c_double), dimension(nlons, nlats) :: temp_grid_tracer
    integer k,trac


    !print *, dtracergdt(:, :, 1, 1)
    phy_trac = ntrac
    !print *, maxval(abs(dvrtgdt(:,:,9))), maxval(abs(ddivgdt)), maxval(abs(dvirtempgdt(:,:,9)))
    !print *, maxval(abs(dtracergdt(:,:,:,1)))
    !$omp parallel do private(k)
    do k=1,nlevs
    call grdtospec(dvrtgdt(:,:,k), spectral_vort_tend(:,k))
    call grdtospec(ddivgdt(:,:,k), spectral_div_tend(:,k))
    call grdtospec(dvirtempgdt(:,:,k), spectral_virtemp_tend(:,k))
    if (k == 1) then
        print_bool = 1
    end if
    print_bool = 0
    !spectral_vort_tend(:,k) = dvrtspecdt(:,k)
    !spectral_div_tend(:,k) = ddivspecdt(:,k)
    !spectral_virtemp_tend(:,k) = dvirtempspecdt(:,k)
    enddo
    !$omp end parallel do
    do trac=1,phy_trac
        do k=1,nlevs
        !print *, ''
        !print *, 'Level: ', k
        temp_grid_tracer(:,:) = dtracergdt(:,:,k,trac)
        !print *, temp_grid_tracer
        !print *,
            call grdtospec(temp_grid_tracer, spectral_tracer_tend(:,k,trac))
        enddo
    enddo

    call grdtospec(dlnpsgdt(:,:), spectral_lnps_tend(:))
    !spectral_lnps_tend(:) = dlnpsspecdt(:)
    !print *, maxval(abs(realpart((spectral_virtemp_tend(:,:)))))
    !print *, maxval(abs(realpart((spectral_tracer_tend(:,:,:)))))
    !print *, maxval(abs(imagpart((spectral_tracer_tend(:,:,1)))))

    !deallocate(spectral_lnps_tend)
    !deallocate(spectral_vort_tend,spectral_div_tend,spectral_virtemp_tend)
    !deallocate(forcingspec)
end subroutine set_tendencies


end module phy_run

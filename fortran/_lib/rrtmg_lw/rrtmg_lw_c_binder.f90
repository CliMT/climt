! This file provides C bindings to the selected RRTMG LW routines.

module rrtmg_lw_c_binder

        use iso_c_binding
        use parkind, only : im => kind_im, rb => kind_rb
        implicit none
        contains

          subroutine rrtmg_lw_set_constants &
              (pi_, grav_, planck_, boltz_, clight_, &
               avogad_, alosmt_, gascon_, sbcnst_, secdy_)&
             bind(c)

            use rrlw_con

            real(kind=rb), intent(in) :: pi_, grav_, planck_
            real(kind=rb), intent(in) :: boltz_, clight_

            real(kind=rb), intent(in) :: avogad_, alosmt_, gascon_
            real(kind=rb), intent(in) :: sbcnst_, secdy_

            pi = pi_
            grav = grav_
            planck = planck_
            boltz = boltz_
            clight = clight_
            avogad = avogad_
            alosmt = alosmt_
            gascon = gascon_
            sbcnst = sbcnst_
            secdy = secdy_

            radcn1 = 2. * planck * clight * clight * 1.e-07
            radcn2 = planck * clight / boltz

          end subroutine rrtmg_lw_set_constants

          subroutine rrtmg_lw_ini_wrapper(cpdair)&
             bind(c)

            use rrtmg_lw_init, only: rrtmg_lw_ini

            real(kind=rb), intent(in) :: cpdair

            call rrtmg_lw_ini(cpdair)

          end subroutine rrtmg_lw_ini_wrapper

          subroutine mcica_subcol_lw_wrapper &
              (iplon, ncol, nlay, icld, permuteseed, irng, play, &
               cldfrac, ciwp, clwp, rei, rel, &
               tauc, &
               cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
               taucmcl)&
             bind(c)

            use mcica_subcol_gen_lw, only: mcica_subcol_lw
            use parrrtm, only : nbndlw, ngptlw

            integer(kind=im), intent(in) :: iplon
            integer(kind=im), intent(in) :: ncol
            integer(kind=im), intent(in) :: nlay
            integer(kind=im), intent(in) :: icld
            integer(kind=im), intent(in) :: permuteseed
            integer(kind=im), intent(inout) :: irng
            real(kind=rb), intent(in) :: play(ncol, nlay)

            real(kind=rb), intent(in) :: cldfrac(ncol, nlay)
            real(kind=rb), intent(in) :: ciwp(ncol, nlay)
            real(kind=rb), intent(in) :: clwp(ncol, nlay)
            real(kind=rb), intent(in) :: rei(ncol, nlay)
            real(kind=rb), intent(in) :: rel(ncol, nlay)

            real(kind=rb), intent(in) :: tauc(nbndlw, ncol, nlay)

            real(kind=rb), intent(out) :: cldfmcl(ngptlw, ncol, nlay)
            real(kind=rb), intent(out) :: ciwpmcl(ngptlw, ncol, nlay)
            real(kind=rb), intent(out) :: clwpmcl(ngptlw, ncol, nlay)
            real(kind=rb), intent(out) :: reicmcl(ncol, nlay)
            real(kind=rb), intent(out) :: relqmcl(ncol, nlay)

            real(kind=rb), intent(out) :: taucmcl(ngptlw, ncol, nlay)

            call mcica_subcol_lw &
                (iplon, ncol, nlay, icld, permuteseed, irng, play, &
                 cldfrac, ciwp, clwp, rei, rel, &
                 tauc, &
                 cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                 taucmcl)

          end subroutine mcica_subcol_lw_wrapper

          subroutine rrtmg_lw_mcica_wrapper &
              (ncol, nlay, icld, idrv, &
               play, plev, tlay, tlev, tsfc, &
               h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
               cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
               inflglw, iceflglw, liqflglw, cldfmcl, &
               taucmcl, &
               ciwpmcl, clwpmcl, reicmcl, relqmcl, &
               tauaer, &
               uflx, dflx, hr, uflxc, dflxc, hrc, &
! optional I/O
               duflx_dt, duflxc_dt)&
             bind(c)

            use rrtmg_lw_rad, only : rrtmg_lw
            use parrrtm, only : ngptlw, nbndlw

            integer(kind=im), intent(in) :: ncol
            integer(kind=im), intent(in) :: nlay
            integer(kind=im), intent(inout) :: icld
            integer(kind=im), intent(in) :: idrv

            real(kind=rb), intent(in) :: play(ncol, nlay)
            real(kind=rb), intent(in) :: plev(ncol, nlay + 1)
            real(kind=rb), intent(in) :: tlay(ncol, nlay)
            real(kind=rb), intent(in) :: tlev(ncol, nlay + 1)
            real(kind=rb), intent(in) :: tsfc(ncol)

            real(kind=rb), intent(in) :: h2ovmr(ncol, nlay)
            real(kind=rb), intent(in) :: o3vmr(ncol, nlay)
            real(kind=rb), intent(in) :: co2vmr(ncol, nlay)
            real(kind=rb), intent(in) :: ch4vmr(ncol, nlay)
            real(kind=rb), intent(in) :: n2ovmr(ncol, nlay)
            real(kind=rb), intent(in) :: o2vmr(ncol, nlay)

            real(kind=rb), intent(in) :: cfc11vmr(ncol, nlay)
            real(kind=rb), intent(in) :: cfc12vmr(ncol, nlay)
            real(kind=rb), intent(in) :: cfc22vmr(ncol, nlay)
            real(kind=rb), intent(in) :: ccl4vmr(ncol, nlay)
            real(kind=rb), intent(in) :: emis(ncol, nlay)

            integer(kind=im), intent(in) :: inflglw
            integer(kind=im), intent(in) :: iceflglw
            integer(kind=im), intent(in) :: liqflglw
            real(kind=rb), intent(in) :: cldfmcl(ngptlw, ncol, nlay)

            real(kind=rb), intent(in) :: taucmcl(ngptlw, ncol, nlay)

            real(kind=rb), intent(in) :: ciwpmcl(ngptlw, ncol, nlay)
            real(kind=rb), intent(in) :: clwpmcl(ngptlw, ncol, nlay)
            real(kind=rb), intent(in) :: reicmcl(ncol, nlay)
            real(kind=rb), intent(in) :: relqmcl(ncol, nlay)

            real(kind=rb), intent(in) :: tauaer(ncol, nlay, nbndlw)

            real(kind=rb), intent(out) :: uflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: dflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: hr(ncol, nlay)
            real(kind=rb), intent(out) :: uflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: dflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: hrc(ncol, nlay)

            real(kind=rb), intent(out), optional :: duflx_dt(ncol, nlay + 1)
            real(kind=rb), intent(out), optional :: duflxc_dt(ncol, nlay + 1)


            call rrtmg_lw &
                (ncol, nlay, icld, idrv, &
                 play, plev, tlay, tlev, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                 cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
                 inflglw, iceflglw, liqflglw, cldfmcl, &
                 taucmcl, &
                 ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                 tauaer, &
                 uflx, dflx, hr, uflxc, dflxc, hrc, &
! optional I/O
                 duflx_dt, duflxc_dt)


          end subroutine rrtmg_lw_mcica_wrapper

          subroutine rrtmg_lw_nomcica_wrapper &
              (ncol, nlay, icld, idrv, &
               play, plev, tlay, tlev, tsfc, &
               h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
               cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
               inflglw, iceflglw, liqflglw, cldfr, &
               taucld, &
               cicewp, cliqwp, reice, reliq, &
               tauaer, &
               uflx, dflx, hr, uflxc, dflxc, hrc, &
! optional I/O
               duflx_dt, duflxc_dt)&
             bind(c)

            use rrtmg_lw_rad_nomcica, only : rrtmg_lw
            use parrrtm, only : nbndlw

            integer(kind=im), intent(in) :: ncol
            integer(kind=im), intent(in) :: nlay
            integer(kind=im), intent(inout) :: icld
            integer(kind=im), intent(in) :: idrv

            real(kind=rb), intent(in) :: play(ncol, nlay)
            real(kind=rb), intent(in) :: plev(ncol, nlay + 1)
            real(kind=rb), intent(in) :: tlay(ncol, nlay)
            real(kind=rb), intent(in) :: tlev(ncol, nlay + 1)
            real(kind=rb), intent(in) :: tsfc(ncol)

            real(kind=rb), intent(in) :: h2ovmr(ncol, nlay)
            real(kind=rb), intent(in) :: o3vmr(ncol, nlay)
            real(kind=rb), intent(in) :: co2vmr(ncol, nlay)
            real(kind=rb), intent(in) :: ch4vmr(ncol, nlay)
            real(kind=rb), intent(in) :: n2ovmr(ncol, nlay)
            real(kind=rb), intent(in) :: o2vmr(ncol, nlay)

            real(kind=rb), intent(in) :: cfc11vmr(ncol, nlay)
            real(kind=rb), intent(in) :: cfc12vmr(ncol, nlay)
            real(kind=rb), intent(in) :: cfc22vmr(ncol, nlay)
            real(kind=rb), intent(in) :: ccl4vmr(ncol, nlay)
            real(kind=rb), intent(in) :: emis(ncol, nlay)

            integer(kind=im), intent(in) :: inflglw
            integer(kind=im), intent(in) :: iceflglw
            integer(kind=im), intent(in) :: liqflglw
            real(kind=rb), intent(in) :: cldfr(ncol, nlay)

            real(kind=rb), intent(in) :: taucld(nbndlw, ncol, nlay)

            real(kind=rb), intent(in) :: cicewp(ncol, nlay)
            real(kind=rb), intent(in) :: cliqwp(ncol, nlay)
            real(kind=rb), intent(in) :: reice(ncol, nlay)
            real(kind=rb), intent(in) :: reliq(ncol, nlay)

            real(kind=rb), intent(in) :: tauaer(ncol, nlay, nbndlw)

            real(kind=rb), intent(out) :: uflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: dflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: hr(ncol, nlay)
            real(kind=rb), intent(out) :: uflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: dflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: hrc(ncol, nlay)

            real(kind=rb), intent(out), optional :: duflx_dt(ncol, nlay + 1)
            real(kind=rb), intent(out), optional :: duflxc_dt(ncol, nlay + 1)


            call rrtmg_lw &
                (ncol, nlay, icld, idrv, &
                 play, plev, tlay, tlev, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                 cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
                 inflglw, iceflglw, liqflglw, cldfr, &
                 taucld, &
                 cicewp, cliqwp, reice, reliq, &
                 tauaer, &
                 uflx, dflx, hr, uflxc, dflxc, hrc, &
! optional I/O
                 duflx_dt, duflxc_dt)


          end subroutine rrtmg_lw_nomcica_wrapper


    end module rrtmg_lw_c_binder

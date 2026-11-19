! This file provides C bindings to the selected RRTMG SW routines.
!
! The main idea to have this separate file is to be able to integrate
! the original RRTM files with minimal changes.
!
! The wrapper functions are required to avoid dummy-dimensioned
! (e.g. arr(:,:)) arrays in the C interface. Arrays declared that way
! are common for the original code but cannot be bound with the
! Fortran 2003 ISO_C_BINDING module. Therefore we declare all the
! arrays with explicit dimensions and pass them to the original routines.

      module rrtmg_sw_c_binder

        use iso_c_binding
        use parkind, only : im => kind_im, rb => kind_rb
        implicit none
        contains

          subroutine rrtmg_sw_set_constants &
              (pi_, grav_, planck_, boltz_, clight_, &
               avogad_, alosmt_, gascon_, sbcnst_, secdy_)&
             bind(c)

            use rrsw_con

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

          end subroutine rrtmg_sw_set_constants

          subroutine rrtmg_sw_ini_wrapper(cpdair)&
             bind(c)

            use rrtmg_sw_init, only: rrtmg_sw_ini

            real(kind=rb), intent(in) :: cpdair

            call rrtmg_sw_ini(cpdair)

          end subroutine rrtmg_sw_ini_wrapper

          subroutine mcica_subcol_sw_wrapper &
              (iplon, ncol, nlay, icld, permuteseed, irng, play, &
               cldfrac, ciwp, clwp, rei, rel, &
               tauc, ssac, asmc, fsfc, &
               cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
               taucmcl, ssacmcl, asmcmcl, fsfcmcl)&
             bind(c)

            use mcica_subcol_gen_sw, only: mcica_subcol_sw
            use parrrsw, only : nbndsw, ngptsw

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

            real(kind=rb), intent(in) :: tauc(nbndsw, ncol, nlay)
            real(kind=rb), intent(in) :: ssac(nbndsw, ncol, nlay)
            real(kind=rb), intent(in) :: asmc(nbndsw, ncol, nlay)
            real(kind=rb), intent(in) :: fsfc(nbndsw, ncol, nlay)

            real(kind=rb), intent(out) :: cldfmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(out) :: ciwpmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(out) :: clwpmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(out) :: reicmcl(ncol, nlay)
            real(kind=rb), intent(out) :: relqmcl(ncol, nlay)

            real(kind=rb), intent(out) :: taucmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(out) :: ssacmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(out) :: asmcmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(out) :: fsfcmcl(ngptsw, ncol, nlay)

            call mcica_subcol_sw &
                (iplon, ncol, nlay, icld, permuteseed, irng, play, &
                 cldfrac, ciwp, clwp, rei, rel, &
                 tauc, ssac, asmc, fsfc, &
                 cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                 taucmcl, ssacmcl, asmcmcl, fsfcmcl)

          end subroutine mcica_subcol_sw_wrapper

          subroutine rrtmg_sw_mcica_wrapper &
              (ncol, nlay, icld, iaer, &
               play, plev, tlay, tlev, tsfc, &
               h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
               asdir, asdif, aldir, aldif, &
               coszen, adjes, dyofyr, scon, isolvar, &
               inflgsw, iceflgsw, liqflgsw, cldfmcl, &
               taucmcl, ssacmcl, asmcmcl, fsfcmcl, &
               ciwpmcl, clwpmcl, reicmcl, relqmcl, &
               tauaer, ssaaer, asmaer, ecaer, &
               swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
! optional I/O
               bndsolvar, indsolvar, solcycfrac)&
             bind(c)

            use rrtmg_sw_rad, only : rrtmg_sw
            use parrrsw, only : ngptsw, nbndsw, naerec

            integer(kind=im), intent(in) :: ncol
            integer(kind=im), intent(in) :: nlay
            integer(kind=im), intent(inout) :: icld
            integer(kind=im), intent(inout) :: iaer

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

            real(kind=rb), intent(in) :: asdir(ncol)
            real(kind=rb), intent(in) :: asdif(ncol)
            real(kind=rb), intent(in) :: aldir(ncol)
            real(kind=rb), intent(in) :: aldif(ncol)

            real(kind=rb), intent(in) :: coszen(ncol)
            real(kind=rb), intent(in) :: adjes
            integer(kind=im), intent(in) :: dyofyr
            real(kind=rb), intent(in) :: scon
            integer(kind=im), intent(in) :: isolvar

            integer(kind=im), intent(in) :: inflgsw
            integer(kind=im), intent(in) :: iceflgsw
            integer(kind=im), intent(in) :: liqflgsw
            real(kind=rb), intent(in) :: cldfmcl(ngptsw, ncol, nlay)

            real(kind=rb), intent(in) :: taucmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(in) :: ssacmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(in) :: asmcmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(in) :: fsfcmcl(ngptsw, ncol, nlay)

            real(kind=rb), intent(in) :: ciwpmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(in) :: clwpmcl(ngptsw, ncol, nlay)
            real(kind=rb), intent(in) :: reicmcl(ncol, nlay)
            real(kind=rb), intent(in) :: relqmcl(ncol, nlay)

            real(kind=rb), intent(in) :: tauaer(ncol, nlay, nbndsw)
            real(kind=rb), intent(in) :: ssaaer(ncol, nlay, nbndsw)
            real(kind=rb), intent(in) :: asmaer(ncol, nlay, nbndsw)
            real(kind=rb), intent(in) :: ecaer(ncol, nlay, naerec)

            real(kind=rb), intent(out) :: swuflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swdflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swhr(ncol, nlay)
            real(kind=rb), intent(out) :: swuflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swdflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swhrc(ncol, nlay)

            real(kind=rb), intent(in), optional :: bndsolvar(nbndsw)
            real(kind=rb), intent(inout), optional :: indsolvar(2)
            real(kind=rb), intent(in), optional :: solcycfrac

            call rrtmg_sw &
                (ncol, nlay, icld, iaer, &
                 play, plev, tlay, tlev, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                 asdir, asdif, aldir, aldif, &
                 coszen, adjes, dyofyr, scon, isolvar, &
                 inflgsw, iceflgsw, liqflgsw, cldfmcl, &
                 taucmcl, ssacmcl, asmcmcl, fsfcmcl, &
                 ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                 tauaer, ssaaer, asmaer, ecaer, &
                 swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
! optional I/O
                 bndsolvar, indsolvar, solcycfrac)

          end subroutine rrtmg_sw_mcica_wrapper

          subroutine rrtmg_sw_nomcica_wrapper &
              (ncol, nlay, icld, iaer, &
               play, plev, tlay, tlev, tsfc, &
               h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
               asdir, asdif, aldir, aldif, &
               coszen, adjes, dyofyr, scon, isolvar, &
               inflgsw, iceflgsw, liqflgsw, cldfr, &
               taucld, ssacld, asmcld, fsfcld, &
               cicewp, cliqwp, reice, reliq, &
               tauaer, ssaaer, asmaer, ecaer, &
               swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
! optional I/O
               bndsolvar, indsolvar, solcycfrac)&
             bind(c)

            use rrtmg_sw_rad_nomcica, only : rrtmg_sw
            use parrrsw, only : nbndsw, naerec

            integer(kind=im), intent(in) :: ncol
            integer(kind=im), intent(in) :: nlay
            integer(kind=im), intent(inout) :: icld
            integer(kind=im), intent(inout) :: iaer

            real(kind=rb), intent(in) :: play(ncol, nlay)
            real(kind=rb), intent(in) :: plev(ncol, nlay + 1)
            real(kind=rb), intent(in) :: tlay(ncol, nlay)
            real(kind=rb), intent(in) :: tlev(ncol, nlay + 1)
            real(kind=rb), intent(in) :: tsfc(ncol)

            real(kind=rb), intent(in) :: h2ovmr(ncol, nlay)
            real(kind=rb), intent(in) :: o3vmr(ncol, nlay)
            real(kind=rb), intent(in) :: co2vmr(ncol, nlay)
            real(kind=rb), intent(in) :: ch4vmr(ncol, nlay)
            real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)
            real(kind=rb), intent(in) :: o2vmr(ncol, nlay)

            real(kind=rb), intent(in) :: asdir(ncol)
            real(kind=rb), intent(in) :: asdif(ncol)
            real(kind=rb), intent(in) :: aldir(ncol)
            real(kind=rb), intent(in) :: aldif(ncol)

            real(kind=rb), intent(in) :: coszen(ncol)
            real(kind=rb), intent(in) :: adjes
            integer(kind=im), intent(in) :: dyofyr
            real(kind=rb), intent(in) :: scon
            integer(kind=im), intent(in) :: isolvar

            integer(kind=im), intent(in) :: inflgsw
            integer(kind=im), intent(in) :: iceflgsw
            integer(kind=im), intent(in) :: liqflgsw
            real(kind=rb), intent(in) :: cldfr(ncol, nlay)

            real(kind=rb), intent(in) :: taucld(nbndsw, ncol, nlay)
            real(kind=rb), intent(in) :: ssacld(nbndsw, ncol, nlay)
            real(kind=rb), intent(in) :: asmcld(nbndsw, ncol, nlay)
            real(kind=rb), intent(in) :: fsfcld(nbndsw, ncol, nlay)

            real(kind=rb), intent(in) :: cicewp(ncol, nlay)
            real(kind=rb), intent(in) :: cliqwp(ncol, nlay)
            real(kind=rb), intent(in) :: reice(ncol, nlay)
            real(kind=rb), intent(in) :: reliq(ncol,nlay)

            real(kind=rb), intent(in) :: tauaer(ncol, nlay, nbndsw)
            real(kind=rb), intent(in) :: ssaaer(ncol, nlay, nbndsw)
            real(kind=rb), intent(in) :: asmaer(ncol, nlay, nbndsw)
            real(kind=rb), intent(in) :: ecaer(ncol, nlay, naerec)

            real(kind=rb), intent(out) :: swuflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swdflx(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swhr(ncol, nlay)
            real(kind=rb), intent(out) :: swuflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swdflxc(ncol, nlay + 1)
            real(kind=rb), intent(out) :: swhrc(ncol,nlay)

            real(kind=rb), intent(in), optional :: bndsolvar(nbndsw)
            real(kind=rb), intent(inout), optional :: indsolvar(2)
            real(kind=rb), intent(in), optional :: solcycfrac

            call rrtmg_sw &
                (ncol, nlay, icld, iaer, &
                 play, plev, tlay, tlev, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                 asdir, asdif, aldir, aldif, &
                 coszen, adjes, dyofyr, scon, isolvar, &
                 inflgsw, iceflgsw, liqflgsw, cldfr, &
                 taucld, ssacld, asmcld, fsfcld, &
                 cicewp, cliqwp, reice, reliq, &
                 tauaer, ssaaer, asmaer, ecaer, &
                 swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
! optional I/O
                 bndsolvar, indsolvar, solcycfrac)
          end subroutine rrtmg_sw_nomcica_wrapper

      end module rrtmg_sw_c_binder

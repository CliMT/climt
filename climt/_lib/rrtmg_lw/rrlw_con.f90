      module rrlw_con

      !use parkind, only : rb => kind_rb
      !JOY use c_double from iso_c_binding
      !instead of kind_rb
      use iso_c_binding

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw constants

! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! fluxfac:  real   : radiance to flux conversion factor 
! heatfac:  real   : flux to heating rate conversion factor
!oneminus:  real   : 1.-1.e-6
! pi     :  real   : pi
! grav   :  real   : acceleration of gravity
! planck :  real   : planck constant
! boltz  :  real   : boltzmann constant
! clight :  real   : speed of light
! avogad :  real   : avogadro constant 
! alosmt :  real   : loschmidt constant
! gascon :  real   : molar gas constant
! radcn1 :  real   : first radiation constant
! radcn2 :  real   : second radiation constant
! sbcnst :  real   : stefan-boltzmann constant
!  secdy :  real   : seconds per day  
!------------------------------------------------------------------

      real(kind=c_double) :: fluxfac, heatfac
      real(kind=c_double) :: oneminus, pi, grav
      real(kind=c_double) :: planck, boltz, clight
      real(kind=c_double) :: avogad, alosmt, gascon
      real(kind=c_double) :: radcn1, radcn2
      real(kind=c_double) :: sbcnst, secdy

      contains

          subroutine set_rrtm_constants(py_pi, py_grav, py_planck, &
                  py_boltz, py_clight, py_avogad, py_alosmt, &
                  py_gascon, py_sbcnst,&
                  py_secdy)bind(c, name='rrtmg_set_constants')

              real(kind=c_double) :: py_pi, py_grav, py_planck
              real(kind=c_double) :: py_boltz, py_clight, py_avogad
              real(kind=c_double) :: py_alosmt, py_gascon, py_sbcnst
              real(kind=c_double) :: py_secdy

              pi = py_pi
              grav = py_grav
              planck = py_planck
              boltz = py_boltz
              clight = py_clight
              avogad = py_avogad
              alosmt = py_alosmt
              gascon = py_gascon

              radcn1 = 2.*planck*clight*clight*1.e-07
              radcn2 = planck*clight/boltz

              sbcnst = py_sbcnst
              secdy = py_secdy

          end subroutine set_rrtm_constants

      end module rrlw_con


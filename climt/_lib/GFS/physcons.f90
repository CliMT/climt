!   this module contains some the most frequently used math and        !
!   physics constatns for gcm models.                                  !
!                                                                      !
!   references:                                                        !
!     as set in NMC handbook from Smithsonian tables.                  !
!                                                                      !
          module physcons                
!
  use kinds, only : r_kind
  use iso_c_binding
!
  implicit none
!
  public

!  --- ...  Math constants

  real(r_kind),parameter:: con_pi     =3.1415926535897931 ! pi
  real(r_kind),parameter:: con_sqrt2  =1.414214e+0        ! square root of 2
  real(r_kind),parameter:: con_sqrt3  =1.732051e+0        ! square root of 3

!  --- ...  Geophysics/Astronomy constants

  real(r_kind):: con_rerth  =6.3712e+6      ! radius of earth   (m)
  real(r_kind):: con_g      =9.80665e+0     ! gravity           (m/s2)
  real(r_kind):: con_omega  =7.2921e-5      ! ang vel of earth  (1/s)
  !real(r_kind),parameter:: con_p0     =1.01325e5      ! std atms pressure (pa)
  !real(r_kind),parameter:: con_solr   =1.3660e+3      ! solar constant    (W/m2)-liu(2002)

!  --- ...  Thermodynamics constants

  real(r_kind):: con_rgas   =8.314472       ! molar gas constant  (J/mol/K)
  real(r_kind):: con_rd     =2.8705e+2      ! gas constant air    (J/kg/K)
  real(r_kind):: con_rv     =4.6150e+2      ! gas constant H2O    (J/kg/K)
  real(r_kind):: con_cp     =1.0046e+3      ! spec heat air @p    (J/kg/K)
!  real(r_kind),parameter:: con_cv     =7.1760e+2      ! spec heat air @v    (J/kg/K)
  real(r_kind):: con_cvap   =1.8460e+3      ! spec heat H2O gas   (J/kg/K)
!  real(r_kind),parameter:: con_cliq   =4.1855e+3      ! spec heat H2O liq   (J/kg/K)
!  real(r_kind),parameter:: con_csol   =2.1060e+3      ! spec heat H2O ice   (J/kg/K)
!  real(r_kind),parameter:: con_hvap   =2.5000e+6      ! lat heat H2O cond   (J/kg)
!  real(r_kind),parameter:: con_hfus   =3.3358e+5      ! lat heat H2O fusion (J/kg)
!  real(r_kind),parameter:: con_psat   =6.1078e+2      ! pres at H2O 3pt     (Pa)  
!  real(r_kind),parameter:: con_t0c    =2.7315e+2      ! temp at 0C          (K)
!  real(r_kind),parameter:: con_ttp    =2.7316e+2      ! temp at H2O 3pt     (K)
!  real(r_kind),parameter:: con_tice   =2.7120e+2      ! temp freezing sea     (K)
!  real(r_kind),parameter:: con_jcal   =4.1855E+0      ! joules per calorie  ()
!  real(r_kind),parameter:: con_rhw0   =1022.0         ! sea water reference density (kg/m^3)
!  real(r_kind),parameter:: con_epsq   =1.0E-12        ! min q for computing precip type

!  Secondary constants

  real(r_kind):: con_rocp
  real(r_kind):: con_cpor
  real(r_kind):: con_rog
  real(r_kind):: con_fvirt
  real(r_kind):: con_eps
  !real(r_kind),parameter:: con_epsm1  =con_rd/con_rv-1.
  !real(r_kind),parameter:: con_dldt   =con_cvap-con_cliq
  !real(r_kind),parameter:: con_xpona  =-con_dldt/con_rv
  !real(r_kind),parameter:: con_xponb  =-con_dldt/con_rv+con_hvap/(con_rv*con_ttp)

!  --- ...  Other Physics/Chemistry constants (source: 2002 CODATA)

!  real(r_kind),parameter:: con_c      =2.99792458e+8  ! speed of light      (m/s)
!  real(r_kind),parameter:: con_plnk   =6.6260693e-34  ! planck constatn     (J/s)
!  real(r_kind),parameter:: con_boltz  =1.3806505e-23  ! boltzmann constant  (J/K)
!  real(r_kind),parameter:: con_sbc    =5.670400e-8    ! stefan-boltzmann    (W/m2/K4)
!  real(r_kind),parameter:: con_avgd   =6.0221415e23   ! avogadro constant   (1/mol)
!  real(r_kind),parameter:: con_gasv   =22413.996e-6   ! vol of ideal gas at 273.15k, 101.325kpa (m3/mol)
!  real(r_kind),parameter:: con_amd    =28.9644        ! molecular wght of dry air (g/mol)
!  real(r_kind),parameter:: con_amw    =18.0154        ! molecular wght of water vapor (g/mol)
!  real(r_kind),parameter:: con_amo3   =47.9982        ! molecular wght of o3  (g/mol)
!  real(r_kind),parameter:: con_amco2  =44.011         ! molecular wght of co2 (g/mol)
!  real(r_kind),parameter:: con_amo2   =31.9999        ! molecular wght of o2  (g/mol)
!  real(r_kind),parameter:: con_amch4  =16.043         ! molecular wght of ch4 (g/mol)
!  real(r_kind),parameter:: con_amn2o  =44.013         ! molecular wght of n2o (g/mol)

  contains
    subroutine set_physical_constants(radius, omega, R, Rd, Rv, g, Cp, C_cond)bind(c, name='gfs_set_constants')
        real(c_double), intent(in):: radius, omega, R, Rd, Rv, g, Cp, C_cond

        con_rerth = radius
        con_omega = omega
        con_rgas = R
        con_rd = Rd
        con_rv = Rv
        con_g = g
        con_cp = Cp
        con_cvap = C_cond

        con_rocp = con_rd/con_cp
        con_cpor = con_cp/con_rd
        con_rog = con_rd/con_g
        con_fvirt = con_rv/con_rd-1.
        con_eps = con_rd/con_rv

    end subroutine

end module physcons

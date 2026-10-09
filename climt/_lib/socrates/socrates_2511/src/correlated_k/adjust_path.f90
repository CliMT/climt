! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Temperature and pressure adjustment for line parameters
!
SUBROUTINE adjust_path ( &
!      Input arguments
       mol, iso, wave_num, S, abhw, lse, abcoeff, sbhw, T, p, gf, &
!      Output arguments
       line_centre, S_adj, alpha_lorentz, alpha_lorentz_air, &
       alpha_lorentz_self, alpha_doppler)

! Description:
!
! Calculates line parameters at appropriate temperature and
! pressure from HITRAN values at reference T,p. 
! Computes Lorentz half width and  Doppler width of line.
! 
! Method:
!
! Uses standard formula
!  
  
! Modules used:
!
  USE realtype_rd, ONLY: RealK
  USE hitran_cnst, ONLY: number_isotopes, q296, qcoeff, iso_mass, &
    speed_of_light, atomic_mass_unit, molar_gas_constant, c2

  IMPLICIT NONE

! Subroutine arguments

! Scalar arguments with intent(in):
  INTEGER, Intent(IN) :: mol ! molecule number
  INTEGER, Intent(IN) :: iso ! isotope number

  REAL  (RealK), Intent(IN) :: T
  REAL  (RealK), Intent(IN) :: p
  REAL  (RealK), Intent(IN) :: gf
  REAL  (RealK), Intent(IN) :: S
  REAL  (RealK), Intent(IN) :: lse
  REAL  (RealK), Intent(IN) :: wave_num
  REAL  (RealK), Intent(IN) :: abhw
  REAL  (RealK), Intent(IN) :: abcoeff
  REAL  (RealK), Intent(IN) :: sbhw

! Scalar arguments with intent(out):
  REAL  (RealK), Intent(OUT) :: line_centre
  REAL  (RealK), Intent(OUT) :: S_adj
  REAL  (RealK), Intent(OUT) :: alpha_lorentz
  REAL  (RealK), Intent(OUT) :: alpha_lorentz_air
  REAL  (RealK), Intent(OUT) :: alpha_lorentz_self
  REAL  (RealK), Intent(OUT) :: alpha_doppler

! Local parameters:
  REAL  (RealK), PARAMETER :: T_ref = 296.0     ! K
  REAL  (RealK), PARAMETER :: p_ref = 1.01325E5 ! Pa

! Local scalars:
  REAL  (RealK) :: QT
  REAL  (RealK) :: QT_ref
  REAL  (RealK) :: boltzmann_fact
  REAL  (RealK) :: stim_emiss_fact
  REAL  (RealK) :: abhw_SI
  REAL  (RealK) :: sbhw_SI
  REAL  (RealK) :: lse_SI
  REAL  (RealK) :: S_SI

  INTEGER :: i
  INTEGER :: reqd_species

!- End of header
  
  reqd_species = 0
  IF (mol == 1) THEN
     reqd_species = iso
  ELSE
     DO i = 1, mol-1
        reqd_species = reqd_species + number_isotopes(i)
     END DO
     reqd_species = reqd_species + iso
  END IF

  line_centre = wave_num * 100.0 ! now in m-1

  S_SI = S / &
       (iso_mass(reqd_species)*atomic_mass_unit * 100.0) ! now in m kg-1

  abhw_SI = abhw * 100.0 ! now in m-1/atm
  sbhw_SI = sbhw * 100.0 ! now in m-1/atm

  lse_SI = lse * 100.0   ! now in m-1

  QT     = qcoeff(INT(T), reqd_species) + (T - REAL(INT(T),RealK)) &
           *(qcoeff(INT(T)+1, reqd_species) - qcoeff(INT(T), reqd_species))

  QT_ref = q296(reqd_species)

  boltzmann_fact  = EXP( c2 * lse_SI * (T-T_ref) / (T * T_ref) )

  stim_emiss_fact = (1 - EXP( -c2 * line_centre / T)) / &
                    (1 - EXP( -c2 * line_centre / T_ref))

  S_adj           = S_SI * (QT_ref / QT) * boltzmann_fact * stim_emiss_fact

  alpha_doppler   = (line_centre / speed_of_light) *  &
                    (2.0*LOG(2.0)*molar_gas_constant*T / &
                    (iso_mass(reqd_species)*1.0e-3))**0.5

! The temperature dependence of the air-broadened width is also used for the
! self-broadened width at present.
  alpha_lorentz_air  = (1.0-gf) * abhw_SI * (p/p_ref) * (T_ref/T) ** abcoeff
  alpha_lorentz_self = gf * sbhw_SI * (p/p_ref) * (T_ref/T) ** abcoeff
  alpha_lorentz      = alpha_lorentz_air + alpha_lorentz_self

END SUBROUTINE adjust_path

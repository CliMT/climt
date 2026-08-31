! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to calculate volume fractions for SRA profiles.
!
! Method:
!   The volume fractions are calculated from extinctions
!   at 550 nm.
!
!- ---------------------------------------------------------------------
      SUBROUTINE calc_volume_fraction(n_component
     &  , fraction_component, ext_550nm_component, ext_550nm_aerosol
     &  , volume_fraction
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
      INTEGER, Intent(IN) ::
     &    n_component
!           Number of components
      REAL  (RealK), Intent(IN) ::
     &    ext_550nm_component(n_component)
!           Component extinc. at 550 nm
     &  , ext_550nm_aerosol
!           Aerosol extinction at 550 nm
     &  , fraction_component(n_component)
!           Fraction of each cmpt.
      REAL  (RealK), Intent(OUT) ::
     &    volume_fraction(n_component)
!           Volume fractions of components
!
!     Local variables
      INTEGER
     &    i
!           Loop variable
      REAL  (RealK) ::
     &    inv_sum_ext
!           1/sum of weight extinctions
     &  , sum_volume_fraction
!           Sum of volume fractions
!
!
      inv_sum_ext=0.0_RealK
      DO i=1, n_component
       inv_sum_ext=inv_sum_ext
     &   +ext_550nm_component(i)*fraction_component(i)
      ENDDO
      sum_volume_fraction=ext_550nm_aerosol/inv_sum_ext
      DO i=1, n_component
        volume_fraction(i)=fraction_component(i)*sum_volume_fraction
      ENDDO
!
!
!
      RETURN
      END

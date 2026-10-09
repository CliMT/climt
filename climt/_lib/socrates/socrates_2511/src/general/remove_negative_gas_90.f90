! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to remove neative gaseous p-scalings.
!
SUBROUTINE remove_negative_gas_90 &
(nd_band, nd_species, nd_k_term, nd_scale_variable, &
 n_band, n_band_absorb, &
 index_absorb, type_absorb, &
 i_scale_fnc, i_band_esft, scale_vector &
)
!
! Method:
!   The exponent of the pressure-scaling of gaseous transmission
!   is examined, and if neagtive it is reset to 0.
!
! Modules
  USE realtype_rd
  USE def_std_io_icf
  USE dimensions_spec_ucf
  USE rad_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
  INTEGER, Intent(IN) :: nd_band
!   Size allocated for spectral bands
  INTEGER, Intent(IN) :: nd_species
!   Size allocated for gaseous species
  INTEGER, Intent(IN) :: nd_k_term
!   Size allocated for k-terms
  INTEGER, Intent(IN) :: nd_scale_variable
!   Size allocated for scaling variables
!
  INTEGER, Intent(IN) :: n_band
!   Number of bands
  INTEGER, Intent(IN), Dimension(nd_band) :: n_band_absorb
!   Number of absorbers in bands
  INTEGER, Intent(IN), Dimension(nd_species, nd_band) :: index_absorb
!   List of active absorbs
  INTEGER, Intent(IN), Dimension(nd_band, nd_species) :: i_scale_fnc
!   Scaling functions
  INTEGER, Intent(IN), Dimension(nd_band, nd_species) :: i_band_esft
!   Number of ESFT terms in band
  INTEGER, Intent(IN), Dimension(nd_species) :: type_absorb
!   Types of absorbers
  REAL  (RealK), Intent(INOUT), Dimension(nd_scale_variable, &
    nd_k_term, nd_band, nd_species) :: scale_vector
!   Scaling parameters for gases
!
! Local variables.
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: k
!   Loop variable
  INTEGER :: l
!   Loop variable
!
!
!
! Loop through the bands checking the scaling data depending on
! the scaling function used.
  DO i=1, n_band
    DO k=1, n_band_absorb(i)
      j=index_absorb(k, i)
      IF (i_scale_fnc(i, j) == IP_scale_power_law) THEN
        DO l=1, i_band_esft(i, j)
          IF (scale_vector(1, l, i, j) < 0.0_RealK) THEN
            WRITE(iu_stdout, '(a, i5, /a, i5, 3x, a, i5)') &
              'Resetting pressure scaling in band ', i, &
              'for gas ', type_absorb(j), &
              'and term ', l
            WRITE(iu_stdout, '(a, 1pe10.3)') &
              'Old scaling was ', scale_vector(1, l, i, j)
            scale_vector(1, l, i, j)=0.0_RealK
          ENDIF
        ENDDO
      ELSE IF (i_scale_fnc(i, j) == IP_scale_power_quad) THEN
        DO l=1, i_band_esft(i, j)
          IF (scale_vector(1, l, i, j) < 0.0_RealK) THEN
            WRITE(iu_stdout, '(a, i5, /a, i5, 3x, a, i5)') &
              'Resetting pressure scaling in band ', i, &
              'for gas ', j, 'and term ', l
            WRITE(iu_stdout, '(a, 1pe10.3)') &
               'Old scaling was ', scale_vector(1, l, i, j)
            scale_vector(1, l, i, j)=0.0_RealK
          ENDIF
        ENDDO
      ELSE IF (i_scale_fnc(i, j) == IP_scale_doppler_quad) THEN
        DO l=1, i_band_esft(i, j)
          IF (scale_vector(1, l, i, j) < 0.0_RealK) THEN
            WRITE(iu_stdout, '(a, i5, /a, i5, 3x, a, i5)') &
              'Resetting pressure scaling in band ', i, &
              'for gas ', j, 'and term ', l
            WRITE(iu_stdout, '(a, 1pe10.3)') &
              'Old scaling was ', scale_vector(1, l, i, j)
            scale_vector(1, l, i, j)=0.0_RealK
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
!
!
!
  RETURN
END SUBROUTINE remove_negative_gas_90

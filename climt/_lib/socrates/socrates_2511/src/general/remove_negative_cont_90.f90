! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to remove neative continuum p-scalings.
!
SUBROUTINE remove_negative_cont_90 &
(nd_band, nd_continuum, nd_scale_variable, &
 n_band, n_band_continuum, &
 index_continuum, i_scale_fnc_cont, scale_continuum &
)
!
! Method:
!   The exponent of the pressure-scaling of the continuum is
!   examined, and if neagtive it is reset to 0.
!
! Modules to set types of variables:
  USE realtype_rd
  USE def_std_io_icf
  USE rad_pcf
!
!
  IMPLICIT NONE
!
!
!
! Include header files
!
! Dummy arguments.
  INTEGER, Intent(in) :: nd_band
!   Size allocated for spectral bands
  INTEGER, Intent(in) :: nd_continuum
!   Size allocated for continua
  INTEGER, Intent(in) :: nd_scale_variable
!   Size allocated for scaling variables
!
  INTEGER, Intent(in) :: n_band
!   Number of bands
  INTEGER, Intent(in), Dimension(nd_band) :: n_band_continuum
!   Number of continua in bands
  INTEGER, Intent(in), Dimension(nd_band, nd_continuum) :: &
    index_continuum
!   List of active continua
  INTEGER, Intent(in), Dimension(nd_band, nd_continuum) :: &
    i_scale_fnc_cont
!   Scaling functions for continua
  REAL  (RealK), Intent(out), &
    Dimension(nd_scale_variable, nd_band, nd_continuum) :: &
    scale_continuum
!   Continuum scaling parameters
!
! Local variables.
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: k
!   Loop variable
!
!
! Loop through the bands checking the scaling data depending on
! the scaling function used.
  DO i=1, n_band
    DO k=1, n_band_continuum(i)
      j=index_continuum(i, k)
      SELECT CASE(i_scale_fnc_cont(i, j))
        CASE(IP_scale_power_law)
          IF (scale_continuum(1, i, j) < 0.0_RealK) THEN
            WRITE(iu_stdout, '(a, i5, /a, i5)') &
              'Resetting pressure scaling in band ', i, &
              'for continuum ', j
            WRITE(iu_stdout, '(a, 1pe10.3)') 'Old scaling was ', &
              scale_continuum(1, i, j)
            scale_continuum(1, i, j)=0.0_RealK
          ENDIF
        CASE(IP_scale_power_quad)
          IF (scale_continuum(1, i, j) < 0.0_RealK) THEN
            WRITE(iu_stdout, '(a, i5, /a, i5)') &
              'Resetting pressure scaling in band ', i, &
              'for continuum ', j
            WRITE(iu_stdout, '(a, 1pe10.3)') 'Old scaling was ', &
              scale_continuum(1, i, j)
            scale_continuum(1, i, j)=0.0_RealK
          ENDIF
        CASE(IP_scale_doppler_quad)
          IF (scale_continuum(1, i, j) < 0.0_RealK) THEN
            WRITE(iu_stdout, '(a, i5, /a, i5)') &
              'Resetting pressure scaling in band ', i, &
              'for continuum ', j
            WRITE(iu_stdout, '(a, 1pe10.3)') 'Old scaling was ', &
              scale_continuum(1, i, j)
            scale_continuum(1, i, j)=0.0_RealK
          ENDIF
      END SELECT
    ENDDO
  ENDDO
!
!
!
  RETURN
END SUBROUTINE remove_negative_cont_90

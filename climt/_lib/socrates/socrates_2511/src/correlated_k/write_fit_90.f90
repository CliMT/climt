! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to write a formatted file containing k-fits.
!
SUBROUTINE write_fit_90 &
(iu_k_out, l_continuum, l_cont_gen, l_self_broadening, &
 i_band, i_gas, i_index, i_index_1, i_index_2, &
 p_fit, t_fit, &
 n_points, amount, transmittance, trans_calc, &
 n_k, k, w_k, i_scale, &
 i_scale_function, scale_vector)
!
! Description:
!   This subroutine constructs a formatted file containing a k-fit.
!
! Method:
!   Straightforward.
!
! Modules used.
  USE realtype_rd
  USE rad_pcf, ONLY: n_scale_variable, ip_scale_fnc_null, ip_scale_null, &
                     ip_scale_t_lookup
  USE file_type_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables
  INTEGER, Intent(IN) :: iu_k_out
!   Number of band
!
  INTEGER, Intent(IN) :: i_band
!   Number of band
  INTEGER, Intent(IN) :: i_gas
!   Identifier for absorber
  INTEGER, Intent(IN) :: i_index
!   Index number of absorber
  INTEGER, Intent(IN) :: i_index_1
!   Index of first continuum gas
  INTEGER, Intent(IN) :: i_index_2
!   Index of second continuum gas
  INTEGER, Intent(IN) :: n_points
!   Number of data points
  INTEGER, Intent(IN) :: n_k
!   Number of terms
  INTEGER, Intent(IN) :: i_scale
!   Type of scaling used.
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function used.
  LOGICAL, Intent(IN) :: l_continuum
!   Continuum flag
  LOGICAL, Intent(IN) :: l_cont_gen
!   Generalised continuum flag
  LOGICAL, Intent(IN) :: l_self_broadening
!   Self-broadening flag
  REAL  (RealK), Intent(IN) :: p_fit
!   Pressure of fit
  REAL  (RealK), Intent(IN) :: t_fit
!   Temperature of fit
  REAL  (RealK), Intent(IN), Dimension(:) :: amount
!   Amounts of absorber
  REAL  (RealK), Intent(IN), Dimension(:) :: transmittance
!   Actual transmittances
  REAL  (RealK), Intent(IN), Dimension(:) :: trans_calc
!   Calculated transmittances
  REAL  (RealK), Intent(IN), Dimension(:) :: k
!   ESFT exponents
  REAL  (RealK), Intent(IN), Dimension(:) :: w_k
!   ESFT weights
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scale_vector
!   Scaling parameters
!
! Local variables.
  INTEGER :: i
!   Loop variable
  INTEGER :: l
!   Loop variable
!
!
!
! Output data for this band.
  IF (l_continuum) THEN
    WRITE(iu_k_out, '(a13, i5, //a)') '*FILE TYPE = ', &
      it_file_cont_fit, &
      'Fitted continuum transmissions and extinction.'
    WRITE(iu_k_out, '(6x, a8, i5, a21, i5)') &
      'in band ', i_band, ': for continuum type ', i_index
  ELSE IF (l_cont_gen) THEN
    WRITE(iu_k_out, '(a13, i5, //a)') '*FILE TYPE = ', &
      it_file_cont_gen_fit, &
      'Fitted continuum transmissions and k-terms.'
    WRITE(iu_k_out, '(6x, a8, i5, a18, i5, a5, i5)') &
      'in band ', i_band, ': for gas indices ', i_index_1, &
       ' and ', i_index_2
  ELSE
    IF (l_self_broadening) THEN
      WRITE(iu_k_out, '(a13, i5, //a)') &
        '*FILE TYPE = ', it_file_line_fit_self_id, &
        'Fitted gaseous transmissions and k-terms.'
    ELSE
      WRITE(iu_k_out, '(a13, i5, //a)') &
        '*FILE TYPE = ', it_file_line_fit_id, &
        'Fitted gaseous transmissions and k-terms.'
    END IF
    WRITE(iu_k_out, '(6x, a8, i5, a26, i5)') &
      'in band ', i_band, ': for absorber identifier ', i_gas
  ENDIF
  IF (i_scale_function == ip_scale_t_lookup) THEN
    WRITE(iu_k_out, '(6x, a15, 1pe10.3, a3)') &
      'at temperature ', t_fit, ' K.'
  ELSE
    WRITE(iu_k_out, '(6x, a12, 1pe10.3, a21, 1pe10.3, a3)') &
      'at pressure ', p_fit, ' Pa, and temperature ', t_fit, ' K.'
  END IF
  WRITE(iu_k_out, '(/)')
  WRITE(iu_k_out, '(5x, a18, 7x, a20, 5x, a20, 5x)' ) &
    'Amount of absorber', 'Actual transmittance', &
    'Fitted transmittance'
  IF (l_continuum) THEN
    WRITE(iu_k_out, '(9x, a11)') '(kg.mol)/m5'
  ELSE IF (l_cont_gen) THEN
    WRITE(iu_k_out, '(9x, a11)') 'kg2/m5'
  ELSE
    WRITE(iu_k_out, '(9x, a5)') 'kg/m2'
  ENDIF
  DO i=1, n_points
    WRITE(iu_k_out, '(3(9x, 1pe16.9))' ) &
      amount(i), transmittance(i), trans_calc(i)
  ENDDO
  WRITE(iu_k_out, '(/)')

  IF (i_scale_function == ip_scale_fnc_null) THEN
    WRITE(iu_k_out, '(5x, a18, i5, a20, i5, a21, i5)') &
      'Number of terms = ', n_k, ': Type of scaling = ', &
      ip_scale_null, ': Scaling function = ', i_scale_function
  ELSE
    WRITE(iu_k_out, '(5x, a18, i5, a20, i5, a21, i5)') &
      'Number of terms = ', n_k, ': Type of scaling = ', &
      i_scale, ': Scaling function = ', i_scale_function
  END IF
  WRITE(iu_k_out, '(4x, a6, 14x, a6, 14x, a18)') &
    'k-term', 'Weight', 'Scaling parameters'
  IF (l_continuum) THEN
    WRITE(iu_k_out, '(5x, a11)') 'm5/(kg.mol)'
  ELSE IF (l_cont_gen) THEN
    WRITE(iu_k_out, '(5x, a11)') 'm5/kg2'
  ELSE
    WRITE(iu_k_out, '(5x, a5)') 'm2/kg'
  ENDIF
  DO i=1, n_k
    WRITE(iu_k_out, '(2(3x, 1pe16.9), (t39, 2(3x, 1pe16.9)))') &
        k(i), w_k(i), &
        (scale_vector(l, i), l=1, &
        n_scale_variable(i_scale_function))
  ENDDO
!
  WRITE(iu_k_out, '(/)')
!
!
!
  RETURN
END SUBROUTINE write_fit_90

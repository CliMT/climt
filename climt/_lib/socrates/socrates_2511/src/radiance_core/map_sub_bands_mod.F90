! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Map the contribution of each k-term to each sub-band
!
!------------------------------------------------------------------------------
MODULE map_sub_bands_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE map_sub_bands(Sp)

  USE realtype_rd,  ONLY: RealK
  USE def_spectrum, ONLY: StrSpecData, allocate_spectrum

  IMPLICIT NONE


! Spectral data:
  TYPE(StrSpecData), INTENT(INOUT) :: Sp

  INTEGER :: i, j, i_sub, i_k, i_band, i_gas
  INTEGER :: j_first, j_last, i_k_min, i_k_max, n_k
  INTEGER :: i_sub_band_gas(Sp%Dim%nd_band)
  REAL (RealK) :: recip_total_energy


  ! Initially loop through sub-bands from block 17 to find
  ! the total number of sub-bands for each k-term
  Sp%Map%n_sub_band_k = 0
  ! and the total number of minor gas k-terms for each sub-band
  Sp%Map%n_k_sub_band = 0
  DO i_sub=1, Sp%Var%n_sub_band
    i_band = Sp%Var%index_sub_band(1, i_sub)
    i_k    = Sp%Var%index_sub_band(2, i_sub)
    IF (i_k == 0) THEN
      IF (Sp%Gas%n_band_absorb(i_band) > 0) THEN
        ! In this case there is only one sub-band for the band
        Sp%Map%n_sub_band_k(:,i_band) = 1
        ! All the major gas k-terms will contribute to this sub-band
        i_gas = Sp%Gas%index_absorb(1, i_band)
        Sp%Map%n_k_sub_band(i_gas, i_sub) = Sp%Gas%i_band_k(i_band, i_gas)
      END IF
    ELSE
      ! Otherwise i_k points to the major gas k-term
      i_gas = Sp%Gas%index_absorb(1, i_band)
      ! Increment the number of sub-bands found for this k-term
      Sp%Map%n_sub_band_k(i_k, i_band) &
        = Sp%Map%n_sub_band_k(i_k, i_band) + 1
      ! For the major gas there is only one k-term in each sub-band
      Sp%Map%n_k_sub_band(i_gas, i_sub) = 1
    END IF
    ! Now loop through the minor gases in this band
    DO i=2, Sp%Gas%n_band_absorb(i_band)
      i_gas = Sp%Gas%index_absorb(i, i_band)
      ! Loop through the minor gas sub-bands which have different limits
      IF (Sp%Gas%sub_band_k(1, i_band, i_gas) == 0) THEN
        ! In this case there is only one gas sub-band for the band
        ! All the gas k-terms will contribute to this sub-band
        Sp%Map%n_k_sub_band(i_gas, i_sub) = Sp%Gas%i_band_k(i_band, i_gas)
      ELSE
        j_first = 0
        DO j=1, Sp%Gas%n_sub_band_gas(i_band, i_gas)
          ! Find range of gas sub-bands that overlap with major gas sub-band
          IF (Sp%Gas%wavelength_sub_band(1, j, i_band, i_gas) >= &
              Sp%Var%wavelength_sub_band(2, i_sub)) EXIT
          IF (Sp%Gas%wavelength_sub_band(2, j, i_band, i_gas) > &
              Sp%Var%wavelength_sub_band(1, i_sub)) THEN
            IF (j_first==0) j_first = j
            j_last = j
          END IF
        END DO
        ! Find the number of unique minor gas k-terms used in this range
        Sp%Map%n_k_sub_band(i_gas, i_sub) = 0
        IF (j_first > 0) THEN
          i_k_min = MINVAL(Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas))
          i_k_max = MAXVAL(Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas))
          DO
            Sp%Map%n_k_sub_band(i_gas, i_sub) &
              = Sp%Map%n_k_sub_band(i_gas, i_sub) + 1
            IF (i_k_min >= i_k_max) EXIT
            i_k_min = MINVAL(Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas), &
              MASK=Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas) > i_k_min)
          END DO
        END IF
      END IF
    END DO
  END DO

  ! Reallocate the spectrum type arrays to hold the mapped values
  Sp%Dim%nd_sub_band_k = MAXVAL(Sp%Map%n_sub_band_k)
  Sp%Dim%nd_k_sub_band = MAXVAL(Sp%Map%n_k_sub_band)
  Sp%Dim%nd_k_term_map = Sp%Dim%nd_k_term

  IF (ALLOCATED(Sp%Map%list_sub_band_k)) THEN
    DEALLOCATE(Sp%Map%list_sub_band_k)
  END IF
  ALLOCATE(Sp%Map%list_sub_band_k(Sp%Dim%nd_sub_band_k, &
                                  Sp%Dim%nd_k_term, &
                                  Sp%Dim%nd_band))

  IF (ALLOCATED(Sp%Map%weight_sub_band_k)) THEN
    DEALLOCATE(Sp%Map%weight_sub_band_k)
  END IF
  ALLOCATE(Sp%Map%weight_sub_band_k(Sp%Dim%nd_sub_band_k, &
                                    Sp%Dim%nd_k_term, &
                                    Sp%Dim%nd_band))

  IF (ALLOCATED(Sp%Map%list_k_sub_band)) THEN
    DEALLOCATE(Sp%Map%list_k_sub_band)
  END IF
  ALLOCATE(Sp%Map%list_k_sub_band(Sp%Dim%nd_k_sub_band, &
                                  Sp%Dim%nd_species, &
                                  Sp%Dim%nd_sub_band))

  IF (ALLOCATED(Sp%Map%weight_k_sub_band)) THEN
    DEALLOCATE(Sp%Map%weight_k_sub_band)
  END IF
  ALLOCATE(Sp%Map%weight_k_sub_band(Sp%Dim%nd_k_sub_band, &
                                    Sp%Dim%nd_species, &
                                    Sp%Dim%nd_sub_band))

  IF (ALLOCATED(Sp%Map%weight_k_major)) THEN
    DEALLOCATE(Sp%Map%weight_k_major)
  END IF
  ALLOCATE(Sp%Map%weight_k_major(Sp%Dim%nd_k_term_map, Sp%Dim%nd_species, &
                                 Sp%Dim%nd_k_term_map, Sp%Dim%nd_band))

  ! Now loop through sub-bands to list the sub-bands that are
  ! contributed to by each k-term and the associated weights
  Sp%Map%n_sub_band_k = 0
  i_sub_band_gas = 0
  Sp%Map%list_sub_band_k = 0
  Sp%Map%weight_sub_band_k = 0.0_RealK
  Sp%Map%list_k_sub_band = 0
  Sp%Map%weight_k_sub_band = 0.0_RealK
  Sp%Map%weight_k_major = 0.0_RealK
  DO i_sub=1, Sp%Var%n_sub_band
    i_band = Sp%Var%index_sub_band(1, i_sub)
    i_k    = Sp%Var%index_sub_band(2, i_sub)
    recip_total_energy = 1.0_RealK &
      / ( 1.0_RealK/Sp%Var%wavelength_sub_band(1, i_sub) &
        - 1.0_RealK/Sp%Var%wavelength_sub_band(2, i_sub) )
    IF (i_k == 0) THEN
      ! In this case there is only one sub-band for the band
      Sp%Map%list_sub_band_k(1, :, i_band) = i_sub
      IF (Sp%Gas%n_band_absorb(i_band) > 0) THEN
        i_gas = Sp%Gas%index_absorb(1, i_band)
        Sp%Map%n_sub_band_k(:, i_band) = 1
        Sp%Map%weight_sub_band_k(1, :, i_band) = Sp%Gas%w(:, i_band, i_gas)
        ! All the major gas k-terms will contribute to this sub-band
        DO n_k=1, Sp%Gas%i_band_k(i_band, i_gas)
          Sp%Map%list_k_sub_band(n_k, i_gas, i_sub) = n_k
          Sp%Map%weight_k_sub_band(n_k, i_gas, i_sub) &
            = Sp%Gas%w(n_k, i_band, i_gas)
        END DO
      END IF
    ELSE
      ! Otherwise i_k points to the major gas k-term
      i_gas = Sp%Gas%index_absorb(1, i_band)
      ! Increment the number of sub-bands found for this band
      i_sub_band_gas(i_band) = i_sub_band_gas(i_band) + 1
      ! Increment the number of sub-bands found for this k-term
      Sp%Map%n_sub_band_k(i_k, i_band) &
        = Sp%Map%n_sub_band_k(i_k, i_band) + 1
      ! Add this sub-band to the list for this k-term
      Sp%Map%list_sub_band_k(Sp%Map%n_sub_band_k(i_k, i_band), i_k, i_band) &
        = i_sub
      ! For the major gas the weight is the fractional contribution of
      ! this sub-band to the band
      IF (Sp%Gas%sub_band_k(1, i_band, i_gas) == 0) THEN
        ! No sub-band data for this gas
        ! so assume each k-term maps to a single sub-band
        Sp%Map%weight_sub_band_k(Sp%Map%n_sub_band_k(i_k,i_band), i_k, i_band) &
          = Sp%Gas%w(i_k, i_band, i_gas)
      ELSE
        ! Otherwise use the gas sub-band weight
        Sp%Map%weight_sub_band_k(Sp%Map%n_sub_band_k(i_k,i_band), i_k, i_band) &
          = Sp%Gas%sub_band_w(i_sub_band_gas(i_band), i_band, i_gas)
      END IF
      ! There is only 1 major gas k-term in the list for the sub-band
      Sp%Map%list_k_sub_band(1, i_gas, i_sub) = i_k
      Sp%Map%weight_k_sub_band(1, i_gas, i_sub) = 1.0_RealK
    END IF
    ! Now loop through the minor gases in this band
    DO i=2, Sp%Gas%n_band_absorb(i_band)
      i_gas = Sp%Gas%index_absorb(i, i_band)
      ! Loop through the minor gas sub-bands which have different limits
      IF (Sp%Gas%sub_band_k(1, i_band, i_gas) == 0 .OR. i_k == 0) THEN
        ! In this case there is only one gas sub-band for the band
        ! All the gas k-terms will contribute to this sub-band
        DO n_k=1, Sp%Gas%i_band_k(i_band, i_gas)
          Sp%Map%list_k_sub_band(n_k, i_gas, i_sub) = n_k
          ! For the minor gases, the weight is the fractional
          ! contribution of the k-term to this sub-band.
          ! Here, that is just the k-term weight.
          Sp%Map%weight_k_sub_band(n_k, i_gas, i_sub) &
            = Sp%Gas%w(n_k, i_band, i_gas)
          ! The fraction of the major gas k-term flux that is overlapped
          ! by each minor gas k-term is also just the k-term weight.
          IF (i_k == 0) THEN
            Sp%Map%weight_k_major(n_k, i, :, i_band) &
              = Sp%Gas%w(n_k, i_band, i_gas)
          ELSE
            Sp%Map%weight_k_major(n_k, i, i_k, i_band) &
              = Sp%Gas%w(n_k, i_band, i_gas)
          END IF
        END DO
      ELSE
        j_first = 0
        DO j=1, Sp%Gas%n_sub_band_gas(i_band, i_gas)
          ! Find range of gas sub-bands that overlap with major gas sub-band
          IF (Sp%Gas%wavelength_sub_band(1, j, i_band, i_gas) >= &
              Sp%Var%wavelength_sub_band(2, i_sub)) EXIT
          IF (Sp%Gas%wavelength_sub_band(2, j, i_band, i_gas) > &
              Sp%Var%wavelength_sub_band(1, i_sub)) THEN
            IF (j_first==0) j_first = j
            j_last = j
          END IF
        END DO
        IF (j_first > 0) THEN
          ! Find the number of unique minor gas k-terms used in this range
          i_k_min = MINVAL(Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas))
          i_k_max = MAXVAL(Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas))
          n_k = 0
          DO
            n_k = n_k + 1
            Sp%Map%list_k_sub_band(n_k, i_gas, i_sub) = i_k_min
            ! Here, the weight is the fraction of the sub-band energy that
            ! is overlapped by the minor gas sub-bands for each k-term.
            DO j=j_first, j_last
              IF (Sp%Gas%sub_band_k(j, i_band, i_gas) == i_k_min) THEN
                Sp%Map%weight_k_sub_band(n_k, i_gas, i_sub) &
                  = Sp%Map%weight_k_sub_band(n_k, i_gas, i_sub) &
                  + ( 1.0_RealK &
                      / MAX(Sp%Gas%wavelength_sub_band(1, j, i_band, i_gas), &
                            Sp%Var%wavelength_sub_band(1, i_sub)) &
                    - 1.0_RealK &
                      / MIN(Sp%Gas%wavelength_sub_band(2, j, i_band, i_gas), &
                            Sp%Var%wavelength_sub_band(2, i_sub)) &
                    ) * recip_total_energy
              END IF
            END DO
            ! Find the fraction of the major gas k-term flux that is
            ! overlapped by each minor gas k-term.
            Sp%Map%weight_k_major(i_k_min, i, i_k, i_band) &
              = Sp%Map%weight_k_major(i_k_min, i, i_k, i_band) &
              + Sp%Map%weight_k_sub_band(n_k, i_gas, i_sub) &
              * Sp%Map%weight_sub_band_k(Sp%Map%n_sub_band_k(i_k, i_band), &
                                         i_k, i_band) &
              / Sp%Gas%w(i_k, i_band, Sp%Gas%index_absorb(1, i_band))
            IF (i_k_min >= i_k_max) EXIT
            i_k_min = MINVAL(Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas), &
              MASK=Sp%Gas%sub_band_k(j_first:j_last, i_band, i_gas) > i_k_min)
          END DO
        END IF
      END IF
    END DO
  END DO

  ! Reallocate the gas sub-band arrays to zero size
  Sp%Dim%nd_sub_band_gas = 0
  DEALLOCATE(Sp%Gas%sub_band_k)
  DEALLOCATE(Sp%Gas%sub_band_w)
  DEALLOCATE(Sp%Gas%wavelength_sub_band)
  CALL allocate_spectrum(Sp)

END SUBROUTINE map_sub_bands
END MODULE map_sub_bands_mod

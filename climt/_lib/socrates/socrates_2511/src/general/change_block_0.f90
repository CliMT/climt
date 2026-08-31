! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change the number of bands in a spectral file
!
SUBROUTINE change_block_0(Sp)

  USE def_spectrum
  USE def_std_io_icf

  IMPLICIT NONE

  TYPE (StrSpecData), INTENT(INOUT) :: Sp
!   Spectral file to be assigned

  INTEGER :: ios
  INTEGER :: fb, lb, nb
  INTEGER :: i, n_sub_band, n_rayleigh_coeff
  
  WRITE(iu_stdout, '(/A/)') 'Enter range of spectral bands to keep:' 
  DO
    READ(iu_stdin, *, IOSTAT=ios) fb, lb
    IF ( (ios == 0) .AND. &
         (fb > 0) .AND. &
         (lb <= Sp%Basic%n_band) ) EXIT
    WRITE(iu_err, '(A)') 'Invalid input: Please re-enter.'
  END DO
  nb = lb-fb+1
  Sp%Basic%n_band = nb

  ! Re-order arrays as necessary:
  ! Basic
  Sp%Basic%wavelength_long(1:nb) = Sp%Basic%wavelength_long(fb:lb)
  Sp%Basic%wavelength_short(1:nb) = Sp%Basic%wavelength_short(fb:lb)
  Sp%Basic%n_band_exclude(1:nb) = Sp%Basic%n_band_exclude(fb:lb)
  Sp%Basic%index_exclude(:,1:nb) = Sp%Basic%index_exclude(:,fb:lb)
  ! Solar
  IF (ALLOCATED(Sp%Solar%solar_flux_band)) &
    Sp%Solar%solar_flux_band(1:nb) = Sp%Solar%solar_flux_band(fb:lb)
  IF (ALLOCATED(Sp%Solar%solar_flux_band_ses)) &
    Sp%Solar%solar_flux_band_ses(:,1:nb) = &
    Sp%Solar%solar_flux_band_ses(:,fb:lb)
  IF (ALLOCATED(Sp%Solar%weight_blue)) &
    Sp%Solar%weight_blue(1:nb) = Sp%Solar%weight_blue(fb:lb)
  ! Rayleigh
  IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff)) &
    Sp%Rayleigh%rayleigh_coeff(1:nb) = Sp%Rayleigh%rayleigh_coeff(fb:lb)
  IF (ALLOCATED(Sp%Rayleigh%rayleigh_coeff_gas)) &
    Sp%Rayleigh%rayleigh_coeff_gas(:,1:nb) = &
    Sp%Rayleigh%rayleigh_coeff_gas(:,fb:lb)
  ! Gas
  IF (ALLOCATED(Sp%Gas%n_band_absorb)) &
    Sp%Gas%n_band_absorb(1:nb) = Sp%Gas%n_band_absorb(fb:lb)
  IF (ALLOCATED(Sp%Gas%index_absorb)) &
    Sp%Gas%index_absorb(:,1:nb) = Sp%Gas%index_absorb(:,fb:lb)
  IF (ALLOCATED(Sp%Gas%n_mix_gas)) &
    Sp%Gas%n_mix_gas(1:nb) = Sp%Gas%n_mix_gas(fb:lb)
  IF (ALLOCATED(Sp%Gas%num_mix)) &
    Sp%Gas%num_mix(1:nb) = Sp%Gas%num_mix(fb:lb)
  IF (ALLOCATED(Sp%Gas%mix_gas_band)) &
    Sp%Gas%mix_gas_band(1:nb) = Sp%Gas%mix_gas_band(fb:lb)
  IF (ALLOCATED(Sp%Gas%num_ref_p)) &
    Sp%Gas%num_ref_p(:,1:nb) = Sp%Gas%num_ref_p(:,fb:lb)
  IF (ALLOCATED(Sp%Gas%num_ref_t)) &
    Sp%Gas%num_ref_t(:,1:nb) = Sp%Gas%num_ref_t(:,fb:lb)
  IF (ALLOCATED(Sp%Gas%i_band_k)) &
    Sp%Gas%i_band_k(1:nb,:) = Sp%Gas%i_band_k(fb:lb,:)
  IF (ALLOCATED(Sp%Gas%i_band_k_ses)) &
    Sp%Gas%i_band_k_ses(1:nb) = Sp%Gas%i_band_k_ses(fb:lb)
  IF (ALLOCATED(Sp%Gas%i_scale_k)) &
    Sp%Gas%i_scale_k(1:nb,:) = Sp%Gas%i_scale_k(fb:lb,:)
  IF (ALLOCATED(Sp%Gas%i_scale_fnc)) &
    Sp%Gas%i_scale_fnc(1:nb,:) = Sp%Gas%i_scale_fnc(fb:lb,:)
  IF (ALLOCATED(Sp%Gas%i_scat)) &
    Sp%Gas%i_scat(:,1:nb,:) = Sp%Gas%i_scat(:,fb:lb,:)
  IF (ALLOCATED(Sp%Gas%k)) &
    Sp%Gas%k(:,1:nb,:) = Sp%Gas%k(:,fb:lb,:)
  IF (ALLOCATED(Sp%Gas%w)) &
    Sp%Gas%w(:,1:nb,:) = Sp%Gas%w(:,fb:lb,:)
  IF (ALLOCATED(Sp%Gas%scale)) &
    Sp%Gas%scale(:,:,1:nb,:) = Sp%Gas%scale(:,:,fb:lb,:)
  IF (ALLOCATED(Sp%Gas%p_ref)) &
    Sp%Gas%p_ref(:,1:nb) = Sp%Gas%p_ref(:,fb:lb)
  IF (ALLOCATED(Sp%Gas%t_ref)) &
    Sp%Gas%t_ref(:,1:nb) = Sp%Gas%t_ref(:,fb:lb)
  IF (ALLOCATED(Sp%Gas%k_lookup)) &
    Sp%Gas%k_lookup(:,:,:,:,1:nb) = Sp%Gas%k_lookup(:,:,:,:,fb:lb)
  IF (ALLOCATED(Sp%Gas%k_lookup_sb)) &
    Sp%Gas%k_lookup_sb(:,:,:,:,:,1:nb) = Sp%Gas%k_lookup_sb(:,:,:,:,:,fb:lb)
  IF (ALLOCATED(Sp%Gas%w_ses)) &
    Sp%Gas%w_ses(:,1:nb) = Sp%Gas%w_ses(:,fb:lb)
  IF (ALLOCATED(Sp%Gas%f_mix)) &
    Sp%Gas%f_mix(1:nb) = Sp%Gas%f_mix(fb:lb)
  IF (ALLOCATED(Sp%Gas%n_sub_band_gas)) &
    Sp%Gas%n_sub_band_gas(1:nb,:) = Sp%Gas%n_sub_band_gas(fb:lb,:)
  IF (ALLOCATED(Sp%Gas%sub_band_k)) &
    Sp%Gas%sub_band_k(:,1:nb,:) = Sp%Gas%sub_band_k(:,fb:lb,:)
  IF (ALLOCATED(Sp%Gas%sub_band_w)) &
    Sp%Gas%sub_band_w(:,1:nb,:) = Sp%Gas%sub_band_w(:,fb:lb,:)
  IF (ALLOCATED(Sp%Gas%wavelength_sub_band)) &
    Sp%Gas%wavelength_sub_band(:,:,1:nb,:) = &
    Sp%Gas%wavelength_sub_band(:,:,fb:lb,:)
  ! Planck
  IF (ALLOCATED(Sp%Planck%thermal_coeff)) &
    Sp%Planck%thermal_coeff(:,1:nb) = Sp%Planck%thermal_coeff(:,fb:lb)
  ! Cont
  IF (ALLOCATED(Sp%Cont%n_band_continuum)) &
    Sp%Cont%n_band_continuum(1:nb) = Sp%Cont%n_band_continuum(fb:lb)
  IF (ALLOCATED(Sp%Cont%index_continuum)) &
    Sp%Cont%index_continuum(1:nb,:) = Sp%Cont%index_continuum(fb:lb,:)
  IF (ALLOCATED(Sp%Cont%i_scale_fnc_cont)) &
    Sp%Cont%i_scale_fnc_cont(1:nb,:) = Sp%Cont%i_scale_fnc_cont(fb:lb,:)
  IF (ALLOCATED(Sp%Cont%k_cont)) &
    Sp%Cont%k_cont(1:nb,:) = Sp%Cont%k_cont(fb:lb,:)
  IF (ALLOCATED(Sp%Cont%scale_cont)) &
    Sp%Cont%scale_cont(:,1:nb,:) = Sp%Cont%scale_cont(:,fb:lb,:)
  IF (ALLOCATED(Sp%Cont%p_ref_cont)) &
    Sp%Cont%p_ref_cont(:,1:nb) = Sp%Cont%p_ref_cont(:,fb:lb)
  IF (ALLOCATED(Sp%Cont%t_ref_cont)) &
    Sp%Cont%t_ref_cont(:,1:nb) = Sp%Cont%t_ref_cont(:,fb:lb)
  IF (ALLOCATED(Sp%Cont%k_cont_ses)) &
    Sp%Cont%k_cont_ses(:,:,1:nb,:) = Sp%Cont%k_cont_ses(:,:,fb:lb,:)
  IF (ALLOCATED(Sp%Cont%k_h2oc)) &
    Sp%Cont%k_h2oc(:,:,:,1:nb) = Sp%Cont%k_h2oc(:,:,:,fb:lb)
  ! Generalised continuum
  IF (ALLOCATED(Sp%ContGen%n_band_cont)) &
    Sp%ContGen%n_band_cont(1:nb) = Sp%ContGen%n_band_cont(fb:lb)
  IF (ALLOCATED(Sp%ContGen%index_cont)) &
    Sp%ContGen%index_cont(:,1:nb) = Sp%ContGen%index_cont(:,fb:lb)
  IF (ALLOCATED(Sp%ContGen%i_band_k_cont)) &
    Sp%ContGen%i_band_k_cont(1:nb,:) = Sp%ContGen%i_band_k_cont(fb:lb,:)
  IF (ALLOCATED(Sp%ContGen%i_cont_overlap_band)) &
    Sp%ContGen%i_cont_overlap_band(1:nb,:) = &
    Sp%ContGen%i_cont_overlap_band(fb:lb,:)
  IF (ALLOCATED(Sp%ContGen%i_scat_cont)) &
    Sp%ContGen%i_scat_cont(:,1:nb,:) = Sp%ContGen%i_scat_cont(:,fb:lb,:)
  IF (ALLOCATED(Sp%ContGen%l_cont_major)) &
    Sp%ContGen%l_cont_major(1:nb) = Sp%ContGen%l_cont_major(fb:lb)
  IF (ALLOCATED(Sp%ContGen%k_cont)) &
    Sp%ContGen%k_cont(:,1:nb,:) = Sp%ContGen%k_cont(:,fb:lb,:)
  IF (ALLOCATED(Sp%ContGen%w_cont)) &
    Sp%ContGen%w_cont(:,1:nb,:) = Sp%ContGen%w_cont(:,fb:lb,:)
  IF (ALLOCATED(Sp%ContGen%k_lookup_cont)) &
    Sp%ContGen%k_lookup_cont(:,:,:,1:nb) = &
    Sp%ContGen%k_lookup_cont(:,:,:,fb:lb)
  ! Drop
  IF (ALLOCATED(Sp%Drop%parm_list)) &
    Sp%Drop%parm_list(:,1:nb,:) = Sp%Drop%parm_list(:,fb:lb,:)
  ! Aerosol
  IF (ALLOCATED(Sp%Aerosol%abs)) &
    Sp%Aerosol%abs(:,:,1:nb) = Sp%Aerosol%abs(:,:,fb:lb)
  IF (ALLOCATED(Sp%Aerosol%scat)) &
    Sp%Aerosol%scat(:,:,1:nb) = Sp%Aerosol%scat(:,:,fb:lb)
  IF (ALLOCATED(Sp%Aerosol%phf_fnc)) &
    Sp%Aerosol%phf_fnc(:,:,:,1:nb) = Sp%Aerosol%phf_fnc(:,:,:,fb:lb)
  ! Ice
  IF (ALLOCATED(Sp%Ice%parm_list)) &
    Sp%Ice%parm_list(:,1:nb,:) = Sp%Ice%parm_list(:,fb:lb,:)
  ! Spectral variability
  IF (ALLOCATED(Sp%Var%index_sub_band)) THEN
    n_sub_band = 0
    n_rayleigh_coeff = Sp%Var%n_rayleigh_coeff
    DO i = 1, Sp%Var%n_sub_band
      IF (Sp%Var%index_sub_band(1,i) < fb) THEN
        n_rayleigh_coeff=MAX(n_rayleigh_coeff-1,0)
      ELSE IF (Sp%Var%index_sub_band(1,i) <= lb) THEN
        n_sub_band=n_sub_band+1
        Sp%Var%index_sub_band(1,n_sub_band) = Sp%Var%index_sub_band(1,i)-fb+1
        Sp%Var%index_sub_band(2,n_sub_band) = Sp%Var%index_sub_band(2,i)
        Sp%Var%wavelength_sub_band(:,n_sub_band)=Sp%Var%wavelength_sub_band(:,i)
        Sp%Var%solar_flux_sub_band(n_sub_band,:)=Sp%Var%solar_flux_sub_band(i,:)
        IF (n_rayleigh_coeff >= n_sub_band) THEN
          Sp%Var%rayleigh_coeff(n_sub_band,:) = Sp%Var%rayleigh_coeff(i,:)
        ELSE
          Sp%Var%rayleigh_coeff(n_sub_band,0) = Sp%Var%rayleigh_coeff(i,0)
        END IF
      END IF
    END DO
    Sp%Var%n_sub_band = n_sub_band
    Sp%Var%n_rayleigh_coeff = MIN(n_rayleigh_coeff,n_sub_band)
  END IF

END SUBROUTINE change_block_0

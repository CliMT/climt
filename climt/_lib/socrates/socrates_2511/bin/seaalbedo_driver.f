! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      subroutine seaalbedo_driver(zen_0, n_band, 
     &     wave_length_short, wave_length_long, sea_albedo)

!     THIS ROUTINE IS AN INTERFACE TO SEAALBEDO
!     THAT CALCULATES SURFACE ALBEDO OF A SEA SURFACE
!     AS A FUNCTION OF WAVELENGTH AND SOLAR ZENITH ANGLE
     
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_spec_ucf

      IMPLICIT NONE

      REAL  (RealK) ::
     &    wave_length_short(npd_band)
!           Shorter wavelength limits
     &  , wave_length_long(npd_band)
!           Longer wavelength limits

      REAL  (RealK), Intent(IN) ::
     &    zen_0
!           Secants or cosines of solar zenith angles

      INTEGER, Intent(IN) ::
     &    n_band      
!           Allowed size for spectral bands

      REAL (RealK), Intent(OUT) ::
     &    sea_albedo(n_band) 

      REAL (RealK) ::
     &   wavelength(n_band), szen

      INTEGER iwv, i_band

         szen = acos(1./zen_0)

         DO i_band = 1, n_band 
!        Calculate the central wavelength of the band
!
            wavelength(i_band) = (wave_length_short(i_band) +
     &           wave_length_long(i_band)) /2.
            
            call  seaalbedo(szen,wavelength(i_band),sea_albedo(i_band))
         Enddo
         
         end
      
      include 'seaalbedo.f'

c     This program makes two files containing the refractive index 
c     required for the CLASSIC fresh and aged biomass burning aerosols.
c     It combines the data from Bond and Bergstrom (2006) (for the BC 
c     component), with the default value of 1.53 + 0i assumed for the OC
c     component. The use of Bond & Bergstrom (2006) replaces WCP (1983)
c     and will make the biomass aerosol species more absorbing.
c
c     Volume fractions of BC and OC are calculated based
c     on the mass fractions traditionally assumed in CLASSIC but with
c     a revised density reduced from 1900 to 1500 kgm-3.
c     This change also increases absorption by 1.9/1.5.
c
c     The output files are to be used by Cscatter in the script
c     make_soot_biomass_specfile_blocks
c     
c     Author = Ben Johnson, Jan 2017.
c
c     This can be compiled with nagfor

      Program refractmaker

      IMPLICIT NONE

      REAL wl                               ! wavelength
      
      REAL rr_mix, ri_mix, ri_soot, rr_soot ! refractive index
      REAL bcf_mass_fr, bcf_mass_ag         ! BC mass fraction
      REAL bcf_vol_fr, bcf_vol_ag           ! BC volume fraction
      REAL density_BC, density_BB           ! Aerosol densities (kgm-3)

c     Refractive index for OC component as in Haywood et al. (2003)
c     and as used in CLASSIC since HadGEM1
      real, parameter :: ri_org = 0.
      real, parameter :: rr_org = 1.53 

      INTEGER i    ! loop variable

c     Open file with RI of the BC component from
c     the Bond & Bergstrom (2006) middle value.
      open(1,file='refract_soot_bond')

c     Output files to be created by this program
      open(2,file='refract_biomass_fresh_GA8')
      open(3,file='refract_biomass_aged_GA8')

c     Write out the necessary header information into those files
      write(2,100)
      write(2,101)
      write(2,102)

      write(3,110)
      write(3,101)
      write(3,102)

c     Calculate the volume fraction of BC in the aerosol
c     based on the mass fraction (assumed in CLASSIC based
c     SAFARI-2000 observations from Haywood et al. 2003)

      bcf_mass_fr = 0.0875 ! as in CLASSIC since HadGEM1
      bcf_mass_ag = 0.0540 ! as in CLASSIC since HadGEM1

      density_BC  = 1500.  ! kgm-3 for BC component now set to be same
                           ! as in GLOMAP-mode (including GA7/7.1)

      density_BB  = 1350.  ! Kgm-3 for overall mixture remains same
                           ! as in CLASSIC since HadGEM1

      bcf_vol_fr = bcf_mass_fr * density_BB / density_BC  
      bcf_vol_ag = bcf_mass_ag * density_BB / density_BC 

      do I = 1, 62
         read(1,*) wl,rr_soot, ri_soot
         rr_mix = bcf_vol_fr * rr_soot + (1.0-bcf_vol_fr) * rr_org
         ri_mix = bcf_vol_fr * ri_soot + (1.0-bcf_vol_fr) * ri_org
         write(2,99) wl, rr_mix, ri_mix
         rr_mix = bcf_vol_ag * rr_soot + (1.0-bcf_vol_ag) * rr_org
         ri_mix = bcf_vol_ag * ri_soot + (1.0-bcf_vol_ag) * ri_org
         write(3,99) wl, rr_mix, ri_mix
      end do

      write(2,103)
      write(2,104)
      write(3,103)
      write(3,104)

 98   Format(5x,'$ (',E11.5,'),' )
 99   Format(3(5x,E11.5))
 100  Format('Refractive Index based of fresh BB - mix of BC and OC')
 110  Format('Refractive Index based of fresh BB - mix of BC and OC')
 101  Format('     Wavelength (m)  Real Part       Imaginary Part')
 102  Format('*BEGIN_DATA')
 103  Format('*END')
 104  Format(' ')

      end

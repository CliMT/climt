! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Simple driver to test the Runes interface to the two-stream code

module runes_driver_mod

  implicit none
  public
  contains

  subroutine runes_driver( &
    g, &
    r_gas_dry, &
    cp_air_dry, &
    seconds_per_day, &
    h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
    cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
    he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
    n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
    hdo, hcl, hf, cosso, tosso, yosos, &
    ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
    hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2, &
    zenith_angle, &
    solar_irrad, &
    grey_albedo_sw, &
    grey_albedo_lw, &
    ground_temp, &
    p_layer, t_layer, p_level, t_level, &
    flux_ground, &
    n_layer, n_band, &
    lw_flux_down, lw_flux_up, lw_heating_rate, &
    sw_flux_down, sw_flux_up, sw_heating_rate, &
    water_soluble, dust_like, oceanic, soot, ash, sulphuric, &
    ammonium_sulphate, saharan_dust, &
    accum_sulphate, aitken_sulphate, &
    fresh_soot, aged_soot, &
    sodium_chloride, seasalt_film, seasalt_jet, &
    dust_div1, dust_div2, dust_div3, &
    dust_div4, dust_div5, dust_div6, &
    biomass_1, biomass_2, &
    biogenic, &
    ocff_fresh, ocff_aged, &
    delta, murk, &
    nitrate, &
    twobindust_1, twobindust_2, &
    cloud_frac, liq_mmr, ice_mmr, liq_dim, ice_dim, &
    mode, l_invert)
  
      use socrates_set_spectrum, only: set_spectrum
      use socrates_runes, only: runes, StrDiag, &
                                ip_source_illuminate, ip_source_thermal, &
                                ip_cloud_representation_off, &
                                ip_cloud_representation_combine_ice_water, &
                                ip_overlap_max_random
      use realtype_rd, only: RealExt
    
      implicit none

      integer, intent(in) :: n_layer, n_band
      integer, parameter :: n_profile = 1

      real(RealExt), intent(in) :: g
      real(RealExt), intent(in) :: r_gas_dry
      real(RealExt), intent(in) :: cp_air_dry
      real(RealExt), parameter :: pi = 4.0*atan(1.0)
      real(RealExt), intent(in) :: seconds_per_day

      integer, intent(in) :: mode
      logical, intent(in) :: l_invert

      type(StrDiag) :: sw_diag, lw_diag
      real(RealExt), target, optional :: sw_heating_rate(n_profile, n_layer)
      real(RealExt), target, optional :: lw_heating_rate(n_profile, n_layer)
      real(RealExt), target, optional :: sw_flux_up(n_profile, 0:n_layer)
      real(RealExt), target, optional :: sw_flux_down(n_profile, 0:n_layer)
      real(RealExt), target, optional :: lw_flux_up(n_profile, 0:n_layer)
      real(RealExt), target, optional :: lw_flux_down(n_profile, 0:n_layer)
    
      real(RealExt) :: p_layer(n_profile, n_layer)    
      real(RealExt) :: t_layer(n_profile, n_layer)    
      real(RealExt) :: p_level(n_profile, 0:n_layer)    
      real(RealExt) :: t_level(n_profile, 0:n_layer)    
      real(RealExt) :: ground_temp
      real(RealExt) :: t_ground(n_profile)

      ! Gases
      real(RealExt) :: h2o(n_profile, n_layer)
      real(RealExt) :: co2(n_profile, n_layer)
      real(RealExt) :: o3(n_profile, n_layer)
      real(RealExt) :: n2o(n_profile, n_layer)
      real(RealExt) :: co(n_profile, n_layer)
      real(RealExt) :: ch4(n_profile, n_layer)
      real(RealExt) :: o2(n_profile, n_layer)
      real(RealExt) :: no(n_profile, n_layer)
      real(RealExt) :: so2(n_profile, n_layer)
      real(RealExt) :: no2(n_profile, n_layer)
      real(RealExt) :: nh3(n_profile, n_layer)
      real(RealExt) :: hno3(n_profile, n_layer)
      real(RealExt) :: n2(n_profile, n_layer)
      real(RealExt) :: cfc11(n_profile, n_layer)
      real(RealExt) :: cfc12(n_profile, n_layer)
      real(RealExt) :: cfc113(n_profile, n_layer)
      real(RealExt) :: hcfc22(n_profile, n_layer)
      real(RealExt) :: hfc125(n_profile, n_layer)
      real(RealExt) :: hfc134a(n_profile, n_layer)
      real(RealExt) :: cfc114(n_profile, n_layer)
      real(RealExt) :: tio(n_profile, n_layer)
      real(RealExt) :: vo(n_profile, n_layer)
      real(RealExt) :: h2(n_profile, n_layer)
      real(RealExt) :: he(n_profile, n_layer)
      real(RealExt) :: ocs(n_profile, n_layer)
      real(RealExt) :: na(n_profile, n_layer)
      real(RealExt) :: k(n_profile, n_layer)
      real(RealExt) :: feh(n_profile, n_layer)
      real(RealExt) :: crh(n_profile, n_layer)
      real(RealExt) :: li(n_profile, n_layer)
      real(RealExt) :: rb(n_profile, n_layer)
      real(RealExt) :: cs(n_profile, n_layer)
      real(RealExt) :: ph3(n_profile, n_layer)
      real(RealExt) :: c2h2(n_profile, n_layer)
      real(RealExt) :: hcn(n_profile, n_layer)
      real(RealExt) :: h2s(n_profile, n_layer)
      real(RealExt) :: ar(n_profile, n_layer)
      real(RealExt) :: o(n_profile, n_layer)
      real(RealExt) :: n(n_profile, n_layer)
      real(RealExt) :: no3(n_profile, n_layer)
      real(RealExt) :: n2o5(n_profile, n_layer)
      real(RealExt) :: hono(n_profile, n_layer)
      real(RealExt) :: ho2no2(n_profile, n_layer)
      real(RealExt) :: h2o2(n_profile, n_layer)
      real(RealExt) :: c2h6(n_profile, n_layer)
      real(RealExt) :: ch3(n_profile, n_layer)
      real(RealExt) :: h2co(n_profile, n_layer)
      real(RealExt) :: ho2(n_profile, n_layer)
      real(RealExt) :: hdo(n_profile, n_layer)
      real(RealExt) :: hcl(n_profile, n_layer)
      real(RealExt) :: hf(n_profile, n_layer)
      real(RealExt) :: cosso(n_profile, n_layer)
      real(RealExt) :: tosso(n_profile, n_layer)
      real(RealExt) :: yosos(n_profile, n_layer)
      real(RealExt) :: ch3cho(n_profile, n_layer)
      real(RealExt) :: ch3ooh(n_profile, n_layer)
      real(RealExt) :: ch3coch3(n_profile, n_layer)
      real(RealExt) :: ch3cocho(n_profile, n_layer)
      real(RealExt) :: chocho(n_profile, n_layer)
      real(RealExt) :: c2h5cho(n_profile, n_layer)
      real(RealExt) :: hoch2cho(n_profile, n_layer)
      real(RealExt) :: c2h5coch3(n_profile, n_layer)
      real(RealExt) :: mvk(n_profile, n_layer)
      real(RealExt) :: macr(n_profile, n_layer)
      real(RealExt) :: pan(n_profile, n_layer)
      real(RealExt) :: ch3ono2(n_profile, n_layer)
      
      ! Clouds
      real(RealExt), intent(in) :: cloud_frac(n_profile, n_layer)
      real(RealExt), intent(in) :: liq_mmr(n_profile, n_layer)
      real(RealExt), intent(in) :: ice_mmr(n_profile, n_layer)
      real(RealExt), intent(in) :: liq_dim(n_profile, n_layer)
      real(RealExt), intent(in) :: ice_dim(n_profile, n_layer)
      
      integer :: cloud_rep
    
      real(RealExt) :: zenith_angle
      real(RealExt) :: cz_angle(n_profile)
      real(RealExt) :: solar_irrad
      real(RealExt) :: irrad(n_profile)
      real(RealExt) :: grey_albedo_sw
      real(RealExt) :: grey_albedo_lw

      ! Aerosols
      real(RealExt), intent(in), dimension(n_profile, n_layer) :: &
        water_soluble, dust_like, oceanic, soot, ash, sulphuric, &
        ammonium_sulphate, saharan_dust, &
        accum_sulphate, aitken_sulphate, &
        fresh_soot, aged_soot, &
        sodium_chloride, seasalt_film, seasalt_jet, &
        dust_div1, dust_div2, dust_div3, &
        dust_div4, dust_div5, dust_div6, &
        biomass_1, biomass_2, &
        biogenic, &
        ocff_fresh, ocff_aged, &
        delta, murk, &
        nitrate, &
        twobindust_1, twobindust_2
            
      real(RealExt), dimension(n_profile, n_layer) :: d_mass, density
      real(RealExt), dimension(n_profile, n_layer) :: layer_heat_capacity
      real(RealExt), dimension(n_profile, n_layer) :: liq_frac, ice_frac
    
      real(RealExt) :: flux_ground(n_profile, n_band)
    
      integer :: i, l

      t_ground = ground_temp
      cz_angle = cos(zenith_angle)
      irrad = solar_irrad

      do i=1, n_layer
        do l=1, n_profile
          d_mass(l, i) = abs(p_level(l, i)-p_level(l, i-1))/g
          density(l, i) = p_layer(l, i)/(r_gas_dry*t_layer(l, i))
          layer_heat_capacity(l, i) = d_mass(l, i)*cp_air_dry
          if (liq_mmr(l, i) + ice_mmr(l, i) > 1.0e-12_RealExt) then
            liq_frac(l, i) = liq_mmr(l, i) / (liq_mmr(l, i) + ice_mmr(l, i))
            ice_frac(l, i) = ice_mmr(l, i) / (liq_mmr(l, i) + ice_mmr(l, i))
          else
            liq_frac(l, i) = 0.0_RealExt
            ice_frac(l, i) = 0.0_RealExt
          end if
        end do
      end do

      cloud_rep = ip_cloud_representation_off
      if (maxval(cloud_frac) > 0.0 .and. maxval(liq_mmr + ice_mmr) > 1.0e-12_RealExt) then
        cloud_rep = ip_cloud_representation_combine_ice_water
      end if

      IF (mode .EQ. 1) THEN
        
        ! LW call

        lw_diag%heating_rate => lw_heating_rate
        lw_diag%flux_up => lw_flux_up
        lw_diag%flux_down => lw_flux_down
        call runes( &
          n_profile = n_profile, &
          n_layer = n_layer, &
          diag = lw_diag, &
          spectrum_name = 'lw', &
          i_source = ip_source_thermal, &
          i_cloud_representation = cloud_rep, &
          i_overlap = ip_overlap_max_random, &
          p_layer = p_layer, &
          t_layer = t_layer, &
          t_level = t_level, &
          mass = d_mass, &
          density = density, &
          layer_heat_capacity = layer_heat_capacity, &
          h2o = h2o, &
          co2 = co2, &
          o3 = o3, &
          n2o = n2o, &
          co = co, &
          ch4 = ch4, &
          o2 = o2, &
          no = no, &
          so2 = so2, &
          no2 = no2, &
          nh3 = nh3, &
          hno3 = hno3, &
          n2 = n2, &
          cfc11 = cfc11, &
          cfc12 = cfc12, &
          cfc113 = cfc113, &
          hcfc22 = hcfc22, &
          hfc125 = hfc125, &
          hfc134a = hfc134a, &
          cfc114 = cfc114, &
          tio = tio, &
          vo = vo, &
          h2 = h2, &
          he = he, &
          ocs = ocs, &
          na = na, &
          k = k, &
          feh = feh, &
          crh = crh, &
          li = li, &
          rb = rb, &
          cs = cs, &
          ph3 = ph3, &
          c2h2 = c2h2, &
          hcn = hcn, &
          h2s = h2s, &
          ar = ar, &
          o = o, &
          n = n, &
          no3 = no3, &
          n2o5 = n2o5, &
          hono = hono, &
          ho2no2 = ho2no2, &
          h2o2 = h2o2, &
          c2h6 = c2h6, &
          ch3 = ch3, &
          h2co = h2co, &
          ho2 = ho2, &
          hdo = hdo, &
          hcl = hcl, &
          hf = hf, &
          cosso = cosso, &
          tosso = tosso, &
          yosos = yosos, &
          ch3cho = ch3cho, &
          ch3ooh = ch3ooh, &
          ch3coch3 = ch3coch3, &
          ch3cocho = ch3cocho, &
          chocho = chocho, &
          c2h5cho = c2h5cho, &
          hoch2cho = hoch2cho, &
          c2h5coch3 = c2h5coch3, &
          mvk = mvk, &
          macr = macr, &
          pan = pan, &
          ch3ono2 = ch3ono2, &
          cloud_frac = cloud_frac, &
          liq_frac = liq_frac, &
          ice_frac = ice_frac, &
          liq_mmr = liq_mmr, &
          ice_mmr = ice_mmr, &
          liq_dim = liq_dim, &
          ice_dim = ice_dim, &
          i_st_water = 5, &
          i_st_ice = 8, &
          l_flux_ground = .false., &
          t_ground = t_ground, &
          flux_ground = flux_ground, &
          cos_zenith_angle = cz_angle, &
          solar_irrad = irrad, &
          l_grey_albedo = .true., &
          grey_albedo = grey_albedo_lw, &
          l_invert = l_invert, &
          l_tile = .false., &
          water_soluble = water_soluble, &
          dust_like = dust_like, &
          oceanic = oceanic, &
          soot = soot, &
          ash = ash, &
          sulphuric = sulphuric, &
          ammonium_sulphate = ammonium_sulphate, &
          saharan_dust = saharan_dust, &
          accum_sulphate = accum_sulphate, &
          aitken_sulphate = aitken_sulphate, &
          fresh_soot = fresh_soot, &
          aged_soot = aged_soot, &
          sodium_chloride = sodium_chloride, &
          seasalt_film = seasalt_film, &
          seasalt_jet = seasalt_jet, &
          dust_div1 = dust_div1, &
          dust_div2 = dust_div2, &
          dust_div3 = dust_div3, &
          dust_div4 = dust_div4, &
          dust_div5 = dust_div5, &
          dust_div6 = dust_div6, &
          biomass_1 = biomass_1, &
          biomass_2 = biomass_2, &
          biogenic = biogenic, &
          ocff_fresh = ocff_fresh, &
          ocff_aged = ocff_aged, &
          delta = delta, &
          murk = murk, &
          nitrate = nitrate, &
          twobindust_1 = twobindust_1, &
          twobindust_2 = twobindust_2, &
          l_water_soluble = .true., &
          l_dust_like = .true., &
          l_oceanic = .true., &
          l_soot = .true., &
          l_ash = .true., &
          l_sulphuric = .true., &
          l_ammonium_sulphate = .true., &
          l_saharan_dust = .true., &
          l_accum_sulphate = .true., &
          l_aitken_sulphate = .true., &
          l_fresh_soot = .true., &
          l_aged_soot = .true., &
          l_sodium_chloride = .true., &
          l_seasalt_film = .true., &
          l_seasalt_jet = .true., &
          l_dust_div1 = .true., &
          l_dust_div2 = .true., &
          l_dust_div3 = .true., &
          l_dust_div4 = .true., &
          l_dust_div5 = .true., &
          l_dust_div6 = .true., &
          l_biomass_1 = .true., &
          l_biomass_2 = .true., &
          l_biogenic = .true., &
          l_ocff_fresh = .true., &
          l_ocff_aged = .true., &
          l_delta = .true., &
          l_murk = .true., &
          l_nitrate = .true., &
          l_twobindust_1 = .true., &
          l_twobindust_2 = .true., &
          l_profile_last = .false., &
          l_debug = .false., &
          i_profile_debug = 1)

      ELSEIF (mode .EQ. 2) THEN
      
        ! SW call
        sw_diag%heating_rate => sw_heating_rate
        sw_diag%flux_up => sw_flux_up
        sw_diag%flux_down => sw_flux_down
        call runes( &
          n_profile = n_profile, &
          n_layer = n_layer, &
          diag = sw_diag, &
          spectrum_name = 'sw', &
          i_source = ip_source_illuminate, &
          i_cloud_representation = cloud_rep, &
          i_overlap = ip_overlap_max_random, &
          p_layer = p_layer, &
          t_layer = t_layer, &
          mass = d_mass, &
          density = density, &
          layer_heat_capacity = layer_heat_capacity, &
          h2o = h2o, &
          co2 = co2, &
          o3 = o3, &
          n2o = n2o, &
          co = co, &
          ch4 = ch4, &
          o2 = o2, &
          no = no, &
          so2 = so2, &
          no2 = no2, &
          nh3 = nh3, &
          hno3 = hno3, &
          n2 = n2, &
          cfc11 = cfc11, &
          cfc12 = cfc12, &
          cfc113 = cfc113, &
          hcfc22 = hcfc22, &
          hfc125 = hfc125, &
          hfc134a = hfc134a, &
          cfc114 = cfc114, &
          tio = tio, &
          vo = vo, &
          h2 = h2, &
          he = he, &
          ocs = ocs, &
          na = na, &
          k = k, &
          feh = feh, &
          crh = crh, &
          li = li, &
          rb = rb, &
          cs = cs, &
          ph3 = ph3, &
          c2h2 = c2h2, &
          hcn = hcn, &
          h2s = h2s, &
          ar = ar, &
          o = o, &
          n = n, &
          no3 = no3, &
          n2o5 = n2o5, &
          hono = hono, &
          ho2no2 = ho2no2, &
          h2o2 = h2o2, &
          c2h6 = c2h6, &
          ch3 = ch3, &
          h2co = h2co, &
          ho2 = ho2, &
          hdo = hdo, &
          hcl = hcl, &
          hf = hf, &
          cosso = cosso, &
          tosso = tosso, &
          yosos = yosos, &
          ch3cho = ch3cho, &
          ch3ooh = ch3ooh, &
          ch3coch3 = ch3coch3, &
          ch3cocho = ch3cocho, &
          chocho = chocho, &
          c2h5cho = c2h5cho, &
          hoch2cho = hoch2cho, &
          c2h5coch3 = c2h5coch3, &
          mvk = mvk, &
          macr = macr, &
          pan = pan, &
          ch3ono2 = ch3ono2, &
          cloud_frac = cloud_frac, &
          liq_frac = liq_frac, &
          ice_frac = ice_frac, &
          liq_mmr = liq_mmr, &
          ice_mmr = ice_mmr, &
          liq_dim = liq_dim, &
          ice_dim = ice_dim, &
          i_st_water = 5, &
          i_st_ice = 8, &
          cos_zenith_angle = cz_angle, &
          solar_irrad = irrad, &
          l_grey_albedo = .true., & 
          grey_albedo = grey_albedo_sw, &
          water_soluble = water_soluble, &
          dust_like = dust_like, &
          oceanic = oceanic, &
          soot = soot, &
          ash = ash, &
          sulphuric = sulphuric, &
          ammonium_sulphate = ammonium_sulphate, &
          saharan_dust = saharan_dust, &
          accum_sulphate = accum_sulphate, &
          aitken_sulphate = aitken_sulphate, &
          fresh_soot = fresh_soot, &
          aged_soot = aged_soot, &
          sodium_chloride = sodium_chloride, &
          seasalt_film = seasalt_film, &
          seasalt_jet = seasalt_jet, &
          dust_div1 = dust_div1, &
          dust_div2 = dust_div2, &
          dust_div3 = dust_div3, &
          dust_div4 = dust_div4, &
          dust_div5 = dust_div5, &
          dust_div6 = dust_div6, &
          biomass_1 = biomass_1, &
          biomass_2 = biomass_2, &
          biogenic = biogenic, &
          ocff_fresh = ocff_fresh, &
          ocff_aged = ocff_aged, &
          delta = delta, &
          murk = murk, &
          nitrate = nitrate, &
          twobindust_1 = twobindust_1, &
          twobindust_2 = twobindust_2, &
          l_water_soluble = .true., &
          l_dust_like = .true., &
          l_oceanic = .true., &
          l_soot = .true., &
          l_ash = .true., &
          l_sulphuric = .true., &
          l_ammonium_sulphate = .true., &
          l_saharan_dust = .true., &
          l_accum_sulphate = .true., &
          l_aitken_sulphate = .true., &
          l_fresh_soot = .true., &
          l_aged_soot = .true., &
          l_sodium_chloride = .true., &
          l_seasalt_film = .true., &
          l_seasalt_jet = .true., &
          l_dust_div1 = .true., &
          l_dust_div2 = .true., &
          l_dust_div3 = .true., &
          l_dust_div4 = .true., &
          l_dust_div5 = .true., &
          l_dust_div6 = .true., &
          l_biomass_1 = .true., &
          l_biomass_2 = .true., &
          l_biogenic = .true., &
          l_ocff_fresh = .true., &
          l_ocff_aged = .true., &
          l_delta = .true., &
          l_murk = .true., &
          l_nitrate = .true., &
          l_twobindust_1 = .true., &
          l_twobindust_2 = .true., &
          l_rayleigh = .true., &
          l_invert = l_invert) ! If true, profiles can be supplied bottom-up

      ENDIF
          
  end subroutine runes_driver
  end module runes_driver_mod
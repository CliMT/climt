! module that contains the global variables

module globals

    use realtype_rd, only: RealExt

    implicit none

    real(RealExt), parameter :: pi = 4.0*atan(1.0)

    ! planetary parameters

    real(RealExt) :: g
    real(RealExt) :: r_gas_dry
    real(RealExt) :: cp_air_dry
    real(RealExt) :: seconds_per_day

    real(RealExt) :: zenith_angle
    real(RealExt) :: solar_irrad
    real(RealExt) :: grey_albedo_sw
    real(RealExt) :: grey_albedo_lw
    
end module globals

! module that contains the wrapper codes

module runes_driver_wrap

    use iso_c_binding, only: c_int, c_double, c_char
    use runes_driver_mod, only: runes_driver
    use socrates_set_spectrum, only: set_spectrum, spectrum_array, spectrum_array_name
    use def_spectrum, only: deallocate_spectrum
    use realtype_rd, only: RealExt
    use globals

    implicit none
    public

    contains

    ! spectrum subroutine

    subroutine input_spectrum_lw(spectral_file_c) bind (c)
        
        character(kind=c_char), intent(in) :: spectral_file_c(*)
        character(len=1024) :: spectral_file
        integer :: i
        integer :: k

        spectral_file = ""
        do i = 1, 1024
            if (spectral_file_c(i) == char(0)) exit
            spectral_file(i:i) = spectral_file_c(i)
        end do

        if (allocated(spectrum_array_name)) then
            do k = 1, size(spectrum_array_name)
                if (spectrum_array_name(k) == 'lw') then
                    call deallocate_spectrum(spectrum_array(k))
                    spectrum_array_name(k) = ''
                    exit
                end if
            end do
        end if

        call set_spectrum( &
            spectrum_name = 'lw', &
            spectral_file = trim(spectral_file), &
            l_all_gases = .true. )
        
    end subroutine input_spectrum_lw

    subroutine input_spectrum_sw(spectral_file_c) bind (c)
        
        character(kind=c_char), intent(in) :: spectral_file_c(*)
        character(len=1024) :: spectral_file
        integer :: i
        integer :: k

        spectral_file = ""
        do i = 1, 1024
            if (spectral_file_c(i) == char(0)) exit
            spectral_file(i:i) = spectral_file_c(i)
        end do

        if (allocated(spectrum_array_name)) then
            do k = 1, size(spectrum_array_name)
                if (spectrum_array_name(k) == 'sw') then
                    call deallocate_spectrum(spectrum_array(k))
                    spectrum_array_name(k) = ''
                    exit
                end if
            end do
        end if

        call set_spectrum( &
        spectrum_name = 'sw', &
        spectral_file = trim(spectral_file), &
        l_all_gases = .true. )
        
    end subroutine input_spectrum_sw

    subroutine set_planetary_parameters(g_c, seconds_per_day_c, r_gas_dry_c, cp_air_dry_c, &
        zenith_angle_c, solar_irrad_c, grey_albedo_lw_c, grey_albedo_sw_c) bind (c)
        
        real(c_double) :: g_c
        real(c_double) :: r_gas_dry_c
        real(c_double) :: cp_air_dry_c
        real(c_double) :: seconds_per_day_c

        real(c_double) :: zenith_angle_c
        real(c_double) :: solar_irrad_c
        real(c_double) :: grey_albedo_sw_c
        real(c_double) :: grey_albedo_lw_c

        ! changing scope
        g = g_c
        r_gas_dry = r_gas_dry_c
        cp_air_dry = cp_air_dry_c
        seconds_per_day = seconds_per_day_c

        zenith_angle = zenith_angle_c
        solar_irrad = solar_irrad_c
        grey_albedo_sw = grey_albedo_sw_c
        grey_albedo_lw = grey_albedo_lw_c

    end subroutine set_planetary_parameters

    subroutine set_fields_lw(&
        n_profile, &
        n_layer, &
        n_band, &
        h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
        cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
        he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
        n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
        hdo, hcl, hf, cosso, tosso, yosos, &
        ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
        hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2, &
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
        p_layer, &
        t_layer, &
        p_level, &
        t_level, &
        ground_temp, &
        flux_ground, &
        lw_flux_down, &
        lw_flux_up, &
        lw_heating_rate, &
        cloud_frac_c, liq_mmr_c, ice_mmr_c, liq_dim_c, ice_dim_c, &
        l_invert_c) bind (c)

        ! profile parameters

        integer(c_int), intent(in), value :: n_profile
        integer(c_int), intent(in), value :: n_layer
        integer(c_int), intent(in), value :: n_band
        integer(c_int), intent(in), value :: l_invert_c

        ! gas fields
        real(c_double), intent(in) :: h2o(n_profile, n_layer)
        real(c_double), intent(in) :: co2(n_profile, n_layer)
        real(c_double), intent(in) :: o3(n_profile, n_layer)
        real(c_double), intent(in) :: n2o(n_profile, n_layer)
        real(c_double), intent(in) :: co(n_profile, n_layer)
        real(c_double), intent(in) :: ch4(n_profile, n_layer)
        real(c_double), intent(in) :: o2(n_profile, n_layer)
        real(c_double), intent(in) :: no(n_profile, n_layer)
        real(c_double), intent(in) :: so2(n_profile, n_layer)
        real(c_double), intent(in) :: no2(n_profile, n_layer)
        real(c_double), intent(in) :: nh3(n_profile, n_layer)
        real(c_double), intent(in) :: hno3(n_profile, n_layer)
        real(c_double), intent(in) :: n2(n_profile, n_layer)
        real(c_double), intent(in) :: cfc11(n_profile, n_layer)
        real(c_double), intent(in) :: cfc12(n_profile, n_layer)
        real(c_double), intent(in) :: cfc113(n_profile, n_layer)
        real(c_double), intent(in) :: hcfc22(n_profile, n_layer)
        real(c_double), intent(in) :: hfc125(n_profile, n_layer)
        real(c_double), intent(in) :: hfc134a(n_profile, n_layer)
        real(c_double), intent(in) :: cfc114(n_profile, n_layer)
        real(c_double), intent(in) :: tio(n_profile, n_layer)
        real(c_double), intent(in) :: vo(n_profile, n_layer)
        real(c_double), intent(in) :: h2(n_profile, n_layer)
        real(c_double), intent(in) :: he(n_profile, n_layer)
        real(c_double), intent(in) :: ocs(n_profile, n_layer)
        real(c_double), intent(in) :: na(n_profile, n_layer)
        real(c_double), intent(in) :: k(n_profile, n_layer)
        real(c_double), intent(in) :: feh(n_profile, n_layer)
        real(c_double), intent(in) :: crh(n_profile, n_layer)
        real(c_double), intent(in) :: li(n_profile, n_layer)
        real(c_double), intent(in) :: rb(n_profile, n_layer)
        real(c_double), intent(in) :: cs(n_profile, n_layer)
        real(c_double), intent(in) :: ph3(n_profile, n_layer)
        real(c_double), intent(in) :: c2h2(n_profile, n_layer)
        real(c_double), intent(in) :: hcn(n_profile, n_layer)
        real(c_double), intent(in) :: h2s(n_profile, n_layer)
        real(c_double), intent(in) :: ar(n_profile, n_layer)
        real(c_double), intent(in) :: o(n_profile, n_layer)
        real(c_double), intent(in) :: n(n_profile, n_layer)
        real(c_double), intent(in) :: no3(n_profile, n_layer)
        real(c_double), intent(in) :: n2o5(n_profile, n_layer)
        real(c_double), intent(in) :: hono(n_profile, n_layer)
        real(c_double), intent(in) :: ho2no2(n_profile, n_layer)
        real(c_double), intent(in) :: h2o2(n_profile, n_layer)
        real(c_double), intent(in) :: c2h6(n_profile, n_layer)
        real(c_double), intent(in) :: ch3(n_profile, n_layer)
        real(c_double), intent(in) :: h2co(n_profile, n_layer)
        real(c_double), intent(in) :: ho2(n_profile, n_layer)
        real(c_double), intent(in) :: hdo(n_profile, n_layer)
        real(c_double), intent(in) :: hcl(n_profile, n_layer)
        real(c_double), intent(in) :: hf(n_profile, n_layer)
        real(c_double), intent(in) :: cosso(n_profile, n_layer)
        real(c_double), intent(in) :: tosso(n_profile, n_layer)
        real(c_double), intent(in) :: yosos(n_profile, n_layer)
        real(c_double), intent(in) :: ch3cho(n_profile, n_layer)
        real(c_double), intent(in) :: ch3ooh(n_profile, n_layer)
        real(c_double), intent(in) :: ch3coch3(n_profile, n_layer)
        real(c_double), intent(in) :: ch3cocho(n_profile, n_layer)
        real(c_double), intent(in) :: chocho(n_profile, n_layer)
        real(c_double), intent(in) :: c2h5cho(n_profile, n_layer)
        real(c_double), intent(in) :: hoch2cho(n_profile, n_layer)
        real(c_double), intent(in) :: c2h5coch3(n_profile, n_layer)
        real(c_double), intent(in) :: mvk(n_profile, n_layer)
        real(c_double), intent(in) :: macr(n_profile, n_layer)
        real(c_double), intent(in) :: pan(n_profile, n_layer)
        real(c_double), intent(in) :: ch3ono2(n_profile, n_layer)

        ! aerosols
        real(c_double), intent(in), dimension(n_profile, n_layer) :: water_soluble, &
        dust_like, &
        oceanic, &
        soot, &
        ash, &
        sulphuric, &
        ammonium_sulphate, &
        saharan_dust, &
        accum_sulphate, &
        aitken_sulphate, &
        fresh_soot, &
        aged_soot, &
        sodium_chloride, &
        seasalt_film, &
        seasalt_jet, &
        dust_div1, &
        dust_div2, &
        dust_div3, &
        dust_div4, &
        dust_div5, &
        dust_div6, &
        biomass_1, &
        biomass_2, &
        biogenic, &
        ocff_fresh, &
        ocff_aged, &
        delta, &
        murk, &
        nitrate, &
        twobindust_1, &
        twobindust_2

        ! pressure and temperature fields

        real(c_double), intent(in) :: p_layer(n_profile, n_layer)
        real(c_double), intent(in) :: t_layer(n_profile, n_layer)
        real(c_double), intent(in) :: p_level(n_profile, n_layer + 1)
        real(c_double), intent(in) :: t_level(n_profile, n_layer + 1)
        real(c_double), intent(in) :: flux_ground(n_profile, n_band)
        real(c_double), intent(in) :: ground_temp    

        ! clouds
        real(c_double), intent(in) :: cloud_frac_c(n_profile, n_layer)
        real(c_double), intent(in) :: liq_mmr_c(n_profile, n_layer)
        real(c_double), intent(in) :: ice_mmr_c(n_profile, n_layer)
        real(c_double), intent(in) :: liq_dim_c(n_profile, n_layer)
        real(c_double), intent(in) :: ice_dim_c(n_profile, n_layer)

        ! fluxes and heating rates

        real(c_double), intent(out) :: lw_flux_down(n_profile, n_layer + 1)
        real(c_double), intent(out) :: lw_flux_up(n_profile, n_layer + 1)
        real(c_double), intent(out) :: lw_heating_rate(n_profile, n_layer)

        integer :: prof_idx

        ! calling the loop

        do prof_idx = 1, n_profile
            call runes_driver( &
            g = g, &
            r_gas_dry = r_gas_dry, &
            cp_air_dry = cp_air_dry, &
            seconds_per_day = seconds_per_day, &
            h2o = h2o(prof_idx:prof_idx, :), &
            co2 = co2(prof_idx:prof_idx, :), &
            o3 = o3(prof_idx:prof_idx, :), &
            n2o = n2o(prof_idx:prof_idx, :), &
            co = co(prof_idx:prof_idx, :), &
            ch4 = ch4(prof_idx:prof_idx, :), &
            o2 = o2(prof_idx:prof_idx, :), &
            no = no(prof_idx:prof_idx, :), &
            so2 = so2(prof_idx:prof_idx, :), &
            no2 = no2(prof_idx:prof_idx, :), &
            nh3 = nh3(prof_idx:prof_idx, :), &
            hno3 = hno3(prof_idx:prof_idx, :), &
            n2 = n2(prof_idx:prof_idx, :), &
            cfc11 = cfc11(prof_idx:prof_idx, :), &
            cfc12 = cfc12(prof_idx:prof_idx, :), &
            cfc113 = cfc113(prof_idx:prof_idx, :), &
            hcfc22 = hcfc22(prof_idx:prof_idx, :), &
            hfc125 = hfc125(prof_idx:prof_idx, :), &
            hfc134a = hfc134a(prof_idx:prof_idx, :), &
            cfc114 = cfc114(prof_idx:prof_idx, :), &
            tio = tio(prof_idx:prof_idx, :), &
            vo = vo(prof_idx:prof_idx, :), &
            h2 = h2(prof_idx:prof_idx, :), &
            he = he(prof_idx:prof_idx, :), &
            ocs = ocs(prof_idx:prof_idx, :), &
            na = na(prof_idx:prof_idx, :), &
            k = k(prof_idx:prof_idx, :), &
            feh = feh(prof_idx:prof_idx, :), &
            crh = crh(prof_idx:prof_idx, :), &
            li = li(prof_idx:prof_idx, :), &
            rb = rb(prof_idx:prof_idx, :), &
            cs = cs(prof_idx:prof_idx, :), &
            ph3 = ph3(prof_idx:prof_idx, :), &
            c2h2 = c2h2(prof_idx:prof_idx, :), &
            hcn = hcn(prof_idx:prof_idx, :), &
            h2s = h2s(prof_idx:prof_idx, :), &
            ar = ar(prof_idx:prof_idx, :), &
            o = o(prof_idx:prof_idx, :), &
            n = n(prof_idx:prof_idx, :), &
            no3 = no3(prof_idx:prof_idx, :), &
            n2o5 = n2o5(prof_idx:prof_idx, :), &
            hono = hono(prof_idx:prof_idx, :), &
            ho2no2 = ho2no2(prof_idx:prof_idx, :), &
            h2o2 = h2o2(prof_idx:prof_idx, :), &
            c2h6 = c2h6(prof_idx:prof_idx, :), &
            ch3 = ch3(prof_idx:prof_idx, :), &
            h2co = h2co(prof_idx:prof_idx, :), &
            ho2 = ho2(prof_idx:prof_idx, :), &
            hdo = hdo(prof_idx:prof_idx, :), &
            hcl = hcl(prof_idx:prof_idx, :), &
            hf = hf(prof_idx:prof_idx, :), &
            cosso = cosso(prof_idx:prof_idx, :), &
            tosso = tosso(prof_idx:prof_idx, :), &
            yosos = yosos(prof_idx:prof_idx, :), &
            ch3cho = ch3cho(prof_idx:prof_idx, :), &
            ch3ooh = ch3ooh(prof_idx:prof_idx, :), &
            ch3coch3 = ch3coch3(prof_idx:prof_idx, :), &
            ch3cocho = ch3cocho(prof_idx:prof_idx, :), &
            chocho = chocho(prof_idx:prof_idx, :), &
            c2h5cho = c2h5cho(prof_idx:prof_idx, :), &
            hoch2cho = hoch2cho(prof_idx:prof_idx, :), &
            c2h5coch3 = c2h5coch3(prof_idx:prof_idx, :), &
            mvk = mvk(prof_idx:prof_idx, :), &
            macr = macr(prof_idx:prof_idx, :), &
            pan = pan(prof_idx:prof_idx, :), &
            ch3ono2 = ch3ono2(prof_idx:prof_idx, :), &
            zenith_angle = zenith_angle, &
            solar_irrad = solar_irrad, &
            grey_albedo_sw = grey_albedo_sw, &
            grey_albedo_lw = grey_albedo_lw, &
            ground_temp = ground_temp, &
            p_layer = p_layer(prof_idx:prof_idx, :), &
            t_layer = t_layer(prof_idx:prof_idx, :), &
            p_level = p_level(prof_idx:prof_idx, :), &
            t_level = t_level(prof_idx:prof_idx, :), &
            flux_ground = flux_ground(prof_idx:prof_idx, :), &
            n_layer = n_layer, &
            n_band = n_band, &
            lw_flux_down = lw_flux_down(prof_idx:prof_idx, :), &
            lw_flux_up = lw_flux_up(prof_idx:prof_idx, :), &
            lw_heating_rate = lw_heating_rate(prof_idx:prof_idx, :), &
            water_soluble = water_soluble(prof_idx:prof_idx, :), &
            dust_like = dust_like(prof_idx:prof_idx, :), &
            oceanic = oceanic(prof_idx:prof_idx, :), &
            soot = soot(prof_idx:prof_idx, :), &
            ash = ash(prof_idx:prof_idx, :), &
            sulphuric = sulphuric(prof_idx:prof_idx, :), &
            ammonium_sulphate = ammonium_sulphate(prof_idx:prof_idx, :), &
            saharan_dust = saharan_dust(prof_idx:prof_idx, :), &
            accum_sulphate = accum_sulphate(prof_idx:prof_idx, :), &
            aitken_sulphate = aitken_sulphate(prof_idx:prof_idx, :), &
            fresh_soot = fresh_soot(prof_idx:prof_idx, :), &
            aged_soot = aged_soot(prof_idx:prof_idx, :), &
            sodium_chloride = sodium_chloride(prof_idx:prof_idx, :), &
            seasalt_film = seasalt_film(prof_idx:prof_idx, :), &
            seasalt_jet = seasalt_jet(prof_idx:prof_idx, :), &
            dust_div1 = dust_div1(prof_idx:prof_idx, :), &
            dust_div2 = dust_div2(prof_idx:prof_idx, :), &
            dust_div3 = dust_div3(prof_idx:prof_idx, :), &
            dust_div4 = dust_div4(prof_idx:prof_idx, :), &
            dust_div5 = dust_div5(prof_idx:prof_idx, :), &
            dust_div6 = dust_div6(prof_idx:prof_idx, :), &
            biomass_1 = biomass_1(prof_idx:prof_idx, :), &
            biomass_2 = biomass_2(prof_idx:prof_idx, :), &
            biogenic = biogenic(prof_idx:prof_idx, :), &
            ocff_fresh = ocff_fresh(prof_idx:prof_idx, :), &
            ocff_aged = ocff_aged(prof_idx:prof_idx, :), &
            delta = delta(prof_idx:prof_idx, :), &
            murk = murk(prof_idx:prof_idx, :), &
            nitrate = nitrate(prof_idx:prof_idx, :), &
            twobindust_1 = twobindust_1(prof_idx:prof_idx, :), &
            twobindust_2 = twobindust_2(prof_idx:prof_idx, :), &
            cloud_frac = cloud_frac_c(prof_idx:prof_idx, :), &
            liq_mmr = liq_mmr_c(prof_idx:prof_idx, :), &
            ice_mmr = ice_mmr_c(prof_idx:prof_idx, :), &
            liq_dim = liq_dim_c(prof_idx:prof_idx, :), &
            ice_dim = ice_dim_c(prof_idx:prof_idx, :), &
            mode = 1, &
            l_invert = (l_invert_c .ne. 0))
        end do

    end subroutine set_fields_lw

    subroutine set_fields_sw(&
        n_profile, &
        n_layer, &
        n_band, &
        h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
        cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
        he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
        n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
        hdo, hcl, hf, cosso, tosso, yosos, &
        ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
        hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2, &
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
        p_layer, &
        t_layer, &
        p_level, &
        t_level, &
        ground_temp, &
        flux_ground, &
        sw_flux_down, &
        sw_flux_up, &
        sw_heating_rate, &
        cloud_frac_c, liq_mmr_c, ice_mmr_c, liq_dim_c, ice_dim_c, &
        l_invert_c) bind (c)

        ! profile parameters

        integer(c_int), intent(in), value :: n_profile
        integer(c_int), intent(in), value :: n_layer
        integer(c_int), intent(in), value :: n_band
        integer(c_int), intent(in), value :: l_invert_c

        ! gas fields
        real(c_double), intent(in) :: h2o(n_profile, n_layer)
        real(c_double), intent(in) :: co2(n_profile, n_layer)
        real(c_double), intent(in) :: o3(n_profile, n_layer)
        real(c_double), intent(in) :: n2o(n_profile, n_layer)
        real(c_double), intent(in) :: co(n_profile, n_layer)
        real(c_double), intent(in) :: ch4(n_profile, n_layer)
        real(c_double), intent(in) :: o2(n_profile, n_layer)
        real(c_double), intent(in) :: no(n_profile, n_layer)
        real(c_double), intent(in) :: so2(n_profile, n_layer)
        real(c_double), intent(in) :: no2(n_profile, n_layer)
        real(c_double), intent(in) :: nh3(n_profile, n_layer)
        real(c_double), intent(in) :: hno3(n_profile, n_layer)
        real(c_double), intent(in) :: n2(n_profile, n_layer)
        real(c_double), intent(in) :: cfc11(n_profile, n_layer)
        real(c_double), intent(in) :: cfc12(n_profile, n_layer)
        real(c_double), intent(in) :: cfc113(n_profile, n_layer)
        real(c_double), intent(in) :: hcfc22(n_profile, n_layer)
        real(c_double), intent(in) :: hfc125(n_profile, n_layer)
        real(c_double), intent(in) :: hfc134a(n_profile, n_layer)
        real(c_double), intent(in) :: cfc114(n_profile, n_layer)
        real(c_double), intent(in) :: tio(n_profile, n_layer)
        real(c_double), intent(in) :: vo(n_profile, n_layer)
        real(c_double), intent(in) :: h2(n_profile, n_layer)
        real(c_double), intent(in) :: he(n_profile, n_layer)
        real(c_double), intent(in) :: ocs(n_profile, n_layer)
        real(c_double), intent(in) :: na(n_profile, n_layer)
        real(c_double), intent(in) :: k(n_profile, n_layer)
        real(c_double), intent(in) :: feh(n_profile, n_layer)
        real(c_double), intent(in) :: crh(n_profile, n_layer)
        real(c_double), intent(in) :: li(n_profile, n_layer)
        real(c_double), intent(in) :: rb(n_profile, n_layer)
        real(c_double), intent(in) :: cs(n_profile, n_layer)
        real(c_double), intent(in) :: ph3(n_profile, n_layer)
        real(c_double), intent(in) :: c2h2(n_profile, n_layer)
        real(c_double), intent(in) :: hcn(n_profile, n_layer)
        real(c_double), intent(in) :: h2s(n_profile, n_layer)
        real(c_double), intent(in) :: ar(n_profile, n_layer)
        real(c_double), intent(in) :: o(n_profile, n_layer)
        real(c_double), intent(in) :: n(n_profile, n_layer)
        real(c_double), intent(in) :: no3(n_profile, n_layer)
        real(c_double), intent(in) :: n2o5(n_profile, n_layer)
        real(c_double), intent(in) :: hono(n_profile, n_layer)
        real(c_double), intent(in) :: ho2no2(n_profile, n_layer)
        real(c_double), intent(in) :: h2o2(n_profile, n_layer)
        real(c_double), intent(in) :: c2h6(n_profile, n_layer)
        real(c_double), intent(in) :: ch3(n_profile, n_layer)
        real(c_double), intent(in) :: h2co(n_profile, n_layer)
        real(c_double), intent(in) :: ho2(n_profile, n_layer)
        real(c_double), intent(in) :: hdo(n_profile, n_layer)
        real(c_double), intent(in) :: hcl(n_profile, n_layer)
        real(c_double), intent(in) :: hf(n_profile, n_layer)
        real(c_double), intent(in) :: cosso(n_profile, n_layer)
        real(c_double), intent(in) :: tosso(n_profile, n_layer)
        real(c_double), intent(in) :: yosos(n_profile, n_layer)
        real(c_double), intent(in) :: ch3cho(n_profile, n_layer)
        real(c_double), intent(in) :: ch3ooh(n_profile, n_layer)
        real(c_double), intent(in) :: ch3coch3(n_profile, n_layer)
        real(c_double), intent(in) :: ch3cocho(n_profile, n_layer)
        real(c_double), intent(in) :: chocho(n_profile, n_layer)
        real(c_double), intent(in) :: c2h5cho(n_profile, n_layer)
        real(c_double), intent(in) :: hoch2cho(n_profile, n_layer)
        real(c_double), intent(in) :: c2h5coch3(n_profile, n_layer)
        real(c_double), intent(in) :: mvk(n_profile, n_layer)
        real(c_double), intent(in) :: macr(n_profile, n_layer)
        real(c_double), intent(in) :: pan(n_profile, n_layer)
        real(c_double), intent(in) :: ch3ono2(n_profile, n_layer)

        ! aerosols
        real(c_double), intent(in), dimension(n_profile, n_layer) :: water_soluble, &
        dust_like, &
        oceanic, &
        soot, &
        ash, &
        sulphuric, &
        ammonium_sulphate, &
        saharan_dust, &
        accum_sulphate, &
        aitken_sulphate, &
        fresh_soot, &
        aged_soot, &
        sodium_chloride, &
        seasalt_film, &
        seasalt_jet, &
        dust_div1, &
        dust_div2, &
        dust_div3, &
        dust_div4, &
        dust_div5, &
        dust_div6, &
        biomass_1, &
        biomass_2, &
        biogenic, &
        ocff_fresh, &
        ocff_aged, &
        delta, &
        murk, &
        nitrate, &
        twobindust_1, &
        twobindust_2

        ! pressure and temperature fields

        real(c_double), intent(in) :: p_layer(n_profile, n_layer)
        real(c_double), intent(in) :: t_layer(n_profile, n_layer)
        real(c_double), intent(in) :: p_level(n_profile, n_layer + 1)
        real(c_double), intent(in) :: t_level(n_profile, n_layer + 1)
        real(c_double), intent(in) :: flux_ground(n_profile, n_band)
        real(c_double), intent(in) :: ground_temp    

        ! clouds
        real(c_double), intent(in) :: cloud_frac_c(n_profile, n_layer)
        real(c_double), intent(in) :: liq_mmr_c(n_profile, n_layer)
        real(c_double), intent(in) :: ice_mmr_c(n_profile, n_layer)
        real(c_double), intent(in) :: liq_dim_c(n_profile, n_layer)
        real(c_double), intent(in) :: ice_dim_c(n_profile, n_layer)

        ! fluxes and heating rates

        real(c_double), intent(out) :: sw_flux_down(n_profile, n_layer + 1)
        real(c_double), intent(out) :: sw_flux_up(n_profile, n_layer + 1)
        real(c_double), intent(out) :: sw_heating_rate(n_profile, n_layer)

        integer :: prof_idx

        ! calling the loop

        do prof_idx = 1, n_profile
            call runes_driver( &
            g = g, &
            r_gas_dry = r_gas_dry, &
            cp_air_dry = cp_air_dry, &
            seconds_per_day = seconds_per_day, &
            h2o = h2o(prof_idx:prof_idx, :), &
            co2 = co2(prof_idx:prof_idx, :), &
            o3 = o3(prof_idx:prof_idx, :), &
            n2o = n2o(prof_idx:prof_idx, :), &
            co = co(prof_idx:prof_idx, :), &
            ch4 = ch4(prof_idx:prof_idx, :), &
            o2 = o2(prof_idx:prof_idx, :), &
            no = no(prof_idx:prof_idx, :), &
            so2 = so2(prof_idx:prof_idx, :), &
            no2 = no2(prof_idx:prof_idx, :), &
            nh3 = nh3(prof_idx:prof_idx, :), &
            hno3 = hno3(prof_idx:prof_idx, :), &
            n2 = n2(prof_idx:prof_idx, :), &
            cfc11 = cfc11(prof_idx:prof_idx, :), &
            cfc12 = cfc12(prof_idx:prof_idx, :), &
            cfc113 = cfc113(prof_idx:prof_idx, :), &
            hcfc22 = hcfc22(prof_idx:prof_idx, :), &
            hfc125 = hfc125(prof_idx:prof_idx, :), &
            hfc134a = hfc134a(prof_idx:prof_idx, :), &
            cfc114 = cfc114(prof_idx:prof_idx, :), &
            tio = tio(prof_idx:prof_idx, :), &
            vo = vo(prof_idx:prof_idx, :), &
            h2 = h2(prof_idx:prof_idx, :), &
            he = he(prof_idx:prof_idx, :), &
            ocs = ocs(prof_idx:prof_idx, :), &
            na = na(prof_idx:prof_idx, :), &
            k = k(prof_idx:prof_idx, :), &
            feh = feh(prof_idx:prof_idx, :), &
            crh = crh(prof_idx:prof_idx, :), &
            li = li(prof_idx:prof_idx, :), &
            rb = rb(prof_idx:prof_idx, :), &
            cs = cs(prof_idx:prof_idx, :), &
            ph3 = ph3(prof_idx:prof_idx, :), &
            c2h2 = c2h2(prof_idx:prof_idx, :), &
            hcn = hcn(prof_idx:prof_idx, :), &
            h2s = h2s(prof_idx:prof_idx, :), &
            ar = ar(prof_idx:prof_idx, :), &
            o = o(prof_idx:prof_idx, :), &
            n = n(prof_idx:prof_idx, :), &
            no3 = no3(prof_idx:prof_idx, :), &
            n2o5 = n2o5(prof_idx:prof_idx, :), &
            hono = hono(prof_idx:prof_idx, :), &
            ho2no2 = ho2no2(prof_idx:prof_idx, :), &
            h2o2 = h2o2(prof_idx:prof_idx, :), &
            c2h6 = c2h6(prof_idx:prof_idx, :), &
            ch3 = ch3(prof_idx:prof_idx, :), &
            h2co = h2co(prof_idx:prof_idx, :), &
            ho2 = ho2(prof_idx:prof_idx, :), &
            hdo = hdo(prof_idx:prof_idx, :), &
            hcl = hcl(prof_idx:prof_idx, :), &
            hf = hf(prof_idx:prof_idx, :), &
            cosso = cosso(prof_idx:prof_idx, :), &
            tosso = tosso(prof_idx:prof_idx, :), &
            yosos = yosos(prof_idx:prof_idx, :), &
            ch3cho = ch3cho(prof_idx:prof_idx, :), &
            ch3ooh = ch3ooh(prof_idx:prof_idx, :), &
            ch3coch3 = ch3coch3(prof_idx:prof_idx, :), &
            ch3cocho = ch3cocho(prof_idx:prof_idx, :), &
            chocho = chocho(prof_idx:prof_idx, :), &
            c2h5cho = c2h5cho(prof_idx:prof_idx, :), &
            hoch2cho = hoch2cho(prof_idx:prof_idx, :), &
            c2h5coch3 = c2h5coch3(prof_idx:prof_idx, :), &
            mvk = mvk(prof_idx:prof_idx, :), &
            macr = macr(prof_idx:prof_idx, :), &
            pan = pan(prof_idx:prof_idx, :), &
            ch3ono2 = ch3ono2(prof_idx:prof_idx, :), &
            zenith_angle = zenith_angle, &
            solar_irrad = solar_irrad, &
            grey_albedo_sw = grey_albedo_sw, &
            grey_albedo_lw = grey_albedo_lw, &
            ground_temp = ground_temp, &
            p_layer = p_layer(prof_idx:prof_idx, :), &
            t_layer = t_layer(prof_idx:prof_idx, :), &
            p_level = p_level(prof_idx:prof_idx, :), &
            t_level = t_level(prof_idx:prof_idx, :), &
            flux_ground = flux_ground(prof_idx:prof_idx, :), &
            n_layer = n_layer, &
            n_band = n_band, &
            sw_flux_down = sw_flux_down(prof_idx:prof_idx, :), &
            sw_flux_up = sw_flux_up(prof_idx:prof_idx, :), &
            sw_heating_rate = sw_heating_rate(prof_idx:prof_idx, :), &
            water_soluble = water_soluble(prof_idx:prof_idx, :), &
            dust_like = dust_like(prof_idx:prof_idx, :), &
            oceanic = oceanic(prof_idx:prof_idx, :), &
            soot = soot(prof_idx:prof_idx, :), &
            ash = ash(prof_idx:prof_idx, :), &
            sulphuric = sulphuric(prof_idx:prof_idx, :), &
            ammonium_sulphate = ammonium_sulphate(prof_idx:prof_idx, :), &
            saharan_dust = saharan_dust(prof_idx:prof_idx, :), &
            accum_sulphate = accum_sulphate(prof_idx:prof_idx, :), &
            aitken_sulphate = aitken_sulphate(prof_idx:prof_idx, :), &
            fresh_soot = fresh_soot(prof_idx:prof_idx, :), &
            aged_soot = aged_soot(prof_idx:prof_idx, :), &
            sodium_chloride = sodium_chloride(prof_idx:prof_idx, :), &
            seasalt_film = seasalt_film(prof_idx:prof_idx, :), &
            seasalt_jet = seasalt_jet(prof_idx:prof_idx, :), &
            dust_div1 = dust_div1(prof_idx:prof_idx, :), &
            dust_div2 = dust_div2(prof_idx:prof_idx, :), &
            dust_div3 = dust_div3(prof_idx:prof_idx, :), &
            dust_div4 = dust_div4(prof_idx:prof_idx, :), &
            dust_div5 = dust_div5(prof_idx:prof_idx, :), &
            dust_div6 = dust_div6(prof_idx:prof_idx, :), &
            biomass_1 = biomass_1(prof_idx:prof_idx, :), &
            biomass_2 = biomass_2(prof_idx:prof_idx, :), &
            biogenic = biogenic(prof_idx:prof_idx, :), &
            ocff_fresh = ocff_fresh(prof_idx:prof_idx, :), &
            ocff_aged = ocff_aged(prof_idx:prof_idx, :), &
            delta = delta(prof_idx:prof_idx, :), &
            murk = murk(prof_idx:prof_idx, :), &
            nitrate = nitrate(prof_idx:prof_idx, :), &
            twobindust_1 = twobindust_1(prof_idx:prof_idx, :), &
            twobindust_2 = twobindust_2(prof_idx:prof_idx, :), &
            cloud_frac = cloud_frac_c(prof_idx:prof_idx, :), &
            liq_mmr = liq_mmr_c(prof_idx:prof_idx, :), &
            ice_mmr = ice_mmr_c(prof_idx:prof_idx, :), &
            liq_dim = liq_dim_c(prof_idx:prof_idx, :), &
            ice_dim = ice_dim_c(prof_idx:prof_idx, :), &
            mode = 2, &
            l_invert = (l_invert_c .ne. 0))
        end do

    end subroutine set_fields_sw

end module runes_driver_wrap
from ._runes_driver_mod import lib as _runes_driver_lib
from cffi import FFI
import numpy as np
from .utils import CDataConverter

ffi = FFI()

class _runesfile(object):
    def _input_spectrum_lw(self, spectral_file):
        _runes_driver_lib.input_spectrum_lw((spectral_file + '\0').encode('utf-8'))

    def _input_spectrum_sw(self, spectral_file):
        _runes_driver_lib.input_spectrum_sw((spectral_file + '\0').encode('utf-8'))

    def _set_planetary_parameters(self, g, seconds_per_day, r_gas_dry, cp_air_dry, zenith_angle, solar_irrad, grey_albedo_lw, grey_albedo_sw):
        to_c = CDataConverter(ffi)

        _runes_driver_lib.set_planetary_parameters(
            to_c(g), to_c(seconds_per_day), to_c(r_gas_dry), to_c(cp_air_dry),
            to_c(zenith_angle), to_c(solar_irrad), to_c(grey_albedo_lw), to_c(grey_albedo_sw)
        )

    def _set_fields_lw(self, n_profile, n_layer, n_band,
                       h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2,
                       cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2,
                       he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3,
                       n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2,
                       hdo, hcl, hf, cosso, tosso, yosos,
                       ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho,
                       hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2,
                       water_soluble, dust_like, oceanic, soot, ash, sulphuric,
                       ammonium_sulphate, saharan_dust,
                       accum_sulphate, aitken_sulphate,
                       fresh_soot, aged_soot,
                       sodium_chloride, seasalt_film, seasalt_jet,
                       dust_div1, dust_div2, dust_div3,
                       dust_div4, dust_div5, dust_div6,
                       biomass_1, biomass_2,
                       biogenic,
                       ocff_fresh, ocff_aged,
                       delta, murk,
                       nitrate,
                       twobindust_1, twobindust_2,
                       p_layer, t_layer, p_level, t_level,
                       ground_temp, flux_ground,
                       lw_flux_down, lw_flux_up, lw_heating_rate,
                       cloud_frac, liq_mmr, ice_mmr, liq_dim, ice_dim, l_invert):
        to_c = CDataConverter(ffi)

        _runes_driver_lib.set_fields_lw(
            n_profile, n_layer, n_band,
            to_c(h2o), to_c(co2), to_c(o3), to_c(n2o), to_c(co), to_c(ch4), to_c(o2), to_c(no), to_c(so2), to_c(no2), to_c(nh3), to_c(hno3), to_c(n2),
            to_c(cfc11), to_c(cfc12), to_c(cfc113), to_c(hcfc22), to_c(hfc125), to_c(hfc134a), to_c(cfc114), to_c(tio), to_c(vo), to_c(h2),
            to_c(he), to_c(ocs), to_c(na), to_c(k), to_c(feh), to_c(crh), to_c(li), to_c(rb), to_c(cs), to_c(ph3), to_c(c2h2), to_c(hcn), to_c(h2s), to_c(ar), to_c(o), to_c(n), to_c(no3),
            to_c(n2o5), to_c(hono), to_c(ho2no2), to_c(h2o2), to_c(c2h6), to_c(ch3), to_c(h2co), to_c(ho2),
            to_c(hdo), to_c(hcl), to_c(hf), to_c(cosso), to_c(tosso), to_c(yosos),
            to_c(ch3cho), to_c(ch3ooh), to_c(ch3coch3), to_c(ch3cocho), to_c(chocho), to_c(c2h5cho),
            to_c(hoch2cho), to_c(c2h5coch3), to_c(mvk), to_c(macr), to_c(pan), to_c(ch3ono2),
            to_c(water_soluble), to_c(dust_like), to_c(oceanic), to_c(soot), to_c(ash), to_c(sulphuric),
            to_c(ammonium_sulphate), to_c(saharan_dust),
            to_c(accum_sulphate), to_c(aitken_sulphate),
            to_c(fresh_soot), to_c(aged_soot),
            to_c(sodium_chloride), to_c(seasalt_film), to_c(seasalt_jet),
            to_c(dust_div1), to_c(dust_div2), to_c(dust_div3),
            to_c(dust_div4), to_c(dust_div5), to_c(dust_div6),
            to_c(biomass_1), to_c(biomass_2),
            to_c(biogenic),
            to_c(ocff_fresh), to_c(ocff_aged),
            to_c(delta), to_c(murk),
            to_c(nitrate),
            to_c(twobindust_1), to_c(twobindust_2),
            to_c(p_layer), to_c(t_layer), to_c(p_level), to_c(t_level),
            to_c(ground_temp), to_c(flux_ground),
            to_c(lw_flux_down), to_c(lw_flux_up), to_c(lw_heating_rate),
            to_c(cloud_frac), to_c(liq_mmr), to_c(ice_mmr), to_c(liq_dim), to_c(ice_dim),
            int(l_invert)
        )

    def _set_fields_sw(self, n_profile, n_layer, n_band,
                       h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2,
                       cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2,
                       he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3,
                       n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2,
                       hdo, hcl, hf, cosso, tosso, yosos,
                       ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho,
                       hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2,
                       water_soluble, dust_like, oceanic, soot, ash, sulphuric,
                       ammonium_sulphate, saharan_dust,
                       accum_sulphate, aitken_sulphate,
                       fresh_soot, aged_soot,
                       sodium_chloride, seasalt_film, seasalt_jet,
                       dust_div1, dust_div2, dust_div3,
                       dust_div4, dust_div5, dust_div6,
                       biomass_1, biomass_2,
                       biogenic,
                       ocff_fresh, ocff_aged,
                       delta, murk,
                       nitrate,
                       twobindust_1, twobindust_2,
                       p_layer, t_layer, p_level, t_level,
                       ground_temp, flux_ground,
                       sw_flux_down, sw_flux_up, sw_heating_rate,
                       cloud_frac, liq_mmr, ice_mmr, liq_dim, ice_dim, l_invert):
        to_c = CDataConverter(ffi)

        _runes_driver_lib.set_fields_sw(
            n_profile, n_layer, n_band,
            to_c(h2o), to_c(co2), to_c(o3), to_c(n2o), to_c(co), to_c(ch4), to_c(o2), to_c(no), to_c(so2), to_c(no2), to_c(nh3), to_c(hno3), to_c(n2),
            to_c(cfc11), to_c(cfc12), to_c(cfc113), to_c(hcfc22), to_c(hfc125), to_c(hfc134a), to_c(cfc114), to_c(tio), to_c(vo), to_c(h2),
            to_c(he), to_c(ocs), to_c(na), to_c(k), to_c(feh), to_c(crh), to_c(li), to_c(rb), to_c(cs), to_c(ph3), to_c(c2h2), to_c(hcn), to_c(h2s), to_c(ar), to_c(o), to_c(n), to_c(no3),
            to_c(n2o5), to_c(hono), to_c(ho2no2), to_c(h2o2), to_c(c2h6), to_c(ch3), to_c(h2co), to_c(ho2),
            to_c(hdo), to_c(hcl), to_c(hf), to_c(cosso), to_c(tosso), to_c(yosos),
            to_c(ch3cho), to_c(ch3ooh), to_c(ch3coch3), to_c(ch3cocho), to_c(chocho), to_c(c2h5cho),
            to_c(hoch2cho), to_c(c2h5coch3), to_c(mvk), to_c(macr), to_c(pan), to_c(ch3ono2),
            to_c(water_soluble), to_c(dust_like), to_c(oceanic), to_c(soot), to_c(ash), to_c(sulphuric),
            to_c(ammonium_sulphate), to_c(saharan_dust),
            to_c(accum_sulphate), to_c(aitken_sulphate),
            to_c(fresh_soot), to_c(aged_soot),
            to_c(sodium_chloride), to_c(seasalt_film), to_c(seasalt_jet),
            to_c(dust_div1), to_c(dust_div2), to_c(dust_div3),
            to_c(dust_div4), to_c(dust_div5), to_c(dust_div6),
            to_c(biomass_1), to_c(biomass_2),
            to_c(biogenic),
            to_c(ocff_fresh), to_c(ocff_aged),
            to_c(delta), to_c(murk),
            to_c(nitrate),
            to_c(twobindust_1), to_c(twobindust_2),
            to_c(p_layer), to_c(t_layer), to_c(p_level), to_c(t_level),
            to_c(ground_temp), to_c(flux_ground),
            to_c(sw_flux_down), to_c(sw_flux_up), to_c(sw_heating_rate),
            to_c(cloud_frac), to_c(liq_mmr), to_c(ice_mmr), to_c(liq_dim), to_c(ice_dim),
            int(l_invert)
        )

import os
import platform
from cffi import FFI

operating_system = platform.system()
dir_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # climt/climt/
lib_dir = os.path.join(dir_path, "_lib") # climt/climt/_lib
socrates_dir = os.path.join(lib_dir, "socrates") # climt/climt/_lib/socrates
platform_lib_dir = os.path.join(lib_dir, operating_system) # climt/climt/_lib/Linux/

ffi_builder = FFI()

ffi_builder.cdef(
    "void input_spectrum_lw(char *spectral_file);"
)
ffi_builder.cdef(
    "void input_spectrum_sw(char *spectral_file);"
)
ffi_builder.cdef(
    "void set_planetary_parameters(double *g_c, double *seconds_per_day_c, double *r_gas_dry_c, double *cp_air_dry_c, double *zenith_angle_c, double *solar_irrad_c, double *grey_albedo_lw_c, double *grey_albedo_sw_c);"
)
ffi_builder.cdef("""
void set_fields_lw(int n_profile, int n_layer, int n_band,
    double *h2o, double *co2, double *o3, double *n2o, double *co, double *ch4, double *o2, double *no, double *so2, double *no2, double *nh3, double *hno3, double *n2,
    double *cfc11, double *cfc12, double *cfc113, double *hcfc22, double *hfc125, double *hfc134a, double *cfc114, double *tio, double *vo, double *h2,
    double *he, double *ocs, double *na, double *k, double *feh, double *crh, double *li, double *rb, double *cs, double *ph3, double *c2h2, double *hcn, double *h2s, double *ar, double *o, double *n, double *no3,
    double *n2o5, double *hono, double *ho2no2, double *h2o2, double *c2h6, double *ch3, double *h2co, double *ho2,
    double *hdo, double *hcl, double *hf, double *cosso, double *tosso, double *yosos,
    double *ch3cho, double *ch3ooh, double *ch3coch3, double *ch3cocho, double *chocho, double *c2h5cho,
    double *hoch2cho, double *c2h5coch3, double *mvk, double *macr, double *pan, double *ch3ono2,
    double *water_soluble, double *dust_like, double *oceanic, double *soot, double *ash, double *sulphuric,
    double *ammonium_sulphate, double *saharan_dust,
    double *accum_sulphate, double *aitken_sulphate,
    double *fresh_soot, double *aged_soot,
    double *sodium_chloride, double *seasalt_film, double *seasalt_jet,
    double *dust_div1, double *dust_div2, double *dust_div3,
    double *dust_div4, double *dust_div5, double *dust_div6,
    double *biomass_1, double *biomass_2,
    double *biogenic,
    double *ocff_fresh, double *ocff_aged,
    double *delta, double *murk,
    double *nitrate,
    double *twobindust_1, double *twobindust_2,
    double *p_layer, double *t_layer, double *p_level, double *t_level,
    double *ground_temp, double *flux_ground,
    double *lw_flux_down, double *lw_flux_up, double *lw_heating_rate,
    double *cloud_frac_c, double *liq_mmr_c, double *ice_mmr_c, double *liq_dim_c, double *ice_dim_c, int l_invert_c);
""")
ffi_builder.cdef("""
void set_fields_sw(int n_profile, int n_layer, int n_band,
    double *h2o, double *co2, double *o3, double *n2o, double *co, double *ch4, double *o2, double *no, double *so2, double *no2, double *nh3, double *hno3, double *n2,
    double *cfc11, double *cfc12, double *cfc113, double *hcfc22, double *hfc125, double *hfc134a, double *cfc114, double *tio, double *vo, double *h2,
    double *he, double *ocs, double *na, double *k, double *feh, double *crh, double *li, double *rb, double *cs, double *ph3, double *c2h2, double *hcn, double *h2s, double *ar, double *o, double *n, double *no3,
    double *n2o5, double *hono, double *ho2no2, double *h2o2, double *c2h6, double *ch3, double *h2co, double *ho2,
    double *hdo, double *hcl, double *hf, double *cosso, double *tosso, double *yosos,
    double *ch3cho, double *ch3ooh, double *ch3coch3, double *ch3cocho, double *chocho, double *c2h5cho,
    double *hoch2cho, double *c2h5coch3, double *mvk, double *macr, double *pan, double *ch3ono2,
    double *water_soluble, double *dust_like, double *oceanic, double *soot, double *ash, double *sulphuric,
    double *ammonium_sulphate, double *saharan_dust,
    double *accum_sulphate, double *aitken_sulphate,
    double *fresh_soot, double *aged_soot,
    double *sodium_chloride, double *seasalt_film, double *seasalt_jet,
    double *dust_div1, double *dust_div2, double *dust_div3,
    double *dust_div4, double *dust_div5, double *dust_div6,
    double *biomass_1, double *biomass_2,
    double *biogenic,
    double *ocff_fresh, double *ocff_aged,
    double *delta, double *murk,
    double *nitrate,
    double *twobindust_1, double *twobindust_2,
    double *p_layer, double *t_layer, double *p_level, double *t_level,
    double *ground_temp, double *flux_ground,
    double *sw_flux_down, double *sw_flux_up, double *sw_heating_rate,
    double *cloud_frac_c, double *liq_mmr_c, double *ice_mmr_c, double *liq_dim_c, double *ice_dim_c, int l_invert_c);
""")

extra_link_args = ["-lgfortran"]
if operating_system == "Linux":
    extra_link_args.extend(["-lm"])

ffi_builder.set_source(
    "climt._components.socrates._runes_driver_mod",
    '#include "runeswrap.h"',
    library_dirs=[platform_lib_dir],
    include_dirs=[os.path.dirname(os.path.abspath(__file__))],
    libraries=["runes"],
    extra_link_args=extra_link_args,
)

if __name__ == "__main__":
    ffi_builder.compile(verbose=True)

import os
import numpy as np


def resolve_spectral_file(filename):
    if os.path.isabs(filename):
        return filename

    # Resolve relative to packaged socrates data directory
    socrates_data_dir = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "_lib",
            "socrates",
            "socrates_2511",
            "data"
        )
    )

    # Walk the directory to find the file
    for root, dirs, files in os.walk(socrates_data_dir):
        if filename in files:
            return os.path.join(root, filename)

    # Return original if not found
    return filename


# Full mapping of Socrates gas parameters to standard CliMT/sympl names
gas_name_mapping = {
    'h2o': 'specific_humidity',
    'co2': 'mole_fraction_of_carbon_dioxide_in_air',
    'o3': 'mole_fraction_of_ozone_in_air',
    'n2o': 'mole_fraction_of_nitrous_oxide_in_air',
    'co': 'mole_fraction_of_carbon_monoxide_in_air',
    'ch4': 'mole_fraction_of_methane_in_air',
    'o2': 'mole_fraction_of_oxygen_in_air',
    'no': 'mole_fraction_of_nitric_oxide_in_air',
    'so2': 'mole_fraction_of_sulphur_dioxide_in_air',
    'no2': 'mole_fraction_of_nitrogen_dioxide_in_air',
    'nh3': 'mole_fraction_of_ammonia_in_air',
    'hno3': 'mole_fraction_of_nitric_acid_in_air',
    'n2': 'mole_fraction_of_nitrogen_in_air',
    'cfc11': 'mole_fraction_of_cfc11_in_air',
    'cfc12': 'mole_fraction_of_cfc12_in_air',
    'cfc113': 'mole_fraction_of_cfc113_in_air',
    'hcfc22': 'mole_fraction_of_hcfc22_in_air',
    'hfc125': 'mole_fraction_of_hfc125_in_air',
    'hfc134a': 'mole_fraction_of_hfc134a_in_air',
    'cfc114': 'mole_fraction_of_cfc114_in_air',
    'tio': 'mole_fraction_of_tio_in_air',
    'vo': 'mole_fraction_of_vo_in_air',
    'h2': 'mole_fraction_of_hydrogen_in_air',
    'he': 'mole_fraction_of_helium_in_air',
    'ocs': 'mole_fraction_of_carbonyl_sulfide_in_air',
    'na': 'mole_fraction_of_sodium_in_air',
    'k': 'mole_fraction_of_potassium_in_air',
    'feh': 'mole_fraction_of_feh_in_air',
    'crh': 'mole_fraction_of_crh_in_air',
    'li': 'mole_fraction_of_lithium_in_air',
    'rb': 'mole_fraction_of_rubidium_in_air',
    'cs': 'mole_fraction_of_caesium_in_air',
    'ph3': 'mole_fraction_of_phosphine_in_air',
    'c2h2': 'mole_fraction_of_acetylene_in_air',
    'hcn': 'mole_fraction_of_hydrogen_cyanide_in_air',
    'h2s': 'mole_fraction_of_hydrogen_sulfide_in_air',
    'ar': 'mole_fraction_of_argon_in_air',
    'o': 'mole_fraction_of_atomic_oxygen_in_air',
    'n': 'mole_fraction_of_atomic_nitrogen_in_air',
    'no3': 'mole_fraction_of_nitrate_radical_in_air',
    'n2o5': 'mole_fraction_of_dinitrogen_pentoxide_in_air',
    'hono': 'mole_fraction_of_nitrous_acid_in_air',
    'ho2no2': 'mole_fraction_of_peroxynitric_acid_in_air',
    'h2o2': 'mole_fraction_of_hydrogen_peroxide_in_air',
    'c2h6': 'mole_fraction_of_ethane_in_air',
    'ch3': 'mole_fraction_of_methyl_radical_in_air',
    'h2co': 'mole_fraction_of_formaldehyde_in_air',
    'ho2': 'mole_fraction_of_hydroperoxyl_radical_in_air',
    'hdo': 'mole_fraction_of_deuterated_water_in_air',
    'hcl': 'mole_fraction_of_hydrogen_chloride_in_air',
    'hf': 'mole_fraction_of_hydrogen_fluoride_in_air',
    'cosso': 'mole_fraction_of_cosso_in_air',
    'tosso': 'mole_fraction_of_tosso_in_air',
    'yosos': 'mole_fraction_of_yosos_in_air',
    'ch3cho': 'mole_fraction_of_acetaldehyde_in_air',
    'ch3ooh': 'mole_fraction_of_methyl_hydroperoxide_in_air',
    'ch3coch3': 'mole_fraction_of_acetone_in_air',
    'ch3cocho': 'mole_fraction_of_methylglyoxal_in_air',
    'chocho': 'mole_fraction_of_glyoxal_in_air',
    'c2h5cho': 'mole_fraction_of_propanal_in_air',
    'hoch2cho': 'mole_fraction_of_glycolaldehyde_in_air',
    'c2h5coch3': 'mole_fraction_of_methyl_ethyl_ketone_in_air',
    'mvk': 'mole_fraction_of_methyl_vinyl_ketone_in_air',
    'macr': 'mole_fraction_of_methacrolein_in_air',
    'pan': 'mole_fraction_of_peroxyacetyl_nitrate_in_air',
    'ch3ono2': 'mole_fraction_of_methyl_nitrate_in_air',
}

# Full mapping of Socrates aerosol parameters to standard CliMT/sympl names
aerosol_name_mapping = {
    'water_soluble': 'mole_fraction_of_water_soluble_in_air',
    'dust_like': 'mole_fraction_of_dust_like_in_air',
    'oceanic': 'mole_fraction_of_oceanic_in_air',
    'soot': 'mole_fraction_of_soot_in_air',
    'ash': 'mole_fraction_of_ash_in_air',
    'sulphuric': 'mole_fraction_of_sulphuric_in_air',
    'ammonium_sulphate': 'mole_fraction_of_ammonium_sulphate_in_air',
    'saharan_dust': 'mole_fraction_of_saharan_dust_in_air',
    'accum_sulphate': 'mole_fraction_of_accum_sulphate_in_air',
    'aitken_sulphate': 'mole_fraction_of_aitken_sulphate_in_air',
    'fresh_soot': 'mole_fraction_of_fresh_soot_in_air',
    'aged_soot': 'mole_fraction_of_aged_soot_in_air',
    'sodium_chloride': 'mole_fraction_of_sodium_chloride_in_air',
    'seasalt_film': 'mole_fraction_of_seasalt_film_in_air',
    'seasalt_jet': 'mole_fraction_of_seasalt_jet_in_air',
    'dust_div1': 'mole_fraction_of_dust_div1_in_air',
    'dust_div2': 'mole_fraction_of_dust_div2_in_air',
    'dust_div3': 'mole_fraction_of_dust_div3_in_air',
    'dust_div4': 'mole_fraction_of_dust_div4_in_air',
    'dust_div5': 'mole_fraction_of_dust_div5_in_air',
    'dust_div6': 'mole_fraction_of_dust_div6_in_air',
    'biomass_1': 'mole_fraction_of_biomass_1_in_air',
    'biomass_2': 'mole_fraction_of_biomass_2_in_air',
    'biogenic': 'mole_fraction_of_biogenic_in_air',
    'ocff_fresh': 'mole_fraction_of_ocff_fresh_in_air',
    'ocff_aged': 'mole_fraction_of_ocff_aged_in_air',
    'delta': 'mole_fraction_of_delta_in_air',
    'murk': 'mole_fraction_of_murk_in_air',
    'nitrate': 'mole_fraction_of_nitrate_in_air',
    'twobindust_1': 'mole_fraction_of_twobindust_1_in_air',
    'twobindust_2': 'mole_fraction_of_twobindust_2_in_air',
}

# Build static input properties dictionary dynamically for module load
base_input_properties = {
    'air_pressure': {
        'dims': ['mid_levels', '*'],
        'units': 'Pa',
    },
    'air_pressure_on_interface_levels': {
        'dims': ['interface_levels', '*'],
        'units': 'Pa',
    },
    'air_temperature': {
        'dims': ['mid_levels', '*'],
        'units': 'degK',
    },
    'air_temperature_on_interface_levels': {
        'dims': ['interface_levels', '*'],
        'units': 'degK',
    },
    'surface_temperature': {
        'dims': ['*'],
        'units': 'degK',
    },
    'zenith_angle': {
        'dims': ['*'],
        'units': 'radians'
    },
    # Cloud fields (optional)
    'cloud_area_fraction_in_atmosphere_layer': {
        'dims': ['mid_levels', '*'],
        'units': 'dimensionless',
    },
    'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {
        'dims': ['mid_levels', '*'],
        'units': 'kg m^-2',
    },
    'mass_content_of_cloud_ice_in_atmosphere_layer': {
        'dims': ['mid_levels', '*'],
        'units': 'kg m^-2',
    },
    'cloud_water_droplet_radius': {
        'dims': ['mid_levels', '*'],
        'units': 'micrometer',
    },
    'cloud_ice_particle_size': {
        'dims': ['mid_levels', '*'],
        'units': 'micrometer',
    },
}

for gas_key, prop_name in gas_name_mapping.items():
    if prop_name != 'specific_humidity':
        base_input_properties[prop_name] = {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        }
    else:
        base_input_properties[prop_name] = {
            'dims': ['mid_levels', '*'],
            'units': 'kg/kg',
        }

for aerosol_key, prop_name in aerosol_name_mapping.items():
    base_input_properties[prop_name] = {
        'dims': ['mid_levels', '*'],
        'units': 'dimensionless',
    }


def get_num_bands(file_path):
    try:
        with open(file_path, 'r', errors='ignore') as f:
            for _ in range(20):
                line = f.readline()
                if not line:
                    break
                if "Number of spectral bands" in line:
                    parts = line.split('=')
                    if len(parts) > 1:
                        return int(parts[1].strip())
    except Exception:
        pass
    return None


class CDataConverter:
    def __init__(self, ffi):
        self.ffi = ffi
        self.temp_arrays = []

    def __call__(self, arr):
        arr_contiguous = np.ascontiguousarray(arr, dtype=np.float64)
        self.temp_arrays.append(arr_contiguous)
        return self.ffi.cast('double *', arr_contiguous.ctypes.data)

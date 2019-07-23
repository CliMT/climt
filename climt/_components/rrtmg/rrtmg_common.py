"""
Dictionary mapping input options of RRTMG radiative
components to integers accepted by the RRTMG fortran
code
"""
rrtmg_cloud_overlap_method_dict = {
    'clear_only': 0,
    'random': 1,
    'maximum_random': 2,
    'maximum': 3
}

"""
Dictionary mapping input options of RRTMG radiative
components to integers accepted by the RRTMG fortran
code
"""
rrtmg_cloud_props_dict = {
    'direct_input': 0,
    'single_cloud_type': 1,
    'liquid_and_ice_clouds': 2
}

"""
Dictionary mapping input options of RRTMG radiative
components to integers accepted by the RRTMG fortran
code
"""
rrtmg_cloud_ice_props_dict = {
    'ebert_curry_one': 0,
    'ebert_curry_two': 1,
    'key_streamer_manual': 2,
    'fu': 3
}

"""
Dictionary mapping input options of RRTMG radiative
components to integers accepted by the RRTMG fortran
code
"""
rrtmg_cloud_liquid_props_dict = {
    'radius_independent_absorption': 0,
    'radius_dependent_absorption': 1
}

"""
Dictionary mapping input options of RRTMG radiative
components to integers accepted by the RRTMG fortran
code
"""
rrtmg_aerosol_input_dict = {
    'no_aerosol': 0,
    'ecmwf': 6,
    'all_aerosol_properties': 10
}

"""
Dictionary mapping input options of RRTMG radiative
components to integers accepted by the RRTMG fortran
code
"""
rrtmg_random_number_dict = {
    'kissvec': 0,
    'mersenne_twister': 1
}

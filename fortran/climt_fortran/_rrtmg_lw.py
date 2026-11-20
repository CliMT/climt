"""RRTMG longwave through cffi.

Same functions, arguments, defaults and in-place behaviour as the 1.0 Cython
module: ``set_constants`` and ``initialise_rrtm_radiation[_mcica]`` set the
Fortran module's constants and this module's flags; the
``rrtm_calculate_longwave_fluxes*`` calls fill the caller's output arrays.

Array shapes are C order, the reverse of the Fortran declarations: with
``L = num_layers`` and ``C = num_cols``, layer fields are ``(L, C)``,
interface fields ``(L + 1, C)``, ``surface_temp`` is ``(C,)``,
``surf_emissivity`` is ``(16, C)`` (one value per band), ``aerosol_tau`` is
``(16, L, C)``, ``cloud_tau`` is ``(L, C, 16)``, and the McICA sub-column
fields are ``(L, C, 140)`` (one value per g-point). Every array must be
float64 of exactly its shape; any memory layout is accepted.

The calls are plain Python over cffi: at one column, checking and converting
the ~35 array arguments adds about 7 us to the Fortran time, about 1.1 times
the 1.0 Cython call time measured alongside it.
"""

import operator

import numpy as np

from . import _rrtmg_lw_cffi as _compiled
from ._arrays import as_double_scalar_pointer, as_int_pointer, copy_back_outputs, input_pointer, output_pointer

ffi, lib = _compiled.ffi, _compiled.lib
NUM_BANDS = 16  # nbndlw in parrrtm.f90
NUM_G_POINTS = 140  # ngptlw in parrrtm.f90

# Flags set by initialise_rrtm_radiation[_mcica] and passed by pointer to every
# call. RRTMG may rewrite them (an out-of-range overlap method becomes 2, a
# nonzero generator choice becomes 1), and the change sticks, as it did in the
# 1.0 module's globals.
_cloud_overlap_method = ffi.new("int32_t *")
_calculate_change_up_flux = ffi.new("int32_t *")
_cloud_props_flag = ffi.new("int32_t *")
_ice_props_flag = ffi.new("int32_t *")
_liq_droplet_flag = ffi.new("int32_t *")
_permute_seed = ffi.new("int32_t *")
_random_number_generator = ffi.new("int32_t *")


def _int32(value):
    """``value`` as an int; TypeError unless an integer, OverflowError unless it fits int32_t."""
    return as_int_pointer(ffi, value)[0]


def initialise_rrtm_radiation_mcica(cp_d_air,
                                    cloud_overlap_method=1,
                                    calculate_change_up_flux=0,
                                    cloud_props_flag=0,
                                    ice_props_flag=0,
                                    liq_droplet_flag=0,
                                    permute_seed=112,
                                    random_number_generator=0):
    """Initialise RRTMG LW and set the flags of the McICA call."""
    seed, generator = _int32(permute_seed), _int32(random_number_generator)
    initialise_rrtm_radiation(cp_d_air, cloud_overlap_method, calculate_change_up_flux,
                              cloud_props_flag, ice_props_flag, liq_droplet_flag)
    _permute_seed[0], _random_number_generator[0] = seed, generator


def initialise_rrtm_radiation(cp_d_air, cloud_overlap_method=1,
                              calculate_change_up_flux=0,
                              cloud_props_flag=0, ice_props_flag=0,
                              liq_droplet_flag=0):
    """Initialise RRTMG LW and set the flags of both calls. Random cloud
    overlap is the default."""
    # Check every argument before any state changes, as 1.0 did.
    flags = [_int32(flag) for flag in (
        cloud_overlap_method, calculate_change_up_flux, cloud_props_flag, ice_props_flag,
        liq_droplet_flag)]
    lib.rrtmg_lw_ini_wrapper(as_double_scalar_pointer(ffi, cp_d_air))
    (_cloud_overlap_method[0], _calculate_change_up_flux[0], _cloud_props_flag[0],
     _ice_props_flag[0], _liq_droplet_flag[0]) = flags


def set_constants(pi, grav, planck, boltz, clight, avogadro, loschmidt, gas_const,
                  stef_boltz, secs_per_day):
    """Set the physical constants in the Fortran module."""
    lib.rrtmg_set_constants(*[as_double_scalar_pointer(ffi, value) for value in (
        pi, grav, planck, boltz, clight, avogadro, loschmidt, gas_const, stef_boltz, secs_per_day)])


def _common_arguments(num_cols, num_layers, layer_pressure, interface_pressure, layer_temp,
                      interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr,
                      o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity):
    """Pointers to the first 20 arguments of both Fortran calls (up to
    ``surf_emissivity``), checked."""
    layers, interfaces = (num_layers, num_cols), (num_layers + 1, num_cols)
    return [
        as_int_pointer(ffi, num_cols), as_int_pointer(ffi, num_layers),
        _cloud_overlap_method, _calculate_change_up_flux,
        input_pointer(ffi, "layer_pressure", layer_pressure, layers),
        input_pointer(ffi, "interface_pressure", interface_pressure, interfaces),
        input_pointer(ffi, "layer_temp", layer_temp, layers),
        input_pointer(ffi, "interface_temp", interface_temp, interfaces),
        input_pointer(ffi, "surface_temp", surface_temp, (num_cols,)),
        input_pointer(ffi, "h2o_vmr", h2o_vmr, layers),
        input_pointer(ffi, "o3_vmr", o3_vmr, layers),
        input_pointer(ffi, "co2_vmr", co2_vmr, layers),
        input_pointer(ffi, "ch4_vmr", ch4_vmr, layers),
        input_pointer(ffi, "n2o_vmr", n2o_vmr, layers),
        input_pointer(ffi, "o2_vmr", o2_vmr, layers),
        input_pointer(ffi, "cfc11_vmr", cfc11_vmr, layers),
        input_pointer(ffi, "cfc12_vmr", cfc12_vmr, layers),
        input_pointer(ffi, "cfc22_vmr", cfc22_vmr, layers),
        input_pointer(ffi, "ccl4_vmr", ccl4_vmr, layers),
        input_pointer(ffi, "surf_emissivity", surf_emissivity, (NUM_BANDS, num_cols)),
    ]


def _cloud_inputs(num_cols, num_layers, cloud_fraction, cloud_tau, cloud_ice_path,
                  cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius):
    """Pointers to the cloud inputs, checked, in the order of the radiation
    calls: fraction, tau, ice path, liquid path, ice size, droplet radius."""
    layers = (num_layers, num_cols)
    return [
        input_pointer(ffi, "cloud_fraction", cloud_fraction, layers),
        input_pointer(ffi, "cloud_tau", cloud_tau, (num_layers, num_cols, NUM_BANDS)),
        input_pointer(ffi, "cloud_ice_path", cloud_ice_path, layers),
        input_pointer(ffi, "cloud_liq_path", cloud_liq_path, layers),
        input_pointer(ffi, "cloud_ice_eff_size", cloud_ice_eff_size, layers),
        input_pointer(ffi, "cloud_droplet_eff_radius", cloud_droplet_eff_radius, layers),
    ]


def _flux_arguments(num_cols, num_layers, aerosol_tau, upward_longwave_flux,
                    downward_longwave_flux, longwave_heating_rate, up_lw_flux_clearsky,
                    down_lw_flux_clearsky, lw_heating_rate_clearsky, duflx_dt, duflxc_dt,
                    copy_back):
    """Pointers to the last 9 arguments of both Fortran calls (from
    ``aerosol_tau``), checked."""
    layers, interfaces = (num_layers, num_cols), (num_layers + 1, num_cols)
    if _calculate_change_up_flux[0] == 0:
        # As in 1.0: RRTMG does not write them, so small placeholders do.
        duflx_dt = duflxc_dt = np.zeros((num_cols, 1))
        change_shape = (num_cols, 1)
    else:
        change_shape = interfaces
    return [
        input_pointer(ffi, "aerosol_tau", aerosol_tau, (NUM_BANDS, num_layers, num_cols)),
        output_pointer(ffi, "upward_longwave_flux", upward_longwave_flux, interfaces, copy_back),
        output_pointer(ffi, "downward_longwave_flux", downward_longwave_flux, interfaces, copy_back),
        output_pointer(ffi, "longwave_heating_rate", longwave_heating_rate, layers, copy_back),
        output_pointer(ffi, "up_lw_flux_clearsky", up_lw_flux_clearsky, interfaces, copy_back),
        output_pointer(ffi, "down_lw_flux_clearsky", down_lw_flux_clearsky, interfaces, copy_back),
        output_pointer(ffi, "lw_heating_rate_clearsky", lw_heating_rate_clearsky, layers, copy_back),
        output_pointer(ffi, "duflx_dt", duflx_dt, change_shape, copy_back),
        output_pointer(ffi, "duflxc_dt", duflxc_dt, change_shape, copy_back),
    ]


def rrtm_calculate_longwave_fluxes(
        num_cols,
        num_layers,
        layer_pressure,
        interface_pressure,
        layer_temp,
        interface_temp,
        surface_temp,
        h2o_vmr,
        o3_vmr,
        co2_vmr,
        ch4_vmr,
        n2o_vmr,
        o2_vmr,
        cfc11_vmr,
        cfc12_vmr,
        cfc22_vmr,
        ccl4_vmr,
        surf_emissivity,
        cloud_fraction,
        aerosol_tau,
        upward_longwave_flux,
        downward_longwave_flux,
        longwave_heating_rate,
        up_lw_flux_clearsky,
        down_lw_flux_clearsky,
        lw_heating_rate_clearsky,
        cloud_tau=None,
        cloud_ice_path=None,
        cloud_liq_path=None,
        cloud_ice_eff_size=None,
        cloud_droplet_eff_radius=None,
        duflx_dt=None,
        duflxc_dt=None):
    """Longwave fluxes and heating rates without McICA.

    Fills ``upward_longwave_flux`` .. ``lw_heating_rate_clearsky`` in place.
    The cloud arguments are required although they default to None (as in
    1.0). ``duflx_dt`` and ``duflxc_dt`` (``(L + 1, C)``) are required, and
    filled, if ``calculate_change_up_flux`` was set; otherwise they are ignored.
    """
    num_cols, num_layers = operator.index(num_cols), operator.index(num_layers)
    copy_back = []
    lib.rrtmg_lw_nomcica_wrapper(
        *_common_arguments(
            num_cols, num_layers, layer_pressure, interface_pressure, layer_temp,
            interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr,
            o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity),
        _cloud_props_flag, _ice_props_flag, _liq_droplet_flag,
        *_cloud_inputs(
            num_cols, num_layers, cloud_fraction, cloud_tau, cloud_ice_path, cloud_liq_path,
            cloud_ice_eff_size, cloud_droplet_eff_radius),
        *_flux_arguments(
            num_cols, num_layers, aerosol_tau, upward_longwave_flux, downward_longwave_flux,
            longwave_heating_rate, up_lw_flux_clearsky, down_lw_flux_clearsky,
            lw_heating_rate_clearsky, duflx_dt, duflxc_dt, copy_back))
    copy_back_outputs(copy_back)


def rrtm_calculate_longwave_fluxes_mcica(
        rrtm_iplon,
        num_cols,
        num_layers,
        layer_pressure,
        interface_pressure,
        layer_temp,
        interface_temp,
        surface_temp,
        h2o_vmr,
        o3_vmr,
        co2_vmr,
        ch4_vmr,
        n2o_vmr,
        o2_vmr,
        cfc11_vmr,
        cfc12_vmr,
        cfc22_vmr,
        ccl4_vmr,
        surf_emissivity,
        cloud_fraction,
        aerosol_tau,
        upward_longwave_flux,
        downward_longwave_flux,
        longwave_heating_rate,
        up_lw_flux_clearsky,
        down_lw_flux_clearsky,
        lw_heating_rate_clearsky,
        cloud_tau,
        cloud_ice_path,
        cloud_liq_path,
        cloud_ice_eff_size,
        cloud_droplet_eff_radius,
        cloud_fraction_mcica,
        cloud_tau_mcica,
        cloud_ice_path_mcica,
        cloud_liq_path_mcica,
        cloud_ice_eff_size_mcica,
        cloud_droplet_eff_radius_mcica,
        duflx_dt=None,
        duflxc_dt=None):
    """Longwave fluxes and heating rates with McICA cloud sub-columns.

    Draws the sub-columns into ``cloud_*_mcica`` (filled in place), then fills
    ``upward_longwave_flux`` .. ``lw_heating_rate_clearsky`` from them.
    ``duflx_dt`` and ``duflxc_dt`` are as in ``rrtm_calculate_longwave_fluxes``.
    """
    num_cols, num_layers = operator.index(num_cols), operator.index(num_layers)
    layers, g_points = (num_layers, num_cols), (num_layers, num_cols, NUM_G_POINTS)
    copy_back = []
    common = _common_arguments(
        num_cols, num_layers, layer_pressure, interface_pressure, layer_temp,
        interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr,
        o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity)
    fluxes = _flux_arguments(
        num_cols, num_layers, aerosol_tau, upward_longwave_flux, downward_longwave_flux,
        longwave_heating_rate, up_lw_flux_clearsky, down_lw_flux_clearsky,
        lw_heating_rate_clearsky, duflx_dt, duflxc_dt, copy_back)
    # McICA sub-columns, in the order of the radiation call.
    sub_columns = [
        output_pointer(ffi, "cloud_fraction_mcica", cloud_fraction_mcica, g_points, copy_back),
        output_pointer(ffi, "cloud_tau_mcica", cloud_tau_mcica, g_points, copy_back),
        output_pointer(ffi, "cloud_ice_path_mcica", cloud_ice_path_mcica, g_points, copy_back),
        output_pointer(ffi, "cloud_liq_path_mcica", cloud_liq_path_mcica, g_points, copy_back),
        output_pointer(ffi, "cloud_ice_eff_size_mcica", cloud_ice_eff_size_mcica, layers, copy_back),
        output_pointer(ffi, "cloud_droplet_eff_radius_mcica", cloud_droplet_eff_radius_mcica, layers,
                       copy_back),
    ]
    fraction, tau, *paths_and_sizes = _cloud_inputs(
        num_cols, num_layers, cloud_fraction, cloud_tau, cloud_ice_path, cloud_liq_path,
        cloud_ice_eff_size, cloud_droplet_eff_radius)
    sub_fraction, sub_tau, *sub_paths_and_sizes = sub_columns
    num_cols_pointer, num_layers_pointer, layer_pressure_pointer = common[0], common[1], common[4]
    # The sub-column generator takes tau last.
    lib.mcica_subcol_lw_wrapper(
        as_int_pointer(ffi, rrtm_iplon), num_cols_pointer, num_layers_pointer, _cloud_overlap_method,
        _permute_seed, _random_number_generator, layer_pressure_pointer,
        fraction, *paths_and_sizes, tau,
        sub_fraction, *sub_paths_and_sizes, sub_tau)
    lib.rrtmg_lw_mcica_wrapper(
        *common, _cloud_props_flag, _ice_props_flag, _liq_droplet_flag,
        *sub_columns, *fluxes)
    copy_back_outputs(copy_back)

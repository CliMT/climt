"""RRTMG shortwave through cffi.

Same functions, arguments, defaults and in-place behaviour as the 1.0 Cython
module: ``set_constants`` and ``initialise_rrtm_radiation[_mcica]`` set the
Fortran module's constants and this module's flags; the
``rrtm_calculate_shortwave_fluxes*`` calls fill the caller's output arrays.

Array shapes are C order, the reverse of the Fortran declarations: with
``L = num_layers`` and ``C = num_cols``, layer fields are ``(L, C)``,
interface fields ``(L + 1, C)``, surface fields (``surface_temp``, the four
albedos, ``cos_zenith_angle``) ``(C,)``, the aerosol optical properties
``(14, L, C)`` (one value per band), ``aerosol_at_55u`` ``(6, L, C)`` (one
value per ECMWF aerosol type), the cloud optical properties ``(L, C, 14)``,
and the McICA sub-column fields ``(L, C, 112)`` (one value per g-point)
except the ice size and droplet radius, ``(L, C)``. Every array must be
float64 of exactly its shape; any memory layout is accepted.

The solar variability arrays given to ``initialise_rrtm_radiation[_mcica]``
are kept, not copied, and passed to every later call, as 1.0 did:
``facular_sunspot_ampl`` (at least 2 values) is rescaled in place by RRTMG at
each call, by the solar cycle fraction, unless both values are 1;
``solar_var_banded`` (at least 14 values; climt passes 16) is read only.

The calls are plain Python over cffi. At one column, checking the 37 array
arguments (46 with McICA, the solar arrays included) and turning them into
pointers adds about 11 us (13 us) to the Fortran time of about 118 us
(154 us): about 1.09 (1.08) times the bare Fortran call, against 1.03
(1.01) for the 1.0 Cython call, inside the speed gate's 10%. A numba kernel,
as for the longwave calls, was not needed.

cffi releases the GIL around each Fortran call (1.0 held it). The Fortran
module state and this module's flags are not protected against concurrent
Python threads: do not call this module from two threads at once.
"""

import operator

import numpy as np

from . import _rrtmg_sw_cffi as _compiled
from ._arrays import as_double_scalar_pointer, check_array, copy_back_outputs

ffi, lib = _compiled.ffi, _compiled.lib
NUM_BANDS = 14  # nbndsw in parrrsw.f90
NUM_G_POINTS = 112  # ngptsw in parrrsw.f90
NUM_AEROSOL_TYPES = 6  # naerec in parrrsw.f90
_FLOAT64 = np.dtype(np.float64)

# Flags set by initialise_rrtm_radiation[_mcica], passed by pointer to every
# call. RRTMG may rewrite them (an out-of-range overlap method becomes 2, an
# unknown aerosol flag 0, a nonzero generator choice 1), and the change
# sticks, as it did in the 1.0 module's globals.
_flags = np.zeros(8, np.int32)
(_CLOUD_OVERLAP_METHOD, _CLOUD_PROPS_FLAG, _ICE_PROPS_FLAG, _LIQ_DROPLET_FLAG,
 _AEROSOL_INPUT_FLAG, _SOLAR_VAR_FLAG, _PERMUTE_SEED, _RANDOM_NUMBER_GENERATOR) = range(8)
_solar_constant = np.zeros(1)
# (facular_sunspot_ampl, solar_var_banded): the caller's arrays; None until
# initialised.
_solar_arrays = None
_flag = [ffi.from_buffer("int32_t[]", _flags) + index for index in range(len(_flags))]
_solar_constant_pointer = ffi.from_buffer("double[]", _solar_constant)

# The array arguments in the order the radiation calls take them, in three
# runs around the flags and the cloud arrays: (name, shape kind, output).
_STATE = (
    ("layer_pressure", "layers", False), ("interface_pressure", "interfaces", False),
    ("layer_temp", "layers", False), ("interface_temp", "interfaces", False),
    ("surface_temp", "columns", False),
    ("h2o_vmr", "layers", False), ("o3_vmr", "layers", False), ("co2_vmr", "layers", False),
    ("ch4_vmr", "layers", False), ("n2o_vmr", "layers", False), ("o2_vmr", "layers", False),
    ("sfc_direct_sw_albedo", "columns", False), ("sfc_diffuse_sw_albedo", "columns", False),
    ("sfc_direct_nir_albedo", "columns", False), ("sfc_diffuse_nir_albedo", "columns", False),
    ("cos_zenith_angle", "columns", False),
)
_CLOUDS = (
    ("cloud_fraction", "layers", False), ("cloud_tau", "cloud_bands", False),
    ("cloud_single_scat_alb", "cloud_bands", False), ("cloud_asym", "cloud_bands", False),
    ("cloud_fwd_scat_frac", "cloud_bands", False), ("cloud_ice_path", "layers", False),
    ("cloud_liq_path", "layers", False), ("cloud_ice_eff_size", "layers", False),
    ("cloud_droplet_eff_radius", "layers", False),
)
_FLUXES = (
    ("aerosol_tau", "aerosol", False), ("aerosol_single_scat_alb", "aerosol", False),
    ("aerosol_asym", "aerosol", False), ("aerosol_at_55u", "aerosol_types", False),
    ("upward_shortwave_flux", "interfaces", True), ("downward_shortwave_flux", "interfaces", True),
    ("shortwave_heating_rate", "layers", True), ("up_sw_flux_clearsky", "interfaces", True),
    ("down_sw_flux_clearsky", "interfaces", True), ("sw_heating_rate_clearsky", "layers", True),
)
# The McICA sub-columns, in the order of _CLOUDS (the radiation call's).
_SUB_COLUMNS = (
    ("cloud_fraction_mcica", "g_points", True), ("cloud_tau_mcica", "g_points", True),
    ("cloud_single_scat_alb_mcica", "g_points", True), ("cloud_asym_mcica", "g_points", True),
    ("cloud_fwd_scat_frac_mcica", "g_points", True), ("cloud_ice_path_mcica", "g_points", True),
    ("cloud_liq_path_mcica", "g_points", True), ("cloud_ice_eff_size_mcica", "layers", True),
    ("cloud_droplet_eff_radius_mcica", "layers", True),
)
_SOLAR = (("solar_var_banded", "banded", False), ("facular_sunspot_ampl", "facular", True))
_MIN_SOLAR_SIZES = {"facular_sunspot_ampl": 2, "solar_var_banded": NUM_BANDS}


def _shapes(num_cols, num_layers):
    """Shape of each kind of array."""
    return {"layers": (num_layers, num_cols), "interfaces": (num_layers + 1, num_cols),
            "columns": (num_cols,), "aerosol": (NUM_BANDS, num_layers, num_cols),
            "aerosol_types": (NUM_AEROSOL_TYPES, num_layers, num_cols),
            "cloud_bands": (num_layers, num_cols, NUM_BANDS),
            "g_points": (num_layers, num_cols, NUM_G_POINTS)}


def _int32(value):
    """``value`` as an int; TypeError unless an integer, OverflowError unless it fits int32_t
    (Fortran stores it in one)."""
    value = operator.index(value)
    if not -2**31 <= value < 2**31:
        raise OverflowError("%d does not fit in an int32_t" % value)
    return value


def _double(value):
    """``value`` as a float; TypeError unless a real number (as 1.0's ``double`` arguments)."""
    return value if type(value) is float else as_double_scalar_pointer(ffi, value)[0]


def _solar_array(name, array, output):
    """``array`` checked as solar variability array ``name``: 1-D float64 of
    at least its minimum size, writeable if ``output``."""
    check_array(name, array, _FLOAT64, output=output)
    if array.ndim != 1 or array.size < _MIN_SOLAR_SIZES[name]:
        raise ValueError("%s: expected a 1-D array of at least %d values, got shape %s"
                         % (name, _MIN_SOLAR_SIZES[name], array.shape))
    return array


def initialise_rrtm_radiation(cp_d_air,
                              solar_constant,
                              facular_sunspot_ampl,
                              solar_var_banded,
                              cloud_overlap_method=1,
                              cloud_props_flag=0,
                              ice_props_flag=0,
                              liq_droplet_flag=0,
                              aerosol_input_flag=0,
                              solar_var_flag=0):
    """Initialise RRTMG SW and set the flags and solar variability of both
    calls. Random cloud overlap is the default."""
    _initialise(cp_d_air, solar_constant, facular_sunspot_ampl, solar_var_banded,
                (cloud_overlap_method, cloud_props_flag, ice_props_flag, liq_droplet_flag,
                 aerosol_input_flag, solar_var_flag))


def initialise_rrtm_radiation_mcica(cp_d_air,
                                    solar_constant,
                                    facular_sunspot_ampl,
                                    solar_var_banded,
                                    cloud_overlap_method=1,
                                    cloud_props_flag=0,
                                    ice_props_flag=0,
                                    liq_droplet_flag=0,
                                    aerosol_input_flag=0,
                                    solar_var_flag=0,
                                    permute_seed=112,
                                    random_number_generator=0):
    """Initialise RRTMG SW and set the flags and solar variability of the McICA call."""
    _initialise(cp_d_air, solar_constant, facular_sunspot_ampl, solar_var_banded,
                (cloud_overlap_method, cloud_props_flag, ice_props_flag, liq_droplet_flag,
                 aerosol_input_flag, solar_var_flag, permute_seed, random_number_generator))


def _initialise(cp_d_air, solar_constant, facular_sunspot_ampl, solar_var_banded, flags):
    """Check every argument, then initialise RRTMG SW and set ``_flags`` (from
    the start, in their order), the solar constant and the solar arrays."""
    global _solar_arrays
    # Check every argument before any state changes, as 1.0 did.
    flags = [_int32(flag) for flag in flags]
    cp_d_air, solar_constant = as_double_scalar_pointer(ffi, cp_d_air), _double(solar_constant)
    solar_arrays = (_solar_array("facular_sunspot_ampl", facular_sunspot_ampl, output=True),
                    _solar_array("solar_var_banded", solar_var_banded, output=False))
    lib.rrtmg_sw_ini_wrapper(cp_d_air)
    _flags[:len(flags)] = flags
    _solar_constant[0] = solar_constant
    _solar_arrays = solar_arrays


def set_constants(pi, grav, planck, boltz, clight, avogadro, loschmidt, gas_const,
                  stef_boltz, secs_per_day):
    """Set the physical constants in the Fortran module."""
    lib.rrtmg_sw_set_constants(*[as_double_scalar_pointer(ffi, value) for value in (
        pi, grav, planck, boltz, clight, avogadro, loschmidt, gas_const, stef_boltz, secs_per_day)])


def _pointers(arrays, table, shapes, copy_back):
    """double* to each of ``arrays``, checked against ``table`` (an ``_*``
    tuple; shapes from ``shapes``): the array itself if C-contiguous, aligned
    and writeable, else a C-contiguous copy (recorded in ``copy_back`` for
    outputs)."""
    pointers = []
    for array, (name, kind, output) in zip(arrays, table):
        shape = shapes[kind]
        if type(array) is not np.ndarray or array.dtype != _FLOAT64 or array.shape != shape:
            check_array(name, array, _FLOAT64, shape, output)  # raises, or accepts a subclass
        flags = array.flags
        if not (flags.c_contiguous and flags.aligned and flags.writeable):
            check_array(name, array, _FLOAT64, shape, output)  # an output must be writeable
            copy = np.array(array, order="C")  # keeps the caller's values
            if output:
                copy_back.append((array, copy))
            array = copy
        pointers.append(ffi.from_buffer("double[]", array, True))
    return pointers


def _solar_pointers(copy_back):
    """``(solar_var_banded, facular_sunspot_ampl)`` pointers, in the radiation calls' order."""
    if _solar_arrays is None:
        raise ValueError("RRTMG SW is not initialised: call initialise_rrtm_radiation[_mcica] first")
    facular_sunspot_ampl, solar_var_banded = _solar_arrays
    return _pointers((solar_var_banded, facular_sunspot_ampl), _SOLAR,
                     {"facular": facular_sunspot_ampl.shape, "banded": solar_var_banded.shape},
                     copy_back)


def _radiation_arguments(copy_back, num_cols, num_layers, day_of_year, solar_cycle_fraction,
                         flux_adj_sun_earth, layer_pressure, interface_pressure, layer_temp,
                         interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr,
                         o2_vmr, sfc_direct_sw_albedo, sfc_direct_nir_albedo, sfc_diffuse_sw_albedo,
                         sfc_diffuse_nir_albedo, cos_zenith_angle, cloud_fraction,
                         upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate,
                         up_sw_flux_clearsky, down_sw_flux_clearsky, sw_heating_rate_clearsky,
                         aerosol_tau, aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u,
                         cloud_tau, cloud_single_scat_alb, cloud_asym, cloud_fwd_scat_frac,
                         cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                         cloud_droplet_eff_radius):
    """The arguments the two calls share (in their order), checked, as
    pointers for ``rrtmg_sw_[no]mcica_wrapper``: ``(before, clouds, after,
    shapes)``, where ``clouds`` are the nine cloud arrays (in the radiation
    call's order), ``before`` and ``after`` the arguments before and after
    them, and ``shapes`` is from ``_shapes``. ``before[0]``, ``before[1]``
    and ``before[4]`` are num_cols, num_layers and layer_pressure.

    Output copies are recorded in ``copy_back``.
    """
    num_cols, num_layers, day_of_year = _int32(num_cols), _int32(num_layers), _int32(day_of_year)
    solar_cycle_fraction, flux_adj_sun_earth = _double(solar_cycle_fraction), _double(flux_adj_sun_earth)
    shapes = _shapes(num_cols, num_layers)
    # The surface albedos go to Fortran in another order than the call's.
    state = _pointers(
        (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp, h2o_vmr,
         o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo, sfc_diffuse_sw_albedo,
         sfc_direct_nir_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle),
        _STATE, shapes, copy_back)
    clouds = _pointers(
        (cloud_fraction, cloud_tau, cloud_single_scat_alb, cloud_asym, cloud_fwd_scat_frac,
         cloud_ice_path, cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius),
        _CLOUDS, shapes, copy_back)
    fluxes = _pointers(
        (aerosol_tau, aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u, upward_shortwave_flux,
         downward_shortwave_flux, shortwave_heating_rate, up_sw_flux_clearsky,
         down_sw_flux_clearsky, sw_heating_rate_clearsky),
        _FLUXES, shapes, copy_back)
    solar = _solar_pointers(copy_back)
    integers = ffi.new("int32_t[3]", (num_cols, num_layers, day_of_year))
    doubles = ffi.new("double[2]", (flux_adj_sun_earth, solar_cycle_fraction))
    before = [integers, integers + 1, _flag[_CLOUD_OVERLAP_METHOD], _flag[_AEROSOL_INPUT_FLAG],
              *state, doubles, integers + 2, _solar_constant_pointer, _flag[_SOLAR_VAR_FLAG],
              _flag[_CLOUD_PROPS_FLAG], _flag[_ICE_PROPS_FLAG], _flag[_LIQ_DROPLET_FLAG]]
    return before, clouds, fluxes + solar + [doubles + 1], shapes


def rrtm_calculate_shortwave_fluxes(
        num_cols,
        num_layers,
        day_of_year,
        solar_cycle_fraction,
        flux_adj_sun_earth,
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
        sfc_direct_sw_albedo,
        sfc_direct_nir_albedo,
        sfc_diffuse_sw_albedo,
        sfc_diffuse_nir_albedo,
        cos_zenith_angle,
        cloud_fraction,
        upward_shortwave_flux,
        downward_shortwave_flux,
        shortwave_heating_rate,
        up_sw_flux_clearsky,
        down_sw_flux_clearsky,
        sw_heating_rate_clearsky,
        aerosol_tau,
        aerosol_single_scat_alb,
        aerosol_asym,
        aerosol_at_55u,
        cloud_tau,
        cloud_single_scat_alb,
        cloud_asym,
        cloud_fwd_scat_frac,
        cloud_ice_path,
        cloud_liq_path,
        cloud_ice_eff_size,
        cloud_droplet_eff_radius):
    """Shortwave fluxes and heating rates without McICA.

    Fills ``upward_shortwave_flux`` .. ``sw_heating_rate_clearsky`` in place.
    """
    copy_back = []
    before, clouds, after, _ = _radiation_arguments(
        copy_back, num_cols, num_layers, day_of_year, solar_cycle_fraction, flux_adj_sun_earth,
        layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp, h2o_vmr,
        o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo, sfc_direct_nir_albedo,
        sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle, cloud_fraction,
        upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate,
        up_sw_flux_clearsky, down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau,
        aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb,
        cloud_asym, cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
        cloud_droplet_eff_radius)
    lib.rrtmg_sw_nomcica_wrapper(*before, *clouds, *after)
    copy_back_outputs(copy_back)


def rrtm_calculate_shortwave_fluxes_mcica(
        rrtm_iplon,
        num_cols,
        num_layers,
        day_of_year,
        solar_cycle_fraction,
        flux_adj_sun_earth,
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
        sfc_direct_sw_albedo,
        sfc_direct_nir_albedo,
        sfc_diffuse_sw_albedo,
        sfc_diffuse_nir_albedo,
        cos_zenith_angle,
        cloud_fraction,
        upward_shortwave_flux,
        downward_shortwave_flux,
        shortwave_heating_rate,
        up_sw_flux_clearsky,
        down_sw_flux_clearsky,
        sw_heating_rate_clearsky,
        aerosol_tau,
        aerosol_single_scat_alb,
        aerosol_asym,
        aerosol_at_55u,
        cloud_tau,
        cloud_single_scat_alb,
        cloud_asym,
        cloud_fwd_scat_frac,
        cloud_ice_path,
        cloud_liq_path,
        cloud_ice_eff_size,
        cloud_droplet_eff_radius,
        cloud_fraction_mcica,
        cloud_tau_mcica,
        cloud_single_scat_alb_mcica,
        cloud_asym_mcica,
        cloud_fwd_scat_frac_mcica,
        cloud_ice_path_mcica,
        cloud_liq_path_mcica,
        cloud_ice_eff_size_mcica,
        cloud_droplet_eff_radius_mcica):
    """Shortwave fluxes and heating rates with McICA cloud sub-columns.

    Draws the sub-columns into ``cloud_*_mcica`` (filled in place), then fills
    ``upward_shortwave_flux`` .. ``sw_heating_rate_clearsky`` from them. As in
    1.0, the radiation call takes ``cloud_ice_eff_size`` and
    ``cloud_droplet_eff_radius``, not their ``_mcica`` sub-columns.
    """
    rrtm_iplon = ffi.new("int32_t *", _int32(rrtm_iplon))
    copy_back = []
    before, clouds, after, shapes = _radiation_arguments(
        copy_back, num_cols, num_layers, day_of_year, solar_cycle_fraction, flux_adj_sun_earth,
        layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp, h2o_vmr,
        o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo, sfc_direct_nir_albedo,
        sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle, cloud_fraction,
        upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate,
        up_sw_flux_clearsky, down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau,
        aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb,
        cloud_asym, cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
        cloud_droplet_eff_radius)
    sub_columns = _pointers(
        (cloud_fraction_mcica, cloud_tau_mcica, cloud_single_scat_alb_mcica, cloud_asym_mcica,
         cloud_fwd_scat_frac_mcica, cloud_ice_path_mcica, cloud_liq_path_mcica,
         cloud_ice_eff_size_mcica, cloud_droplet_eff_radius_mcica),
        _SUB_COLUMNS, shapes, copy_back)
    # The sub-column generator takes the optical properties (tau .. forward
    # scattering fraction) last.
    lib.mcica_subcol_sw_wrapper(
        rrtm_iplon, before[0], before[1], _flag[_CLOUD_OVERLAP_METHOD], _flag[_PERMUTE_SEED],
        _flag[_RANDOM_NUMBER_GENERATOR], before[4],
        clouds[0], *clouds[5:], *clouds[1:5], sub_columns[0], *sub_columns[5:], *sub_columns[1:5])
    lib.rrtmg_sw_mcica_wrapper(*before, *sub_columns[:7], *clouds[7:], *after)
    copy_back_outputs(copy_back)

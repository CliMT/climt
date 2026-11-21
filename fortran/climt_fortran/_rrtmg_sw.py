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

Each call goes to a numba kernel compiled (and cached) for one signature
only: C-contiguous, writeable float64 arrays, the solar variability arrays
included. The kernel checks the shapes and calls Fortran directly, as for
the longwave calls. Any other argument (a strided or read-only array,
another dtype, None) fails numba's dispatch with TypeError before anything
runs; the call then checks every argument in Python, raising with its name,
and passes C-contiguous copies, writing outputs (``facular_sunspot_ampl``
included) back. The first call in a process with a cold numba cache compiles
its kernel (a second or two for each call), saying so on stderr; later
processes load it from the cache.

The kernels run with the GIL released, both Fortran calls of the McICA call
included (1.0 held the GIL throughout). The Fortran module state and this
module's flags are not protected against concurrent Python threads: do not
call this module from two threads at once.
"""

import operator

import numba
import numpy as np
from numba import types

from . import _rrtmg_sw_cffi as _compiled
from ._arrays import as_double_scalar_pointer, check_array
from ._numba_cffi import announce_compiles, external_function, keep_alive, run_kernel

ffi, lib = _compiled.ffi, _compiled.lib
_mcica_subcol = external_function(_compiled, "mcica_subcol_sw_wrapper")
_rrtmg_sw_mcica = external_function(_compiled, "rrtmg_sw_mcica_wrapper")
_rrtmg_sw_nomcica = external_function(_compiled, "rrtmg_sw_nomcica_wrapper")
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

# The array arguments, in call order: (name, shape kind, output). The arrays
# up to cloud_droplet_eff_radius are the same in both calls; the solar
# variability arrays, kept from the initialise call, come last.
_COMMON = (
    ("layer_pressure", "layers", False), ("interface_pressure", "interfaces", False),
    ("layer_temp", "layers", False), ("interface_temp", "interfaces", False),
    ("surface_temp", "columns", False),
    ("h2o_vmr", "layers", False), ("o3_vmr", "layers", False), ("co2_vmr", "layers", False),
    ("ch4_vmr", "layers", False), ("n2o_vmr", "layers", False), ("o2_vmr", "layers", False),
    ("sfc_direct_sw_albedo", "columns", False), ("sfc_direct_nir_albedo", "columns", False),
    ("sfc_diffuse_sw_albedo", "columns", False), ("sfc_diffuse_nir_albedo", "columns", False),
    ("cos_zenith_angle", "columns", False), ("cloud_fraction", "layers", False),
    ("upward_shortwave_flux", "interfaces", True), ("downward_shortwave_flux", "interfaces", True),
    ("shortwave_heating_rate", "layers", True), ("up_sw_flux_clearsky", "interfaces", True),
    ("down_sw_flux_clearsky", "interfaces", True), ("sw_heating_rate_clearsky", "layers", True),
    ("aerosol_tau", "aerosol", False), ("aerosol_single_scat_alb", "aerosol", False),
    ("aerosol_asym", "aerosol", False), ("aerosol_at_55u", "aerosol_types", False),
    ("cloud_tau", "cloud_bands", False), ("cloud_single_scat_alb", "cloud_bands", False),
    ("cloud_asym", "cloud_bands", False), ("cloud_fwd_scat_frac", "cloud_bands", False),
    ("cloud_ice_path", "layers", False), ("cloud_liq_path", "layers", False),
    ("cloud_ice_eff_size", "layers", False), ("cloud_droplet_eff_radius", "layers", False),
)
_SUB_COLUMNS = (
    ("cloud_fraction_mcica", "g_points", True), ("cloud_tau_mcica", "g_points", True),
    ("cloud_single_scat_alb_mcica", "g_points", True), ("cloud_asym_mcica", "g_points", True),
    ("cloud_fwd_scat_frac_mcica", "g_points", True), ("cloud_ice_path_mcica", "g_points", True),
    ("cloud_liq_path_mcica", "g_points", True), ("cloud_ice_eff_size_mcica", "layers", True),
    ("cloud_droplet_eff_radius_mcica", "layers", True),
)
_SOLAR = (("solar_var_banded", "banded", False), ("facular_sunspot_ampl", "facular", True))
_NOMCICA_ARRAYS = _COMMON + _SOLAR
_MCICA_ARRAYS = _COMMON + _SUB_COLUMNS + _SOLAR
# The array arguments of a public function, by name, from its locals(): those
# of _*_ARRAYS without the solar variability arrays, which come from _solar().
_nomcica_call_arrays = operator.itemgetter(*[name for name, _, _ in _COMMON])
_mcica_call_arrays = operator.itemgetter(*[name for name, _, _ in _COMMON + _SUB_COLUMNS])
_NDIM = {"layers": 2, "interfaces": 2, "columns": 1, "aerosol": 3, "aerosol_types": 3,
         "cloud_bands": 3, "g_points": 3, "banded": 1, "facular": 1}
_MIN_SOLAR_SIZES = {"facular_sunspot_ampl": 2, "solar_var_banded": NUM_BANDS}


def _shapes(num_cols, num_layers):
    """Shape of each kind of array (the kernels check the same shapes; the
    solar variability arrays keep the shapes they were given with)."""
    facular_sunspot_ampl, solar_var_banded = _solar_arrays
    return {"layers": (num_layers, num_cols), "interfaces": (num_layers + 1, num_cols),
            "columns": (num_cols,), "aerosol": (NUM_BANDS, num_layers, num_cols),
            "aerosol_types": (NUM_AEROSOL_TYPES, num_layers, num_cols),
            "cloud_bands": (num_layers, num_cols, NUM_BANDS),
            "g_points": (num_layers, num_cols, NUM_G_POINTS),
            "banded": solar_var_banded.shape, "facular": facular_sunspot_ampl.shape}


def _signature(leading_integers, arrays):
    """The one signature a kernel is compiled for: ``leading_integers`` int64
    scalars, solar_cycle_fraction and flux_adj_sun_earth (float64), the flags,
    the solar constant (a 1-value array), then C-contiguous writeable float64
    ``arrays``."""
    return types.boolean(
        *[types.int64] * leading_integers, types.float64, types.float64,
        types.Array(types.int32, 1, "C"), types.Array(types.float64, 1, "C"),
        *[types.Array(types.float64, _NDIM[kind], "C") for _, kind, _ in arrays])


_NOMCICA_SIGNATURE = _signature(3, _NOMCICA_ARRAYS)  # num_cols, num_layers, day_of_year
_MCICA_SIGNATURE = _signature(4, _MCICA_ARRAYS)  # rrtm_iplon, num_cols, num_layers, day_of_year


def _int32(value):
    """``value`` as an int; TypeError unless an integer, OverflowError unless it fits int32_t
    (the kernels and Fortran store it in one)."""
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


def _solar():
    """``(solar_var_banded, facular_sunspot_ampl)``, the kernels' last arguments."""
    if _solar_arrays is None:
        raise ValueError("RRTMG SW is not initialised: call initialise_rrtm_radiation[_mcica] first")
    facular_sunspot_ampl, solar_var_banded = _solar_arrays
    return solar_var_banded, facular_sunspot_ampl


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
    num_cols, num_layers = _int32(num_cols), _int32(num_layers)
    leading = (num_cols, num_layers, _int32(day_of_year), _double(solar_cycle_fraction),
               _double(flux_adj_sun_earth), _flags, _solar_constant)
    # The parameter names must match the table's (test_argument_order checks).
    run_kernel(_nomcica_kernel, _NOMCICA_SIGNATURE, leading, _nomcica_call_arrays(locals()) + _solar(),
               _NOMCICA_ARRAYS, _shapes, num_cols, num_layers)


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
    num_cols, num_layers = _int32(num_cols), _int32(num_layers)
    leading = (_int32(rrtm_iplon), num_cols, num_layers, _int32(day_of_year),
               _double(solar_cycle_fraction), _double(flux_adj_sun_earth), _flags, _solar_constant)
    # The parameter names must match the table's (test_argument_order checks).
    run_kernel(_mcica_kernel, _MCICA_SIGNATURE, leading, _mcica_call_arrays(locals()) + _solar(),
               _MCICA_ARRAYS, _shapes, num_cols, num_layers)


# The kernels take the arrays in the order of the Python calls (the
# _*_ARRAYS tables) and pass them to Fortran in its order. Integers and the
# two double scalars go to Fortran through small arrays; the flags and the
# solar constant through _flags and _solar_constant themselves, so Fortran's
# rewrites stick. Each kernel gathers the arrays of _COMMON into one tuple,
# ``common``, and passes it on with ``*common``.
# The order of _COMMON is written out again in the parameters of
# _shapes_ok, _flux_pointers and each kernel, and in each ``common``;
# test_argument_order checks every Fortran call against 1.0's.

@numba.njit(cache=True)
def _shapes_ok(num_cols, num_layers, layer_pressure, interface_pressure, layer_temp, interface_temp,
               surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo,
               sfc_direct_nir_albedo, sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle,
               cloud_fraction, upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate,
               up_sw_flux_clearsky, down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau,
               aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb,
               cloud_asym, cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
               cloud_droplet_eff_radius, solar_var_banded, facular_sunspot_ampl):
    """Whether the arrays of _COMMON have the shapes of ``_shapes``, and the
    solar variability arrays at least their minimum sizes."""
    layers = (num_layers, num_cols)
    interfaces = (num_layers + 1, num_cols)
    columns = (num_cols,)
    aerosol = (NUM_BANDS, num_layers, num_cols)
    cloud_bands = (num_layers, num_cols, NUM_BANDS)
    for array in (layer_pressure, layer_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr,
                  cloud_fraction, shortwave_heating_rate, sw_heating_rate_clearsky, cloud_ice_path,
                  cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius):
        if array.shape != layers:
            return False
    for array in (interface_pressure, interface_temp, upward_shortwave_flux, downward_shortwave_flux,
                  up_sw_flux_clearsky, down_sw_flux_clearsky):
        if array.shape != interfaces:
            return False
    for array in (surface_temp, sfc_direct_sw_albedo, sfc_direct_nir_albedo, sfc_diffuse_sw_albedo,
                  sfc_diffuse_nir_albedo, cos_zenith_angle):
        if array.shape != columns:
            return False
    for array in (aerosol_tau, aerosol_single_scat_alb, aerosol_asym):
        if array.shape != aerosol:
            return False
    for array in (cloud_tau, cloud_single_scat_alb, cloud_asym, cloud_fwd_scat_frac):
        if array.shape != cloud_bands:
            return False
    return (aerosol_at_55u.shape == (NUM_AEROSOL_TYPES, num_layers, num_cols)
            and solar_var_banded.size >= NUM_BANDS and facular_sunspot_ampl.size >= 2)


@numba.njit(cache=True)
def _flux_pointers(sizes, doubles, flags, solar_constant, layer_pressure, interface_pressure, layer_temp,
                   interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr,
                   sfc_direct_sw_albedo, sfc_direct_nir_albedo, sfc_diffuse_sw_albedo,
                   sfc_diffuse_nir_albedo, cos_zenith_angle, cloud_fraction, upward_shortwave_flux,
                   downward_shortwave_flux, shortwave_heating_rate, up_sw_flux_clearsky,
                   down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau, aerosol_single_scat_alb,
                   aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb, cloud_asym,
                   cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                   cloud_droplet_eff_radius, solar_var_banded, facular_sunspot_ampl):
    """Pointers for ``rrtmg_sw_[no]mcica_wrapper`` as ``(before, after)``: the
    arguments before and after its nine cloud arrays (which differ between
    the two). ``sizes`` holds num_cols, num_layers and day_of_year;
    ``doubles`` flux_adj_sun_earth and solar_cycle_fraction. The caller keeps
    both alive past the Fortran call."""
    # The surface albedos go to Fortran in another order than the call's.
    before = (
        ffi.from_buffer(sizes[0:]), ffi.from_buffer(sizes[1:]),
        ffi.from_buffer(flags[_CLOUD_OVERLAP_METHOD:]), ffi.from_buffer(flags[_AEROSOL_INPUT_FLAG:]),
        ffi.from_buffer(layer_pressure), ffi.from_buffer(interface_pressure),
        ffi.from_buffer(layer_temp), ffi.from_buffer(interface_temp), ffi.from_buffer(surface_temp),
        ffi.from_buffer(h2o_vmr), ffi.from_buffer(o3_vmr), ffi.from_buffer(co2_vmr),
        ffi.from_buffer(ch4_vmr), ffi.from_buffer(n2o_vmr), ffi.from_buffer(o2_vmr),
        ffi.from_buffer(sfc_direct_sw_albedo), ffi.from_buffer(sfc_diffuse_sw_albedo),
        ffi.from_buffer(sfc_direct_nir_albedo), ffi.from_buffer(sfc_diffuse_nir_albedo),
        ffi.from_buffer(cos_zenith_angle), ffi.from_buffer(doubles[0:]), ffi.from_buffer(sizes[2:]),
        ffi.from_buffer(solar_constant), ffi.from_buffer(flags[_SOLAR_VAR_FLAG:]),
        ffi.from_buffer(flags[_CLOUD_PROPS_FLAG:]), ffi.from_buffer(flags[_ICE_PROPS_FLAG:]),
        ffi.from_buffer(flags[_LIQ_DROPLET_FLAG:]))
    after = (
        ffi.from_buffer(aerosol_tau), ffi.from_buffer(aerosol_single_scat_alb),
        ffi.from_buffer(aerosol_asym), ffi.from_buffer(aerosol_at_55u),
        ffi.from_buffer(upward_shortwave_flux), ffi.from_buffer(downward_shortwave_flux),
        ffi.from_buffer(shortwave_heating_rate), ffi.from_buffer(up_sw_flux_clearsky),
        ffi.from_buffer(down_sw_flux_clearsky), ffi.from_buffer(sw_heating_rate_clearsky),
        ffi.from_buffer(solar_var_banded), ffi.from_buffer(facular_sunspot_ampl),
        ffi.from_buffer(doubles[1:]))
    return before, after


@numba.njit(cache=True)
def _cloud_pointers(fraction, tau, single_scat_alb, asym, fwd_scat_frac, ice_path, liq_path, ice_eff_size,
                    droplet_eff_radius):
    """Pointers to nine cloud arrays, in the order of ``rrtmg_sw_[no]mcica_wrapper``."""
    return (ffi.from_buffer(fraction), ffi.from_buffer(tau), ffi.from_buffer(single_scat_alb),
            ffi.from_buffer(asym), ffi.from_buffer(fwd_scat_frac), ffi.from_buffer(ice_path),
            ffi.from_buffer(liq_path), ffi.from_buffer(ice_eff_size), ffi.from_buffer(droplet_eff_radius))


def _nomcica_kernel(num_cols, num_layers, day_of_year, solar_cycle_fraction, flux_adj_sun_earth, flags,
                    solar_constant, layer_pressure, interface_pressure, layer_temp, interface_temp,
                    surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo,
                    sfc_direct_nir_albedo, sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle,
                    cloud_fraction, upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate,
                    up_sw_flux_clearsky, down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau,
                    aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb,
                    cloud_asym, cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                    cloud_droplet_eff_radius, solar_var_banded, facular_sunspot_ampl):
    common = (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp, h2o_vmr, o3_vmr,
              co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo, sfc_direct_nir_albedo,
              sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle, cloud_fraction,
              upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate, up_sw_flux_clearsky,
              down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau, aerosol_single_scat_alb,
              aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb, cloud_asym,
              cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
              cloud_droplet_eff_radius)
    if not _shapes_ok(num_cols, num_layers, *common, solar_var_banded, facular_sunspot_ampl):
        return False
    sizes = np.empty(3, np.int32)
    sizes[0] = num_cols
    sizes[1] = num_layers
    sizes[2] = day_of_year
    doubles = np.empty(2)
    doubles[0] = flux_adj_sun_earth
    doubles[1] = solar_cycle_fraction
    before, after = _flux_pointers(sizes, doubles, flags, solar_constant, *common, solar_var_banded,
                                   facular_sunspot_ampl)
    clouds = _cloud_pointers(cloud_fraction, cloud_tau, cloud_single_scat_alb, cloud_asym,
                             cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                             cloud_droplet_eff_radius)
    _rrtmg_sw_nomcica(*before, *clouds, *after)
    keep_alive(sizes)
    keep_alive(doubles)
    return True


def _mcica_kernel(rrtm_iplon, num_cols, num_layers, day_of_year, solar_cycle_fraction, flux_adj_sun_earth,
                  flags, solar_constant, layer_pressure, interface_pressure, layer_temp, interface_temp,
                  surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo,
                  sfc_direct_nir_albedo, sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle,
                  cloud_fraction, upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate,
                  up_sw_flux_clearsky, down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau,
                  aerosol_single_scat_alb, aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb,
                  cloud_asym, cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                  cloud_droplet_eff_radius, cloud_fraction_mcica, cloud_tau_mcica,
                  cloud_single_scat_alb_mcica, cloud_asym_mcica, cloud_fwd_scat_frac_mcica,
                  cloud_ice_path_mcica, cloud_liq_path_mcica, cloud_ice_eff_size_mcica,
                  cloud_droplet_eff_radius_mcica, solar_var_banded, facular_sunspot_ampl):
    common = (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp, h2o_vmr, o3_vmr,
              co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, sfc_direct_sw_albedo, sfc_direct_nir_albedo,
              sfc_diffuse_sw_albedo, sfc_diffuse_nir_albedo, cos_zenith_angle, cloud_fraction,
              upward_shortwave_flux, downward_shortwave_flux, shortwave_heating_rate, up_sw_flux_clearsky,
              down_sw_flux_clearsky, sw_heating_rate_clearsky, aerosol_tau, aerosol_single_scat_alb,
              aerosol_asym, aerosol_at_55u, cloud_tau, cloud_single_scat_alb, cloud_asym,
              cloud_fwd_scat_frac, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
              cloud_droplet_eff_radius)
    g_points = (num_layers, num_cols, NUM_G_POINTS)
    layers = (num_layers, num_cols)
    if not _shapes_ok(num_cols, num_layers, *common, solar_var_banded, facular_sunspot_ampl):
        return False
    for array in (cloud_fraction_mcica, cloud_tau_mcica, cloud_single_scat_alb_mcica, cloud_asym_mcica,
                  cloud_fwd_scat_frac_mcica, cloud_ice_path_mcica, cloud_liq_path_mcica):
        if array.shape != g_points:
            return False
    if cloud_ice_eff_size_mcica.shape != layers or cloud_droplet_eff_radius_mcica.shape != layers:
        return False
    sizes = np.empty(4, np.int32)
    sizes[0] = rrtm_iplon
    sizes[1] = num_cols
    sizes[2] = num_layers
    sizes[3] = day_of_year
    doubles = np.empty(2)
    doubles[0] = flux_adj_sun_earth
    doubles[1] = solar_cycle_fraction
    # The sub-column generator takes the optical properties (tau .. forward
    # scattering fraction) last.
    _mcica_subcol(
        ffi.from_buffer(sizes[0:]), ffi.from_buffer(sizes[1:]), ffi.from_buffer(sizes[2:]),
        ffi.from_buffer(flags[_CLOUD_OVERLAP_METHOD:]), ffi.from_buffer(flags[_PERMUTE_SEED:]),
        ffi.from_buffer(flags[_RANDOM_NUMBER_GENERATOR:]), ffi.from_buffer(layer_pressure),
        ffi.from_buffer(cloud_fraction), ffi.from_buffer(cloud_ice_path), ffi.from_buffer(cloud_liq_path),
        ffi.from_buffer(cloud_ice_eff_size), ffi.from_buffer(cloud_droplet_eff_radius),
        ffi.from_buffer(cloud_tau), ffi.from_buffer(cloud_single_scat_alb), ffi.from_buffer(cloud_asym),
        ffi.from_buffer(cloud_fwd_scat_frac), ffi.from_buffer(cloud_fraction_mcica),
        ffi.from_buffer(cloud_ice_path_mcica), ffi.from_buffer(cloud_liq_path_mcica),
        ffi.from_buffer(cloud_ice_eff_size_mcica), ffi.from_buffer(cloud_droplet_eff_radius_mcica),
        ffi.from_buffer(cloud_tau_mcica), ffi.from_buffer(cloud_single_scat_alb_mcica),
        ffi.from_buffer(cloud_asym_mcica), ffi.from_buffer(cloud_fwd_scat_frac_mcica))
    before, after = _flux_pointers(sizes[1:], doubles, flags, solar_constant, *common, solar_var_banded,
                                   facular_sunspot_ampl)
    # As in 1.0, the radiation call takes the non-McICA ice size and droplet
    # radius.
    clouds = _cloud_pointers(cloud_fraction_mcica, cloud_tau_mcica, cloud_single_scat_alb_mcica,
                             cloud_asym_mcica, cloud_fwd_scat_frac_mcica, cloud_ice_path_mcica,
                             cloud_liq_path_mcica, cloud_ice_eff_size, cloud_droplet_eff_radius)
    _rrtmg_sw_mcica(*before, *clouds, *after)
    keep_alive(sizes)
    keep_alive(doubles)
    return True


announce_compiles(_nomcica_kernel, "the RRTMG shortwave wrapper")
announce_compiles(_mcica_kernel, "the RRTMG shortwave McICA wrapper")

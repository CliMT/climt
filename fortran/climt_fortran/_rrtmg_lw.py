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

Each call goes to a numba kernel compiled (and cached) for one signature
only: C-contiguous, writeable float64 arrays. The kernel checks the shapes
and calls Fortran directly. Any other argument (a strided or read-only array,
another dtype, None) fails numba's dispatch with TypeError before anything
runs; the call then checks every argument in Python, raising with its name,
and passes C-contiguous copies, writing outputs back. Converting the ~35
arguments to cffi pointers in Python instead made the calls about 1.10 times
slower than the Cython version. The first call in a process with a cold
numba cache compiles its kernel (about 0.7 s and 0.4 s for the two calls),
saying so on stderr; later processes load it from the cache.

The kernels run with the GIL released, both Fortran calls of the McICA call
included (1.0 held the GIL for the sub-column generator). The Fortran module
state is not protected against concurrent Python threads, as with 1.0's
``with nogil`` main call: do not call this module from two threads at once.
"""

import operator

import numba
import numpy as np
from numba import types

from . import _rrtmg_lw_cffi as _compiled
from ._arrays import as_double_scalar_pointer, contiguous_arguments, copy_back_outputs
from ._numba_cffi import announce_compiles, compiled_kernel, external_function, keep_alive

ffi, lib = _compiled.ffi, _compiled.lib
_mcica_subcol = external_function(_compiled, "mcica_subcol_lw_wrapper")
_rrtmg_lw_mcica = external_function(_compiled, "rrtmg_lw_mcica_wrapper")
_rrtmg_lw_nomcica = external_function(_compiled, "rrtmg_lw_nomcica_wrapper")
NUM_BANDS = 16  # nbndlw in parrrtm.f90
NUM_G_POINTS = 140  # ngptlw in parrrtm.f90

# Flags set by initialise_rrtm_radiation[_mcica], passed by pointer to every
# call. RRTMG may rewrite them (an out-of-range overlap method becomes 2, a
# nonzero generator choice becomes 1), and the change sticks, as it did in the
# 1.0 module's globals.
_flags = np.zeros(7, np.int32)
(_CLOUD_OVERLAP_METHOD, _CALCULATE_CHANGE_UP_FLUX, _CLOUD_PROPS_FLAG, _ICE_PROPS_FLAG,
 _LIQ_DROPLET_FLAG, _PERMUTE_SEED, _RANDOM_NUMBER_GENERATOR) = range(7)

# The array arguments, in call order: (name, shape kind, output). The arrays
# up to cloud_droplet_eff_radius are the same in both calls.
_COMMON = (
    ("layer_pressure", "layers", False), ("interface_pressure", "interfaces", False),
    ("layer_temp", "layers", False), ("interface_temp", "interfaces", False),
    ("surface_temp", "columns", False),
    ("h2o_vmr", "layers", False), ("o3_vmr", "layers", False), ("co2_vmr", "layers", False),
    ("ch4_vmr", "layers", False), ("n2o_vmr", "layers", False), ("o2_vmr", "layers", False),
    ("cfc11_vmr", "layers", False), ("cfc12_vmr", "layers", False),
    ("cfc22_vmr", "layers", False), ("ccl4_vmr", "layers", False),
    ("surf_emissivity", "bands", False), ("cloud_fraction", "layers", False),
    ("aerosol_tau", "aerosol", False),
    ("upward_longwave_flux", "interfaces", True), ("downward_longwave_flux", "interfaces", True),
    ("longwave_heating_rate", "layers", True), ("up_lw_flux_clearsky", "interfaces", True),
    ("down_lw_flux_clearsky", "interfaces", True), ("lw_heating_rate_clearsky", "layers", True),
    ("cloud_tau", "cloud_tau", False), ("cloud_ice_path", "layers", False),
    ("cloud_liq_path", "layers", False), ("cloud_ice_eff_size", "layers", False),
    ("cloud_droplet_eff_radius", "layers", False),
)
_SUB_COLUMNS = (
    ("cloud_fraction_mcica", "g_points", True), ("cloud_tau_mcica", "g_points", True),
    ("cloud_ice_path_mcica", "g_points", True), ("cloud_liq_path_mcica", "g_points", True),
    ("cloud_ice_eff_size_mcica", "layers", True), ("cloud_droplet_eff_radius_mcica", "layers", True),
)
_CHANGE_UP_FLUX = (("duflx_dt", "change", True), ("duflxc_dt", "change", True))
_NOMCICA_ARRAYS = _COMMON + _CHANGE_UP_FLUX
_MCICA_ARRAYS = _COMMON + _SUB_COLUMNS + _CHANGE_UP_FLUX
_NDIM = {"layers": 2, "interfaces": 2, "columns": 1, "bands": 2, "aerosol": 3, "cloud_tau": 3,
         "g_points": 3, "change": 2}


def _shapes(num_cols, num_layers):
    """Shape of each kind of array (the kernels check the same shapes)."""
    change = (num_layers + 1, num_cols) if _flags[_CALCULATE_CHANGE_UP_FLUX] else (num_cols, 1)
    return {"layers": (num_layers, num_cols), "interfaces": (num_layers + 1, num_cols),
            "columns": (num_cols,), "bands": (NUM_BANDS, num_cols),
            "aerosol": (NUM_BANDS, num_layers, num_cols),
            "cloud_tau": (num_layers, num_cols, NUM_BANDS),
            "g_points": (num_layers, num_cols, NUM_G_POINTS), "change": change}


def _signature(leading_integers, arrays):
    """The one signature a kernel is compiled for: ``leading_integers`` int64
    scalars, the flags, then C-contiguous writeable float64 ``arrays``."""
    return types.boolean(
        *[types.int64] * leading_integers, types.Array(types.int32, 1, "C"),
        *[types.Array(types.float64, _NDIM[kind], "C") for _, kind, _ in arrays])


_NOMCICA_SIGNATURE = _signature(2, _NOMCICA_ARRAYS)  # num_cols, num_layers
_MCICA_SIGNATURE = _signature(3, _MCICA_ARRAYS)  # rrtm_iplon, num_cols, num_layers


def _int32(value):
    """``value`` as an int; TypeError unless an integer, OverflowError unless it fits int32_t
    (the kernels and Fortran store it in one)."""
    value = operator.index(value)
    if not -2**31 <= value < 2**31:
        raise OverflowError("%d does not fit in an int32_t" % value)
    return value


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
    _flags[_PERMUTE_SEED], _flags[_RANDOM_NUMBER_GENERATOR] = seed, generator


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
    _flags[_CLOUD_OVERLAP_METHOD:_LIQ_DROPLET_FLAG + 1] = flags


def set_constants(pi, grav, planck, boltz, clight, avogadro, loschmidt, gas_const,
                  stef_boltz, secs_per_day):
    """Set the physical constants in the Fortran module."""
    lib.rrtmg_set_constants(*[as_double_scalar_pointer(ffi, value) for value in (
        pi, grav, planck, boltz, clight, avogadro, loschmidt, gas_const, stef_boltz, secs_per_day)])


def _call(kernel, signature, leading, arrays, table, num_cols, num_layers):
    """Run the kernel compiled for ``signature`` on ``(*leading, *arrays)``.

    If it does not take the arguments (TypeError from numba's dispatch, or
    False: a wrong shape), check every array against ``table`` (an
    ``_*_ARRAYS`` tuple), raising with the name of the bad one, and run it on
    C-contiguous writeable copies, writing outputs back.
    """
    kernel = compiled_kernel(kernel, signature)
    try:
        if kernel(*leading, *arrays):
            return
    except TypeError:
        pass
    shapes = _shapes(num_cols, num_layers)
    copy_back = []
    arrays = contiguous_arguments(
        arrays, [(name, shapes[kind], output) for name, kind, output in table], copy_back)
    if not kernel(*leading, *arrays):
        raise RuntimeError("the arguments passed the checks but not the kernel's")
    copy_back_outputs(copy_back)


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
    num_cols, num_layers = _int32(num_cols), _int32(num_layers)
    if _flags[_CALCULATE_CHANGE_UP_FLUX] == 0:
        # As in 1.0: RRTMG does not write them, so small placeholders do.
        duflx_dt = duflxc_dt = np.zeros((num_cols, 1))
    _call(_nomcica_kernel, _NOMCICA_SIGNATURE, (num_cols, num_layers, _flags),
          (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp,
           h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr,
           ccl4_vmr, surf_emissivity, cloud_fraction, aerosol_tau, upward_longwave_flux,
           downward_longwave_flux, longwave_heating_rate, up_lw_flux_clearsky,
           down_lw_flux_clearsky, lw_heating_rate_clearsky, cloud_tau, cloud_ice_path,
           cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius, duflx_dt, duflxc_dt),
          _NOMCICA_ARRAYS, num_cols, num_layers)


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
    rrtm_iplon = _int32(rrtm_iplon)
    num_cols, num_layers = _int32(num_cols), _int32(num_layers)
    if _flags[_CALCULATE_CHANGE_UP_FLUX] == 0:
        duflx_dt = duflxc_dt = np.zeros((num_cols, 1))
    _call(_mcica_kernel, _MCICA_SIGNATURE, (rrtm_iplon, num_cols, num_layers, _flags),
          (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp,
           h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr,
           ccl4_vmr, surf_emissivity, cloud_fraction, aerosol_tau, upward_longwave_flux,
           downward_longwave_flux, longwave_heating_rate, up_lw_flux_clearsky,
           down_lw_flux_clearsky, lw_heating_rate_clearsky, cloud_tau, cloud_ice_path,
           cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius, cloud_fraction_mcica,
           cloud_tau_mcica, cloud_ice_path_mcica, cloud_liq_path_mcica, cloud_ice_eff_size_mcica,
           cloud_droplet_eff_radius_mcica, duflx_dt, duflxc_dt),
          _MCICA_ARRAYS, num_cols, num_layers)


# The kernels take the arrays in the order of the Python calls (the
# _*_ARRAYS tables) and pass them to Fortran in its order. Integers go to
# Fortran through a small int32 array; the flags through _flags itself, so
# Fortran's rewrites stick. Each kernel gathers the arrays of _COMMON into
# one tuple, ``common``, and passes it on with ``*common``.

@numba.njit(cache=True)
def _shapes_ok(num_cols, num_layers, flags, layer_pressure, interface_pressure, layer_temp,
               interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr,
               cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity, cloud_fraction,
               aerosol_tau, upward_longwave_flux, downward_longwave_flux, longwave_heating_rate,
               up_lw_flux_clearsky, down_lw_flux_clearsky, lw_heating_rate_clearsky, cloud_tau,
               cloud_ice_path, cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius,
               duflx_dt, duflxc_dt):
    """Whether the arrays of _COMMON and _CHANGE_UP_FLUX have the shapes of ``_shapes``."""
    layers = (num_layers, num_cols)
    interfaces = (num_layers + 1, num_cols)
    change = interfaces if flags[_CALCULATE_CHANGE_UP_FLUX] != 0 else (num_cols, 1)
    for array in (layer_pressure, layer_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr,
                  cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, cloud_fraction, longwave_heating_rate,
                  lw_heating_rate_clearsky, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                  cloud_droplet_eff_radius):
        if array.shape != layers:
            return False
    for array in (interface_pressure, interface_temp, upward_longwave_flux, downward_longwave_flux,
                  up_lw_flux_clearsky, down_lw_flux_clearsky):
        if array.shape != interfaces:
            return False
    return (surface_temp.shape == (num_cols,) and surf_emissivity.shape == (NUM_BANDS, num_cols)
            and aerosol_tau.shape == (NUM_BANDS, num_layers, num_cols)
            and cloud_tau.shape == (num_layers, num_cols, NUM_BANDS)
            and duflx_dt.shape == change and duflxc_dt.shape == change)


@numba.njit(cache=True)
def _flux_pointers(sizes, flags, layer_pressure, interface_pressure, layer_temp, interface_temp,
                   surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, cfc11_vmr,
                   cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity, cloud_fraction, aerosol_tau,
                   upward_longwave_flux, downward_longwave_flux, longwave_heating_rate,
                   up_lw_flux_clearsky, down_lw_flux_clearsky, lw_heating_rate_clearsky,
                   cloud_tau, cloud_ice_path, cloud_liq_path, cloud_ice_eff_size,
                   cloud_droplet_eff_radius, duflx_dt, duflxc_dt):
    """Pointers for ``rrtmg_lw_[no]mcica_wrapper`` as ``(before, after)``: the
    arguments before and after its six cloud arrays (which differ between the
    two). ``sizes`` holds num_cols then num_layers. The caller keeps ``sizes``
    alive past the Fortran call."""
    before = (
        ffi.from_buffer(sizes[0:]), ffi.from_buffer(sizes[1:]),
        ffi.from_buffer(flags[_CLOUD_OVERLAP_METHOD:]), ffi.from_buffer(flags[_CALCULATE_CHANGE_UP_FLUX:]),
        ffi.from_buffer(layer_pressure), ffi.from_buffer(interface_pressure),
        ffi.from_buffer(layer_temp), ffi.from_buffer(interface_temp), ffi.from_buffer(surface_temp),
        ffi.from_buffer(h2o_vmr), ffi.from_buffer(o3_vmr), ffi.from_buffer(co2_vmr),
        ffi.from_buffer(ch4_vmr), ffi.from_buffer(n2o_vmr), ffi.from_buffer(o2_vmr),
        ffi.from_buffer(cfc11_vmr), ffi.from_buffer(cfc12_vmr), ffi.from_buffer(cfc22_vmr),
        ffi.from_buffer(ccl4_vmr), ffi.from_buffer(surf_emissivity),
        ffi.from_buffer(flags[_CLOUD_PROPS_FLAG:]), ffi.from_buffer(flags[_ICE_PROPS_FLAG:]),
        ffi.from_buffer(flags[_LIQ_DROPLET_FLAG:]))
    after = (
        ffi.from_buffer(aerosol_tau),
        ffi.from_buffer(upward_longwave_flux), ffi.from_buffer(downward_longwave_flux),
        ffi.from_buffer(longwave_heating_rate), ffi.from_buffer(up_lw_flux_clearsky),
        ffi.from_buffer(down_lw_flux_clearsky), ffi.from_buffer(lw_heating_rate_clearsky),
        ffi.from_buffer(duflx_dt), ffi.from_buffer(duflxc_dt))
    return before, after


@numba.njit(cache=True)
def _cloud_pointers(fraction, tau, ice_path, liq_path, ice_eff_size, droplet_eff_radius):
    """Pointers to six cloud arrays, in the order of ``rrtmg_lw_[no]mcica_wrapper``."""
    return (ffi.from_buffer(fraction), ffi.from_buffer(tau), ffi.from_buffer(ice_path),
            ffi.from_buffer(liq_path), ffi.from_buffer(ice_eff_size),
            ffi.from_buffer(droplet_eff_radius))


def _nomcica_kernel(num_cols, num_layers, flags, layer_pressure, interface_pressure, layer_temp,
                    interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr,
                    o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity,
                    cloud_fraction, aerosol_tau, upward_longwave_flux, downward_longwave_flux,
                    longwave_heating_rate, up_lw_flux_clearsky, down_lw_flux_clearsky,
                    lw_heating_rate_clearsky, cloud_tau, cloud_ice_path, cloud_liq_path,
                    cloud_ice_eff_size, cloud_droplet_eff_radius, duflx_dt, duflxc_dt):
    common = (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp,
              h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, cfc11_vmr, cfc12_vmr,
              cfc22_vmr, ccl4_vmr, surf_emissivity, cloud_fraction, aerosol_tau,
              upward_longwave_flux, downward_longwave_flux, longwave_heating_rate,
              up_lw_flux_clearsky, down_lw_flux_clearsky, lw_heating_rate_clearsky, cloud_tau,
              cloud_ice_path, cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius)
    if not _shapes_ok(num_cols, num_layers, flags, *common, duflx_dt, duflxc_dt):
        return False
    sizes = np.empty(2, np.int32)
    sizes[0] = num_cols
    sizes[1] = num_layers
    before, after = _flux_pointers(sizes, flags, *common, duflx_dt, duflxc_dt)
    clouds = _cloud_pointers(cloud_fraction, cloud_tau, cloud_ice_path, cloud_liq_path,
                             cloud_ice_eff_size, cloud_droplet_eff_radius)
    _rrtmg_lw_nomcica(*before, *clouds, *after)
    keep_alive(sizes)
    return True


def _mcica_kernel(rrtm_iplon, num_cols, num_layers, flags, layer_pressure, interface_pressure,
                  layer_temp, interface_temp, surface_temp, h2o_vmr, o3_vmr, co2_vmr, ch4_vmr,
                  n2o_vmr, o2_vmr, cfc11_vmr, cfc12_vmr, cfc22_vmr, ccl4_vmr, surf_emissivity,
                  cloud_fraction, aerosol_tau, upward_longwave_flux, downward_longwave_flux,
                  longwave_heating_rate, up_lw_flux_clearsky, down_lw_flux_clearsky,
                  lw_heating_rate_clearsky, cloud_tau, cloud_ice_path, cloud_liq_path,
                  cloud_ice_eff_size, cloud_droplet_eff_radius, cloud_fraction_mcica,
                  cloud_tau_mcica, cloud_ice_path_mcica, cloud_liq_path_mcica,
                  cloud_ice_eff_size_mcica, cloud_droplet_eff_radius_mcica, duflx_dt, duflxc_dt):
    common = (layer_pressure, interface_pressure, layer_temp, interface_temp, surface_temp,
              h2o_vmr, o3_vmr, co2_vmr, ch4_vmr, n2o_vmr, o2_vmr, cfc11_vmr, cfc12_vmr,
              cfc22_vmr, ccl4_vmr, surf_emissivity, cloud_fraction, aerosol_tau,
              upward_longwave_flux, downward_longwave_flux, longwave_heating_rate,
              up_lw_flux_clearsky, down_lw_flux_clearsky, lw_heating_rate_clearsky, cloud_tau,
              cloud_ice_path, cloud_liq_path, cloud_ice_eff_size, cloud_droplet_eff_radius)
    g_points = (num_layers, num_cols, NUM_G_POINTS)
    layers = (num_layers, num_cols)
    if not (_shapes_ok(num_cols, num_layers, flags, *common, duflx_dt, duflxc_dt)
            and cloud_fraction_mcica.shape == g_points and cloud_tau_mcica.shape == g_points
            and cloud_ice_path_mcica.shape == g_points and cloud_liq_path_mcica.shape == g_points
            and cloud_ice_eff_size_mcica.shape == layers
            and cloud_droplet_eff_radius_mcica.shape == layers):
        return False
    sizes = np.empty(3, np.int32)
    sizes[0] = rrtm_iplon
    sizes[1] = num_cols
    sizes[2] = num_layers
    # The sub-column generator takes tau last.
    _mcica_subcol(
        ffi.from_buffer(sizes[0:]), ffi.from_buffer(sizes[1:]), ffi.from_buffer(sizes[2:]),
        ffi.from_buffer(flags[_CLOUD_OVERLAP_METHOD:]), ffi.from_buffer(flags[_PERMUTE_SEED:]),
        ffi.from_buffer(flags[_RANDOM_NUMBER_GENERATOR:]), ffi.from_buffer(layer_pressure),
        ffi.from_buffer(cloud_fraction), ffi.from_buffer(cloud_ice_path),
        ffi.from_buffer(cloud_liq_path), ffi.from_buffer(cloud_ice_eff_size),
        ffi.from_buffer(cloud_droplet_eff_radius), ffi.from_buffer(cloud_tau),
        ffi.from_buffer(cloud_fraction_mcica), ffi.from_buffer(cloud_ice_path_mcica),
        ffi.from_buffer(cloud_liq_path_mcica), ffi.from_buffer(cloud_ice_eff_size_mcica),
        ffi.from_buffer(cloud_droplet_eff_radius_mcica), ffi.from_buffer(cloud_tau_mcica))
    before, after = _flux_pointers(sizes[1:], flags, *common, duflx_dt, duflxc_dt)
    clouds = _cloud_pointers(cloud_fraction_mcica, cloud_tau_mcica, cloud_ice_path_mcica,
                             cloud_liq_path_mcica, cloud_ice_eff_size_mcica,
                             cloud_droplet_eff_radius_mcica)
    _rrtmg_lw_mcica(*before, *clouds, *after)
    keep_alive(sizes)
    return True


announce_compiles(_nomcica_kernel, "the RRTMG longwave wrapper")
announce_compiles(_mcica_kernel, "the RRTMG longwave McICA wrapper")

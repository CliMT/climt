"""cffi build script: one compiled module per Fortran library.

Run by setuptools through ``cffi_modules`` in ``setup.py`` (API, out-of-line
mode). A build-time file only: it lives in ``fortran/`` next to
``_build_config.py``, outside the installed package. Each ``*_ffi`` object below compiles ``climt_fortran._<name>_cffi``,
a small extension linked against the static library ``_lib/<system>/lib*.a``.
The prototypes are those of the routines' ``bind(c)`` names; every argument is
passed by pointer, as Fortran expects.
"""

import os

import cffi
from _build_config import build_config  # fortran/, put on sys.path by setup.py

_CONFIG = build_config()


def _ffi(module_name, prototypes, static_library, extra_flags=()):
    """An FFI declaring ``prototypes`` that builds ``climt_fortran.<module_name>``
    linked against ``static_library`` from ``_lib/<system>/``. ``extra_flags``
    go to both the compiler and the linker (``-fopenmp`` for RRTMG, as the
    Cython extensions had)."""
    ffi = cffi.FFI()
    ffi.cdef(prototypes)
    ffi.set_source(
        "climt_fortran." + module_name,
        "#include <stdint.h>\n" + prototypes,
        libraries=_CONFIG["libraries"],
        library_dirs=_CONFIG["library_dirs"],
        include_dirs=_CONFIG["include_dirs"],
        extra_compile_args=_CONFIG["extra_compile_args"] + list(extra_flags),
        extra_link_args=([os.path.join(_CONFIG["lib_dir"], static_library)] + list(extra_flags)
                         + _CONFIG["extra_link_args"]),
    )
    return ffi


EMANUEL_PROTOTYPES = """
void init_emanuel_convection_fortran(
    int32_t *pbl, int32_t *least_conv_level,
    double *thresh_water_level,
    double *crit_temp, double *entrain_coeff,
    double *downdraft_frac_area,
    double *precip_frac_outside_cloud,
    double *rain_speed, double *snow_speed,
    double *rain_evap_coeff, double *snow_evap_coeff,
    double *mom_tran_coeff, double *max_neg_temp_pert,
    double *beta, double *alpha, double *damp_amp,
    double *Cpd, double *Cpv, double *Cl,
    double *gas_const_vapour, double *gas_const_air,
    double *lat_heat, double *grav, double *density_water,
    double *reference_mass_flux_timescale);

void emanuel_convection(
    double *temp, double *q, double *qs,
    double *u, double *v,
    double *pmid, double *pint,
    int32_t *nlevs,
    int32_t *max_conv_lev,
    int32_t *num_tracers,
    double *dt, int32_t *conv_state,
    double *dtemp, double *dq, double *du,
    double *dv, double *precip,
    double *downdraft_vel_scale,
    double *downdraft_temp_scale,
    double *downdraft_q_scale,
    double *cloud_base_mass_flux,
    double *cape, double *tracers,
    double *dtracers);
"""

emanuel_ffi = _ffi("_emanuel_convection_cffi", EMANUEL_PROTOTYPES, "libemanuel.a")

# climt-fortran 2.0: cffi instead of Cython — short design

Date: 2026-11-19. Size: small change. Branch: `feature/scm-build`.
Second of three projects; comes after `2026-11-19-numba-ports-design.md`.

## Intent and why

Remove Cython from climt entirely. After the numba ports, climt-fortran holds
RRTMG longwave, RRTMG shortwave and Emanuel convection. Each is reached today
through a Cython `.pyx` file whose `cdef extern` blocks declare Fortran
routines that already export C names with `bind(c)` and take every argument
by pointer. cffi declares the same routines with no extra language in the
build, and numba can call cffi functions directly, which keeps approach A of
the build step open.

## Goals

- climt-fortran 2.0.0: no Cython at build or run time. Its modules
  `_rrtmg_lw`, `_rrtmg_sw` and `_emanuel_convection` keep their names and the
  Python functions climt calls today, now written in plain Python over cffi
  calls (array allocation, Fortran-order layout and output reshaping move from
  the `.pyx` files into Python).
- `_simple_physics` and `_dcmip` removed (ported to numba in the first
  project, which must have landed in climt first).
- `INTERFACE_VERSION` becomes 2 in climt-fortran and in climt's
  `EXPECTED_INTERFACE_VERSION`; climt's requirement becomes
  `climt-fortran>=2,<3`, and `MODULES` in `climt/_core/fortran.py` drops the
  two removed modules. The existing mismatch message tells users to upgrade
  both.
- Same numbers as 1.0: the existing RRTMG and Emanuel tests pass unchanged.
- `cffi` becomes a run-time dependency of climt-fortran (its compiled module
  needs `_cffi_backend`).
- Release order: the numba ports land in climt first (climt, still pinned
  `<2`, then simply stops importing the two modules from 1.0); climt-fortran
  2.0.0 is released; then climt raises its pin and interface version.
- Wheels built by the existing cibuildwheel workflow for the same platforms
  and Python versions.

## Non-goals and anti-goals

- No change to the Fortran itself, beyond adding `bind(c)` to any entry point
  that lacks one.
- climt-fortran 1.x stays as released. climt's `<2` pin means current climt
  releases never pick up 2.0.
- No slower calls: per-call time through the new modules is no worse than
  the Cython version, measured at one column.

## Approach

cffi **API (out-of-line) mode**: a build script declares the C prototypes
(copied from the existing `cdef extern` blocks) and compiles a small
extension linked against the Fortran objects, as the Cython extensions are
linked today. Chosen over ABI mode (`ffi.dlopen` of a shared library) because
API mode checks the prototypes against a C compiler at build time and needs
no shared-library search path at run time. Alternatives: ctypes (no
build-time checking, slower calls); keeping Cython (what the user wants
gone).

Fortran module state set by the initialisation routines (`rrtmg_set_constants`,
`rrtmg_lw_ini_wrapper`, `init_emanuel_convection_fortran`) stays as it is:
set once per process, as now.

## Left to the builder

The build script layout under `fortran/`; whether `.pyx` helper logic is
ported to one Python module per library or a shared helper; release notes
wording.

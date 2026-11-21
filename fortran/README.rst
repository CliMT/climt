=============
climt-fortran
=============

Compiled Fortran extensions for `climt <https://github.com/CliMT/climt>`_.

This package holds the three compiled modules that climt's Fortran-backed
components use: ``_emanuel_convection`` (``EmanuelConvection``),
``_rrtmg_lw`` (``RRTMGLongwave``) and ``_rrtmg_sw`` (``RRTMGShortwave``).
It has no user-facing API of its own: climt finds and loads these modules
itself.

Installing
----------

Most users never install this package directly. ``pip install climt`` pulls it
in as a dependency on every platform that has pre-built wheels. Without it,
climt still imports and every pure-Python component works; the components
above raise ``climt.FortranUnavailableError`` explaining how to install it.

``climt.has_fortran_extensions()`` tells you whether the extensions loaded.

The first time climt uses Emanuel convection or RRTMG after climt-fortran is
installed or upgraded, numba compiles the small wrappers that call the
Fortran code, which takes a few seconds; a line such as ``climt-fortran:
compiling the RRTMG longwave wrapper (this happens once after install; later
runs load it from cache)`` on stderr says so. numba caches the result, so
later runs start straight away.

Building from source
--------------------

You need ``gfortran``, a C compiler and ``make``. From a climt checkout::

    pip install ./fortran

or, for development, an editable install that builds the extensions in place::

    pip install -e ./fortran

``setup.py`` first builds the static Fortran libraries in ``_lib/`` with
``make``, then builds the cffi modules (``_cffi_build.py``) against them; no
Cython is needed.

To run the tests of a compiled install, use
``CLIMT_REQUIRE_FORTRAN=1 pytest fortran/tests``. The speed gates against the
1.0 Cython timings are skipped by default, because those timings were recorded
on one machine; maintainers run them with ``CLIMT_SPEED_GATES=1`` on their own
machine before tagging a release.
The compilers can be chosen with ``FC`` and ``CC`` and the optimisation flags
with ``CLIMT_OPT_FLAGS`` (default ``-O3``). On macOS the build looks for a
Homebrew ``gcc`` if ``FC`` and ``CC`` are not set; with Homebrew gcc you may
also need ``export SDKROOT="$(xcrun --show-sdk-path)"``.

The stub
--------

Setting ``CLIMT_FORTRAN_STUB=1`` builds a stub instead: a pure-Python
``py3-none-any`` wheel containing the Python wrappers in ``climt_fortran/``
but none of the compiled ``_*_cffi`` modules. It needs no compilers, and
needs neither cffi nor numpy at run time (the build still installs cffi, a
build requirement)::

    CLIMT_FORTRAN_STUB=1 pip install ./fortran

It satisfies climt's dependency on platforms where the extensions cannot be
built. climt reports it as a stub, and Fortran-backed components raise
``climt.FortranUnavailableError`` explaining how to build from source.

Licences
--------

The package itself is distributed under climt's BSD licence. The Fortran
sources under ``_lib/`` keep their original licences and copyright notices,
which are reproduced in the source files themselves (for example, RRTMG is
copyright Atmospheric & Environmental Research, Inc.).

Changes in 2.0.0 (unreleased)
-----------------------------

* Cython is gone: the modules are built with cffi, and building from source
  no longer needs Cython. ``INTERFACE_VERSION`` is now 2.
* ``_simple_physics`` and ``_dcmip`` were removed; climt now has numba
  versions of ``SimplePhysics`` and ``DcmipInitialConditions``.
* numba is a run-time dependency. The wrappers call Fortran from numba
  kernels, compiled once per install with numba's bundled LLVM (no system
  compiler needed) and cached; the first compile is announced by a one-line
  notice on stderr.
* RRTMG longwave and shortwave calls now release the GIL; Emanuel convection
  calls still hold it. climt-fortran is not thread-safe: do not call one
  module from two threads at once.
* Emanuel ``convect`` now accepts arrays in any memory layout (1.0 required
  C-contiguous arrays).
* RRTMG shortwave now raises an error if called before
  ``initialise_rrtm_radiation``.
* Solar variability arrays passed to the shortwave module are now validated.
* Strided ``facular_sunspot_ampl`` input now works (1.0 misread it).

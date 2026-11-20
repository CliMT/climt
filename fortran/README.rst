=============
climt-fortran
=============

Compiled Fortran extensions for `climt <https://github.com/CliMT/climt>`_.

This package holds the five compiled extension modules that climt's
Fortran-backed components use: ``_dcmip`` (``DcmipInitialConditions``),
``_emanuel_convection`` (``EmanuelConvection``), ``_rrtmg_lw``
(``RRTMGLongwave``), ``_rrtmg_sw`` (``RRTMGShortwave``) and
``_simple_physics`` (``SimplePhysics``). It has no user-facing API of its own:
climt finds and loads these modules itself.

Installing
----------

Most users never install this package directly. ``pip install climt`` pulls it
in as a dependency on every platform that has pre-built wheels. Without it,
climt still imports and every pure-Python component works; the components
above raise ``climt.FortranUnavailableError`` explaining how to install it.

``climt.has_fortran_extensions()`` tells you whether the extensions loaded.

Building from source
--------------------

You need ``gfortran``, a C compiler and ``make``. From a climt checkout::

    pip install ./fortran

or, for development, an editable install that builds the extensions in place::

    pip install -e ./fortran

``setup.py`` first builds the static Fortran libraries in ``_lib/`` with
``make``, then compiles the Cython wrappers in ``climt_fortran/`` against them.
The compilers can be chosen with ``FC`` and ``CC`` and the optimisation flags
with ``CLIMT_OPT_FLAGS`` (default ``-O3``). On macOS the build looks for a
Homebrew ``gcc`` if ``FC`` and ``CC`` are not set; with Homebrew gcc you may
also need ``export SDKROOT="$(xcrun --show-sdk-path)"``.

The stub
--------

Setting ``CLIMT_FORTRAN_STUB=1`` builds a stub instead: a pure-Python
``py3-none-any`` wheel containing only ``climt_fortran/__init__.py``, which
needs no compilers, Cython or numpy::

    CLIMT_FORTRAN_STUB=1 pip install ./fortran

It satisfies climt's dependency on platforms where the extensions cannot be
built. climt treats it like a missing package: Fortran-backed components raise
``climt.FortranUnavailableError``.

Licences
--------

The package itself is distributed under climt's BSD licence. The Fortran
sources under ``_lib/`` keep their original licences and copyright notices,
which are reproduced in the source files themselves (for example, RRTMG is
copyright Atmospheric & Environmental Research, Inc.).

Changes in 2.0.0 (unreleased)
-----------------------------

* climt-fortran now requires numba at run time. Emanuel convection calls
  Fortran through cffi from a numba kernel, so its first ``convect`` call
  JIT-compiles: about 2.4 s once per machine, about 0.16 s in later processes
  (from numba's cache).

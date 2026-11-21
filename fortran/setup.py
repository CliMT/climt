#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build script for climt-fortran, the compiled extensions used by climt.

The compiled modules are cffi modules (``_cffi_build.py``, next to this file)
linked against the static Fortran libraries built in ``_lib/``.

Set ``CLIMT_FORTRAN_STUB=1`` to build a stub: a pure-Python wheel with no
compiled extensions, which needs no compilers, cffi or numpy.
"""

import atexit
import os
import shutil
import sys
import tempfile

from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
STUB = os.environ.get("CLIMT_FORTRAN_STUB") == "1"

with open(os.path.join(HERE, "README.rst")) as readme_file:
    readme = readme_file.read()

metadata = dict(
    name="climt-fortran",
    version="2.0.0",
    description="Compiled Fortran extensions for climt",
    long_description=readme,
    long_description_content_type="text/x-rst",
    author="Rodrigo Caballero",
    author_email="rodrigo.caballero@misu.su.se",
    url="https://github.com/CliMT/climt",
    packages=["climt_fortran"],
    python_requires=">=3.10",
    license="BSD license",
    zip_safe=False,
    keywords="climt",
)


def fresh_build_base():
    """A new, empty build directory, removed when this process exits.

    The stub must not pick up products of an earlier compiled build: setuptools
    copies the whole of ``build/lib`` (and ``build/bdist.*``) into the wheel, so
    a stale ``.so`` left there would end up in a wheel tagged ``py3-none-any``.
    """
    path = tempfile.mkdtemp(prefix="climt-fortran-stub-build-")
    atexit.register(shutil.rmtree, path, True)
    return path


def compiled_build():
    """setup() arguments for the real build: the cffi modules, and the
    commands that build _lib first."""
    from setuptools.command.build_ext import build_ext as native_build_ext

    try:
        from setuptools.command.bdist_wheel import bdist_wheel as native_bdist_wheel
    except ImportError:
        from wheel.bdist_wheel import bdist_wheel as native_bdist_wheel

    from _build_config import build_static_libraries

    class climt_build_ext(native_build_ext):
        def run(self):
            build_static_libraries()
            native_build_ext.run(self)

    class climt_bdist_wheel(native_bdist_wheel):
        def run(self):
            self.run_command("build")
            native_bdist_wheel.run(self)

    return dict(
        # numba: the wrappers call Fortran from njit kernels (speed).
        install_requires=["numpy>=1.16.0", "cffi>=1.17", "numba>=0.59.0"],
        # Each entry builds climt_fortran._<name>_cffi. cffi wraps the
        # build_ext below, so _lib is still built first.
        cffi_modules=["_cffi_build.py:emanuel_ffi", "_cffi_build.py:rrtmg_lw_ffi",
                      "_cffi_build.py:rrtmg_sw_ffi"],
        cmdclass={
            "build_ext": climt_build_ext,
            "bdist_wheel": climt_bdist_wheel,
        },
    )


if STUB:
    # No cffi, no numpy, no extensions, no custom commands. Only Python
    # modules of the package go into the wheel: package data is off, so
    # in-place .so files next to __init__.py are not collected.
    setup(
        **metadata,
        install_requires=[],
        include_package_data=False,
        options={"build": {"build_base": fresh_build_base()}},
    )
else:
    setup(**metadata, **compiled_build())

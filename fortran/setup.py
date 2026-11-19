#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build script for climt-fortran, the compiled extensions used by climt.

Set ``CLIMT_FORTRAN_STUB=1`` to build a stub: a pure-Python wheel with no
compiled extensions, which needs no compilers, Cython or numpy.
"""

import atexit
import glob
import os
import platform
import re
import shutil
import subprocess
import tempfile

from setuptools import Extension, setup

HERE = os.path.dirname(os.path.abspath(__file__))
STUB = os.environ.get("CLIMT_FORTRAN_STUB") == "1"

with open(os.path.join(HERE, "README.rst")) as readme_file:
    readme = readme_file.read()

metadata = dict(
    name="climt-fortran",
    version="1.0.0",
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


# Find first gcc directory
def find_homebrew_gcc():
    # Check both standard Homebrew locations
    paths = ["/opt/homebrew/Cellar/gcc*", "/usr/local/Cellar/gcc*"]
    candidates = []
    for path in paths:
        candidates.extend(glob.glob(path))

    if candidates:
        # Prefer the newest version if multiple are found
        candidates.sort(reverse=True)
        print(f"Found Homebrew GCC candidates: {candidates}")
        return candidates[0]

    print("Warning: Could not find Homebrew GCC in standard locations.")
    return None


# Platform specific settings
def guess_compiler_name(env_name):

    search_string = ""
    if env_name == "FC":
        search_string = r"gfortran-\d+$"
    if env_name == "CC":
        search_string = r"gcc-\d+$"

    gcc_dir = find_homebrew_gcc()
    if gcc_dir:
        # Search deeply in the found gcc directory for bin/
        for root, dirs, files in os.walk(gcc_dir):
            if not root.endswith("/bin"):  # Optimization: Only look in bin directories
                continue

            for line in files:
                if re.match(search_string, line):
                    print(f"Guessing {env_name} = {line} from {root}")
                    os.environ[env_name] = os.path.join(root, line)
                    return  # Stop after finding one


def compiled_build():
    """setup() arguments for the real build: the five extensions and the commands that build _lib first."""
    import numpy as np
    from Cython.Build.Distutils import build_ext as native_build_ext

    try:
        from setuptools.command.bdist_wheel import bdist_wheel as native_bdist_wheel
    except ImportError:
        from wheel.bdist_wheel import bdist_wheel as native_bdist_wheel

    operating_system = platform.system()

    libraries = ["m", "gfortran"]
    default_link_args = []
    default_compile_args = []

    if operating_system == "Linux":
        libraries = ["m", "gfortran", "rt"]
        default_link_args = ["-lgfortran", "-lm"]

    compiled_path = os.path.join(HERE, "_lib")
    lib_path_list = [os.path.join(compiled_path, operating_system + "/")]
    inc_path = os.path.join(compiled_path, "include")
    include_dirs = [np.get_include(), inc_path]

    # Attempt to guess compilers on Darwin if not set
    if operating_system == "Darwin":
        if "FC" not in os.environ:
            guess_compiler_name("FC")
        if "CC" not in os.environ:
            guess_compiler_name("CC")

    if "FC" not in os.environ:
        os.environ["FC"] = "gfortran"
        os.environ["F77"] = "gfortran"

    if "CC" not in os.environ:
        os.environ["CC"] = "gcc"

    if "CLIMT_OPT_FLAGS" not in os.environ:
        os.environ["CLIMT_OPT_FLAGS"] = "-O3"

    if operating_system == "Windows":
        os.environ["CC"] = "gcc.exe"
        os.environ["FC"] = "gfortran.exe"
        os.environ["AR"] = "gcc-ar.exe"
        libraries = []
        default_link_args = ["-l:libgfortran.a", "-l:libquadmath.a", "-l:libm.a"]
        default_compile_args = ["-DMS_WIN64"]

    os.environ["FFLAGS"] = "-fPIC -fno-range-check " + os.environ["CLIMT_OPT_FLAGS"]
    os.environ["CFLAGS"] = "-fPIC " + os.environ["CLIMT_OPT_FLAGS"]

    if operating_system == "Darwin":
        gcc_dir = find_homebrew_gcc()
        if gcc_dir:
            # Find libgfortran.a so the extensions can link against it
            for root, dirs, files in os.walk(gcc_dir):
                if "lib" not in root:
                    continue  # Optimization
                for line in files:
                    if re.match(r"libgfortran\.a", line):
                        # Avoid i386 libs on 64bit systems
                        if "i386" not in root:
                            if root not in lib_path_list:
                                lib_path_list.append(root)

        if "MACOSX_DEPLOYMENT_TARGET" not in os.environ:
            os.environ["FFLAGS"] += " -mmacosx-version-min=10.9"
            os.environ["CFLAGS"] += " -mmacosx-version-min=10.9"
        default_link_args = []
        os.environ["LDSHARED"] = os.environ["CC"] + " -bundle -undefined dynamic_lookup"

    print("Compilers: ", os.environ.get("CC"), os.environ.get("FC"))
    print("LDSHARED: ", os.environ.get("LDSHARED"))
    print("FFLAGS: ", os.environ.get("FFLAGS"))
    print("Lib Paths: ", lib_path_list)

    def build_libraries():
        """Build the static Fortran libraries in _lib with make."""
        curr_dir = os.getcwd()
        os.chdir(compiled_path)
        os.environ["PWD"] = compiled_path

        cmd = ["make", "CLIMT_ARCH=" + operating_system]
        print(f"Building libraries in {compiled_path}")
        print(f"Command: {cmd}")

        ret = subprocess.call(cmd)

        os.chdir(curr_dir)
        os.environ["PWD"] = curr_dir

        if ret != 0:
            raise RuntimeError(f"Library build failed with exit code {ret}")

    class climt_build_ext(native_build_ext):
        def run(self):
            build_libraries()
            native_build_ext.run(self)

    class climt_bdist_wheel(native_bdist_wheel):
        def run(self):
            self.run_command("build")
            native_bdist_wheel.run(self)

    def fortran_extension(name, library, openmp=False):
        """An extension climt_fortran.<name> built from climt_fortran/<name>.pyx and _lib/<system>/<library>."""
        openmp_flags = ["-fopenmp"] if openmp else []
        return Extension(
            "climt_fortran." + name,
            sources=["climt_fortran/" + name + ".pyx"],
            libraries=libraries,
            include_dirs=include_dirs,
            extra_compile_args=default_compile_args + openmp_flags,
            library_dirs=lib_path_list,
            extra_link_args=[os.path.join(lib_path_list[0], library)] + openmp_flags + default_link_args,
        )

    ext_modules = [
        fortran_extension("_simple_physics", "libsimple_physics.a"),
        fortran_extension("_emanuel_convection", "libemanuel.a"),
        fortran_extension("_rrtmg_lw", "librrtmg_lw.a", openmp=True),
        fortran_extension("_rrtmg_sw", "librrtmg_sw.a", openmp=True),
        fortran_extension("_dcmip", "libdcmip.a"),
    ]

    return dict(
        install_requires=["numpy>=1.16.0"],
        ext_modules=ext_modules,
        cmdclass={
            "build_ext": climt_build_ext,
            "bdist_wheel": climt_bdist_wheel,
        },
    )


if STUB:
    # No Cython, no numpy, no extensions, no custom commands. Only Python
    # modules of the package go into the wheel: package data is off, so
    # in-place .so/.c/.pyx files next to __init__.py are not collected.
    setup(
        **metadata,
        install_requires=[],
        include_package_data=False,
        options={"build": {"build_base": fresh_build_base()}},
    )
else:
    setup(**metadata, **compiled_build())

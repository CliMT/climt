"""Platform and compiler settings for building climt-fortran.

Shared by ``setup.py`` and the extension build scripts: ``build_config()`` sets
the compiler environment (FC, CC, FFLAGS, CFLAGS, LDSHARED) and returns the
link and include settings; ``build_static_libraries()`` builds ``_lib`` with make.
"""

import glob
import os
import platform
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
COMPILED_PATH = os.path.join(HERE, "_lib")


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


def build_config():
    """Set the compiler environment and return the build settings.

    Returns a dict with ``libraries``, ``library_dirs``, ``extra_link_args``,
    ``extra_compile_args``, ``include_dirs`` (``_lib/include``; no NumPy) and
    ``lib_dir`` (``_lib/<system>/``, where the static libraries are).
    """
    operating_system = platform.system()

    libraries = ["m", "gfortran"]
    default_link_args = []
    default_compile_args = []

    if operating_system == "Linux":
        libraries = ["m", "gfortran", "rt"]
        default_link_args = ["-lgfortran", "-lm"]

    lib_path_list = [os.path.join(COMPILED_PATH, operating_system + "/")]
    include_dirs = [os.path.join(COMPILED_PATH, "include")]

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

    return dict(
        libraries=libraries,
        library_dirs=lib_path_list,
        extra_link_args=default_link_args,
        extra_compile_args=default_compile_args,
        include_dirs=include_dirs,
        lib_dir=lib_path_list[0],
    )


def build_static_libraries():
    """Build the static Fortran libraries in _lib with make."""
    curr_dir = os.getcwd()
    os.chdir(COMPILED_PATH)
    os.environ["PWD"] = COMPILED_PATH

    cmd = ["make", "CLIMT_ARCH=" + platform.system()]
    print(f"Building libraries in {COMPILED_PATH}")
    print(f"Command: {cmd}")

    ret = subprocess.call(cmd)

    os.chdir(curr_dir)
    os.environ["PWD"] = curr_dir

    if ret != 0:
        raise RuntimeError(f"Library build failed with exit code {ret}")

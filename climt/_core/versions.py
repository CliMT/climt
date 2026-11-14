"""Version report for bug reports, built without importing reported packages."""
import platform
import sys
from importlib.metadata import PackageNotFoundError, version

REPORTED_PACKAGES = (
    "sympl", "numpy", "scipy", "xarray", "unyt", "pint", "numba", "jax", "jaxlib",
    "netCDF4", "zarr", "pyzmq", "panel", "holoviews", "gfs_dynamical_core",
)


def _installed_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def show_versions(file=None):
    """Print and return versions of climt, Python, the platform and key dependencies.

    Package versions come from installed metadata, so reported packages are not imported.
    Missing packages map to None and print as "not installed".
    """
    import climt

    # A broken compiled extension can raise more than ImportError (e.g. the
    # ValueError of a numpy ABI mismatch, climt#174); report it, don't raise it.
    fortran_error = None
    try:
        fortran = bool(climt.has_fortran_extensions())
    except Exception as err:
        fortran, fortran_error = False, err

    info = {
        "climt": climt.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "fortran_extensions": fortran,
    }
    for name in REPORTED_PACKAGES:
        info[name] = _installed_version(name)

    if file is None:
        file = sys.stdout
    for key, value in info.items():
        text = "not installed" if value is None else str(value)
        if key == "fortran_extensions" and fortran_error is not None:
            text += " ({!r})".format(fortran_error)
        print("{}: {}".format(key, text), file=file)
    return info

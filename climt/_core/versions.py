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

    info = {
        "climt": climt.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "fortran_extensions": bool(climt.has_fortran_extensions()),
    }
    for name in REPORTED_PACKAGES:
        info[name] = _installed_version(name)

    if file is None:
        file = sys.stdout
    for key, value in info.items():
        print("{}: {}".format(key, "not installed" if value is None else value), file=file)
    return info

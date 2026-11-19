"""Version report for bug reports, built without importing reported packages."""
import platform
import sys
from importlib.metadata import PackageNotFoundError, version

REPORTED_PACKAGES = (
    "sympl", "numpy", "scipy", "xarray", "unyt", "pint", "numba", "jax", "jaxlib",
    "netCDF4", "zarr", "pyzmq", "panel", "holoviews", "gfs_dynamical_core",
    "climt-fortran",
)


def _installed_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _fortran_note():
    """Return "<state>: <detail>" for the climt-fortran companion, "" when it is
    available, or a description of the error if its status cannot be read."""
    try:
        from .fortran import fortran_status
        status = fortran_status()
    except Exception as err:
        return "status unavailable: {!r}".format(err)
    if status["state"] == "available":
        return ""
    if status["detail"]:
        return "{}: {}".format(status["state"], status["detail"])
    return status["state"]


def show_versions(file=None):
    """Print and return versions of climt, Python, the platform and key dependencies.

    Package versions come from installed metadata, so reported packages are not imported.
    Missing packages map to None and print as "not installed".
    """
    import climt

    # A broken compiled extension can raise more than ImportError (e.g. the
    # ValueError of a numpy ABI mismatch, climt#174). The loader records it in
    # the companion's status; report that here, never raise it.
    fortran_error = None
    try:
        fortran = bool(climt.has_fortran_extensions())
    except Exception as err:
        fortran, fortran_error = False, err
    fortran_note = _fortran_note()

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
        if key == "fortran_extensions":
            if fortran_error is not None:
                text += " ({!r})".format(fortran_error)
            elif fortran_note:
                text += " ({})".format(fortran_note)
        print("{}: {}".format(key, text), file=file)
    return info

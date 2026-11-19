"""Loader for the optional ``climt_fortran`` companion package.

The compiled Fortran/Cython extensions live in the separate ``climt-fortran``
distribution. This module is the only place the core imports them from.
"""
import importlib
import importlib.machinery

EXPECTED_INTERFACE_VERSION = 1
MODULES = ("_dcmip", "_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw", "_simple_physics")

_PACKAGE = "climt_fortran"
_cache = None


class FortranUnavailableError(ImportError):
    """Raised when a component needs a Fortran extension that cannot be loaded."""


def _reset_cache() -> None:
    """Forget cached import results. For tests."""
    global _cache
    _cache = None


def _load():
    """Import the companion and its modules once; return the cached result."""
    global _cache
    if _cache is not None:
        return _cache
    result = {
        "state": "missing", "installed": False, "version": None,
        "interface_version": None, "detail": "",
        "modules": {name: False for name in MODULES}, "loaded": {},
    }
    try:
        package = importlib.import_module(_PACKAGE)
    except ImportError:
        _cache = result
        return result
    except Exception as error:
        result["installed"] = True
        result["state"] = "broken"
        result["detail"] = "%s: %s: %s" % (_PACKAGE, type(error).__name__, error)
        _cache = result
        return result
    result["installed"] = True
    result["version"] = getattr(package, "__version__", None)
    interface = getattr(package, "INTERFACE_VERSION", None)
    result["interface_version"] = interface
    if interface != EXPECTED_INTERFACE_VERSION:
        result["state"] = "interface_mismatch"
        result["detail"] = "%s interface version %r, climt expects %d" % (
            _PACKAGE, interface, EXPECTED_INTERFACE_VERSION)
        _cache = result
        return result
    # Look for the extensions only inside the package just imported, whose
    # interface version was checked. importlib.util.find_spec would also ask
    # every finder on sys.meta_path, and setuptools' editable-install finder
    # resolves climt_fortran.<name> to its checkout whatever __path__ says.
    search_path = getattr(package, "__path__", None)
    failures = []
    found_any = False
    for name in MODULES:
        full = "%s.%s" % (_PACKAGE, name)
        try:
            spec = None if search_path is None else importlib.machinery.PathFinder.find_spec(full, list(search_path))
        except Exception as error:
            failures.append("%s: %s" % (name, error))
            continue
        if spec is None:
            failures.append("%s: module not found" % name)
            continue
        found_any = True
        try:
            result["loaded"][name] = importlib.import_module(full)
            result["modules"][name] = True
        except Exception as error:
            failures.append("%s: %s: %s" % (name, type(error).__name__, error))
    if not found_any:
        result["state"] = "stub"
        result["detail"] = "%s contains no compiled extensions" % _PACKAGE
    elif failures:
        result["state"] = "broken"
        result["detail"] = "; ".join(failures)
    else:
        result["state"] = "available"
    _cache = result
    return result


def fortran_status() -> dict:
    """Describe the companion: keys 'state', 'installed', 'version',
    'interface_version', 'modules' and 'detail'."""
    result = _load()
    status = {key: value for key, value in result.items() if key != "loaded"}
    status["modules"] = dict(result["modules"])
    return status


def load_extension(name: str):
    """Return the imported extension module ``name``, or None if unavailable."""
    if name not in MODULES:
        raise ValueError("unknown Fortran extension %r; expected one of %s" % (name, ", ".join(MODULES)))
    return _load()["loaded"].get(name)


def unavailable_error(component: str, alternative: str | None = None) -> FortranUnavailableError:
    """Build (do not raise) the error for a component whose extension is missing."""
    result = _load()
    state = result["state"]
    if state == "missing":
        body = "The climt-fortran package is not installed. Install it with: pip install climt-fortran"
    elif state == "stub":
        body = ("The installed climt-fortran package contains no compiled extensions. "
                "Build it from source (needs gfortran and a C compiler) with: "
                "pip install --no-binary climt-fortran climt-fortran")
    elif state == "interface_mismatch":
        body = ("The installed climt-fortran package does not match this version of climt (%s). "
                "Upgrade both with: pip install -U climt climt-fortran" % result["detail"])
    else:
        body = ("The installed climt-fortran package failed to load: %s. "
                "Reinstall it with: pip install --force-reinstall climt-fortran" % result["detail"])
    message = "%s needs a compiled Fortran extension that is unavailable. %s" % (component, body)
    if alternative is not None:
        message += " Use %s instead." % alternative
    return FortranUnavailableError(message)

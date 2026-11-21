"""Loader for the optional ``climt_fortran`` companion package.

The compiled Fortran extensions live in the separate ``climt-fortran``
distribution. This module is the only place the core imports them from.

Each module ``climt_fortran._<name>`` counts as installed when its compiled
code is present: the cffi module ``_<name>_cffi`` that a Python wrapper
``_<name>`` calls, or ``_<name>`` itself as an extension module. A wrapper
alone, as in a stub install, does not count.
"""
import importlib
import importlib.machinery

EXPECTED_INTERFACE_VERSION = 2
MODULES = ("_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw")

_PACKAGE = "climt_fortran"
_cache = None


class FortranUnavailableError(ImportError):
    """Raised when a component needs a Fortran extension that cannot be loaded."""


def _reset_cache() -> None:
    """Forget cached import results. For tests."""
    global _cache
    _cache = None


def _compiled_code_found(full, search_path) -> bool:
    """Whether the compiled code behind module ``full`` is in ``search_path``."""
    # Keep in step with _installed in fortran/climt_fortran/__init__.py.
    finder = importlib.machinery.PathFinder
    if finder.find_spec(full + "_cffi", search_path) is not None:
        return True
    spec = finder.find_spec(full, search_path)
    return spec is not None and str(spec.origin).endswith(tuple(importlib.machinery.EXTENSION_SUFFIXES))


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
            found = search_path is not None and _compiled_code_found(full, list(search_path))
        except Exception as error:
            failures.append("%s: %s" % (name, error))
            continue
        if not found:
            failures.append("%s: compiled module not found" % name)
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
    """Build (do not raise) the error for a component whose extension is missing.

    The message states why the extension is unavailable, gives the command
    that fixes it on its own indented line and, if ``alternative`` is given,
    ends with "Use <alternative> instead."
    """
    result = _load()
    state = result["state"]
    detail = result["detail"].rstrip(". ")
    if state == "missing":
        reason = "The climt-fortran package is not installed. Install it with:"
        command = "pip install climt-fortran"
    elif state == "stub":
        reason = ("The installed climt-fortran package contains no compiled extensions. "
                  "Build it from source (needs gfortran and a C compiler) with:")
        command = "pip install --force-reinstall --no-deps --no-binary climt-fortran climt-fortran"
    elif state == "interface_mismatch":
        reason = ("The installed climt-fortran package does not match this version of climt (%s). "
                  "Upgrade both with:" % detail)
        command = "pip install -U climt climt-fortran"
    else:
        reason = "The installed climt-fortran package failed to load%s. Reinstall it with:" % (
            " (%s)" % detail if detail else "")
        command = "pip install --force-reinstall climt-fortran"
    lines = ["%s needs a compiled Fortran extension that is unavailable. %s" % (component, reason),
             "", "    " + command]
    if alternative is not None:
        lines += ["", "Use %s instead." % alternative]
    return FortranUnavailableError("\n".join(lines))

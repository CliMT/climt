"""Compiled Fortran extensions for climt.

climt loads these modules itself (see ``climt._core.fortran``); importing this
package does not import any of them.

Each module ``_<name>`` is a Python wrapper over a compiled cffi module
``_<name>_cffi``.
"""
import importlib.machinery

__version__ = "2.0.0"
INTERFACE_VERSION = 2
MODULES = ("_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw")


def _installed(name: str) -> bool:
    """Whether the compiled code behind module ``name`` is present: its cffi
    module. A wrapper alone (a stub install) does not count."""
    # Keep in step with _compiled_code_found in climt/_core/fortran.py.
    # Search only this package's directory: importlib.util.find_spec would
    # also ask setuptools' editable-install finder, which maps any
    # climt_fortran.<name> to a source checkout.
    full = __name__ + "." + name + "_cffi"
    return importlib.machinery.PathFinder.find_spec(full, __path__) is not None


def available() -> dict:
    """Map each module name to whether its compiled code is installed, without importing it."""
    return {name: _installed(name) for name in MODULES}

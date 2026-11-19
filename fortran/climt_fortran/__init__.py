"""Compiled Fortran extensions for climt.

climt loads these modules itself (see ``climt._core.fortran``); importing this
package does not import any of them.

Each module ``_<name>`` is either a Python wrapper over a compiled cffi module
``_<name>_cffi``, or (until all have moved to cffi) a Cython extension.
"""
import importlib.machinery

__version__ = "1.0.0"
INTERFACE_VERSION = 1
MODULES = ("_dcmip", "_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw", "_simple_physics")


def _installed(name: str) -> bool:
    """Whether the compiled code behind module ``name`` is present: its cffi
    module, or the module itself as an extension. A wrapper alone (a stub
    install) does not count."""
    # Search only this package's directory: importlib.util.find_spec would
    # also ask setuptools' editable-install finder, which maps any
    # climt_fortran.<name> to a source checkout.
    full = __name__ + "." + name
    if importlib.machinery.PathFinder.find_spec(full + "_cffi", __path__) is not None:
        return True
    spec = importlib.machinery.PathFinder.find_spec(full, __path__)
    return spec is not None and str(spec.origin).endswith(tuple(importlib.machinery.EXTENSION_SUFFIXES))


def available() -> dict:
    """Map each module name to whether its compiled code is installed, without importing it."""
    return {name: _installed(name) for name in MODULES}

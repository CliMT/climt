"""Compiled Fortran extensions for climt.

climt loads these modules itself (see ``climt._core.fortran``); importing this
package does not import any of them.
"""
import importlib.util

__version__ = "1.0.0"
INTERFACE_VERSION = 1
MODULES = ("_dcmip", "_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw", "_simple_physics")


def available() -> dict:
    """Map each extension module name to whether it is installed, without importing it."""
    return {name: importlib.util.find_spec(__name__ + "." + name) is not None for name in MODULES}

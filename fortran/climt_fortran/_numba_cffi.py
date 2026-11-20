"""Call the Fortran routines of a compiled cffi module from cached numba code.

A cffi function called from ``@numba.njit`` code is compiled in as a constant
address, so numba will not cache that code ("uses dynamic globals") and
recompiles it in every process. Calling it instead by a symbol name that is
registered with LLVM at import time keeps the code cacheable.

numba's on-disk cache does not see changes to globals: if a cffi prototype in
``fortran/_cffi_build.py`` changes while the wrapper ``.py`` holding the njit
kernel does not, a stale cached kernel could call Fortran with the old
signature. The Fortran sources are frozen and a reinstall rewrites the ``.py``
(invalidating the cache), but after editing a prototype clear
``climt_fortran/__pycache__`` (and numba's user cache directory if
site-packages was read-only).
"""

import llvmlite.binding
from numba import types
from numba.core.typing import cffi_utils


def external_function(module, name):
    """A numba ``ExternalFunction`` calling ``module.lib.<name>`` from njit code.

    ``module`` is a compiled cffi module (``climt_fortran._<lib>_cffi``). Its
    prototype gives the signature, so it is declared in one place only. The
    module is also registered with numba, which types ``module.ffi.from_buffer``
    in njit code.
    """
    cffi_utils.register_module(module)
    function = getattr(module.lib, name)
    symbol = "%s_%s" % (module.__name__.replace(".", "_"), name)
    llvmlite.binding.add_symbol(symbol, cffi_utils.get_pointer(function))
    return types.ExternalFunction(symbol, cffi_utils.make_function_type(function).sig)

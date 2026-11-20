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
import numba
from numba import types
from numba.core.typing import cffi_utils

_KERNELS = {}


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


def compiled_kernel(function, signature):
    """``function`` compiled by numba for ``signature`` only, on first use.

    Cached on disk, with the GIL released. Called with arguments of other
    types (an array of another dtype, layout or ndim, a read-only array,
    None), the result raises TypeError before running anything, and nothing
    is recompiled. Later calls with the same function and signature return
    the same compiled function (an equal but distinct signature object
    compiles again, from numba's disk cache).
    """
    # Keyed on the signature object's id: hashing a numba signature of ~40
    # arguments takes ~6 us, as long as some calls. The entry holds the
    # signature, so its id is not reused while the entry exists.
    key = (function, id(signature))
    try:
        return _KERNELS[key][1]
    except KeyError:
        kernel = numba.njit(signature, cache=True, nogil=True)(function)
        _KERNELS[key] = (signature, kernel)
        return kernel


@numba.njit(cache=True)
def keep_alive(array):
    """Do nothing with ``array``, so that it is still in use here.

    numba frees an array created in a kernel right after its last use, and a
    pointer from ``ffi.from_buffer`` does not keep it alive: call this after a
    Fortran call given a pointer into such an array.
    """
    return array.size >= 0

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

Compiling a kernel takes from about 0.4 s to a few seconds. A function given
to ``announce_compiles`` has its first compile in each process announced on
stderr, so that the pause is explained; loading it from numba's cache prints
nothing. The listener behind this is registered with numba for the whole
process, so it runs at every numba compile, of any package; for functions
not given to ``announce_compiles`` its cost is a dict lookup.
"""

import sys

import llvmlite.binding
import numba
from numba import types
from numba.core.typing import cffi_utils

from ._arrays import contiguous_arguments, copy_back_outputs

_KERNELS = {}
ONCE_AFTER_INSTALL = "this happens once after install; later runs load it from cache"
# Python function -> its notice, for the functions given to announce_compiles.
_NOTICES = {}
_ANNOUNCED = set()


def announce_compiles(function, what, when=ONCE_AFTER_INSTALL):
    """Announce on stderr the first time in this process that numba compiles
    ``function`` (a Python function, or the numba dispatcher of one).

    The notice is one line, ``climt-fortran: compiling <what> (<when>)``,
    written and flushed before the compile starts. It is not printed when
    numba loads ``function`` from its cache, nor for later compiles in the
    same process (e.g. a lazily compiled function given new argument types).
    """
    function = getattr(function, "py_func", function)
    _NOTICES[function] = "climt-fortran: compiling %s (%s)\n" % (what, when)


def _announce(function):
    """Write ``function``'s notice to stderr unless it was written before.
    Never raises: a missing or broken stderr only loses the notice."""
    try:
        if function in _ANNOUNCED or function not in _NOTICES:
            return
        _ANNOUNCED.add(function)
        sys.stderr.write(_NOTICES[function])
        sys.stderr.flush()
    except Exception:
        pass


def _listen_for_compiles():
    """Have numba call ``_announce`` just before it compiles a function, and
    return True; return False if that is not possible.

    This is the cache-miss probe, and the only code here that depends on
    numba internals: numba's event API (``numba.core.event``). In numba 0.68
    ``Dispatcher.compile`` first looks for the signature in the dispatcher,
    then in the on-disk cache, and only if both miss fires a
    ``"numba:compile"`` event, with the dispatcher in ``event.data``, around
    the compile itself. tests/test_numba_cffi.py pins this.
    """
    try:
        from numba.core import event

        class CompileListener(event.Listener):
            def on_start(self, compile_event):
                try:
                    _announce(compile_event.data["dispatcher"].py_func)
                except Exception:
                    pass

            def on_end(self, compile_event):
                pass

        event.register("numba:compile", CompileListener())
        return True
    except Exception:
        return False


_COMPILES_DETECTED = _listen_for_compiles()


def compile_may_follow(function):
    """Call before a call or decoration that may compile ``function``. If
    compiles cannot be detected, announce one now (once per process): a
    redundant notice is better than an unexplained pause."""
    if not _COMPILES_DETECTED:
        _announce(getattr(function, "py_func", function))


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
    compiles again, from numba's disk cache). Give ``function`` to
    ``announce_compiles`` to have a compile announced.
    """
    # Keyed on the signature object's id: hashing a numba signature of ~40
    # arguments takes ~6 us, as long as some calls. The entry holds the
    # signature, so its id is not reused while the entry exists.
    key = (function, id(signature))
    try:
        return _KERNELS[key][1]
    except KeyError:
        compile_may_follow(function)
        kernel = numba.njit(signature, cache=True, nogil=True)(function)
        _KERNELS[key] = (signature, kernel)
        return kernel


def run_kernel(function, signature, leading, arrays, table, shapes, num_cols, num_layers):
    """Run ``function`` compiled for ``signature`` (see ``compiled_kernel``)
    on ``(*leading, *arrays)``; it returns whether it took the arguments.

    If it does not take them (TypeError from numba's dispatch, or False: a
    wrong shape), check every array against ``table`` (``(name, shape kind,
    output)`` for each, in order; ``shapes(num_cols, num_layers)`` maps each
    kind to its shape), raising with the name of the bad one, and run it on
    C-contiguous writeable copies, writing outputs back.
    """
    kernel = compiled_kernel(function, signature)
    try:
        if kernel(*leading, *arrays):
            return
    except TypeError:
        pass
    shape = shapes(num_cols, num_layers)
    copy_back = []
    arrays = contiguous_arguments(
        arrays, [(name, shape[kind], output) for name, kind, output in table], copy_back)
    if not kernel(*leading, *arrays):
        raise RuntimeError("the arguments passed the checks but not the kernel's")
    copy_back_outputs(copy_back)


@numba.njit(cache=True)
def keep_alive(array):
    """Do nothing with ``array``, so that it is still in use here.

    numba frees an array created in a kernel right after its last use, and a
    pointer from ``ffi.from_buffer`` does not keep it alive: call this after a
    Fortran call given a pointer into such an array.
    """
    return array.size >= 0

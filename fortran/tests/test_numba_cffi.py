"""The numba helpers: the compiled-kernel cache and the first-compile notice."""

import os
import subprocess
import sys
import types

import numba
import pytest

from climt_fortran import _numba_cffi

HERE = os.path.dirname(os.path.abspath(__file__))
NOTICE = "climt-fortran: compiling "


@pytest.fixture
def fresh_notices(monkeypatch):
    """No notice registered or written yet in this test."""
    monkeypatch.setattr(_numba_cffi, "_NOTICES", {})
    monkeypatch.setattr(_numba_cffi, "_ANNOUNCED", set())


def test_compiled_kernel_cache_key(monkeypatch):
    # Keyed on (function, id(signature)): the same function and signature
    # object give the same kernel, compiled once; another signature object
    # (equal or not) or another function compiles again.
    compiles = []

    def njit(signature, **options):
        assert options == {"cache": True, "nogil": True}
        return lambda function: compiles.append((function, signature)) or object()

    monkeypatch.setattr(_numba_cffi, "numba", types.SimpleNamespace(njit=njit))
    monkeypatch.setattr(_numba_cffi, "_KERNELS", {})

    def kernel(x):
        return x

    def other(x):
        return x

    signature = numba.types.int64(numba.types.int64)
    first = _numba_cffi.compiled_kernel(kernel, signature)
    assert _numba_cffi.compiled_kernel(kernel, signature) is first
    assert compiles == [(kernel, signature)]
    equal = numba.types.int64(numba.types.int64)
    assert equal == signature and equal is not signature
    assert _numba_cffi.compiled_kernel(kernel, equal) is not first
    different = numba.types.float64(numba.types.float64)
    assert _numba_cffi.compiled_kernel(kernel, different) is not first
    assert _numba_cffi.compiled_kernel(other, signature) is not first
    assert compiles[1:] == [(kernel, equal), (kernel, different), (other, signature)]
    assert _numba_cffi.compiled_kernel(kernel, different) is _numba_cffi.compiled_kernel(kernel, different)
    assert len(compiles) == 4


def test_compiles_are_detected():
    # The cache-miss probe installs with numba 0.68 (the fallback is not in use).
    assert _numba_cffi._COMPILES_DETECTED


def test_first_compile_announced_once(fresh_notices, capfd):
    def announced(x):
        return x + 1

    def silent(x):
        return x + 2

    _numba_cffi.announce_compiles(announced, "a test function", "just a test")
    function, other = numba.njit(announced), numba.njit(silent)
    assert function(1) == 2  # compiles: announced
    assert function(1.5) == 2.5  # compiles again, for float64: not announced again
    assert other(1) == 3  # not given to announce_compiles
    _numba_cffi.compile_may_follow(function)  # compiles are detected: nothing to do
    out, err = capfd.readouterr()
    assert out == "" and err == "climt-fortran: compiling a test function (just a test)\n"


def test_fallback_announces_before_first_call(fresh_notices, monkeypatch, capfd):
    # Without the probe, compile_may_follow announces a possible compile, once.
    def kernel(x):
        return x

    _numba_cffi.announce_compiles(kernel, "a test kernel")
    monkeypatch.setattr(_numba_cffi, "_COMPILES_DETECTED", False)
    _numba_cffi.compile_may_follow(kernel)
    _numba_cffi.compile_may_follow(numba.njit(kernel))  # a dispatcher counts as its function
    assert capfd.readouterr().err == NOTICE + "a test kernel (%s)\n" % _numba_cffi.ONCE_AFTER_INSTALL


def test_probe_failure_falls_back(monkeypatch):
    def broken(kind, listener):
        raise RuntimeError("event API changed")

    monkeypatch.setattr(numba.core.event, "register", broken)
    assert _numba_cffi._listen_for_compiles() is False


@pytest.mark.parametrize("stderr", [None, "raises"])
def test_notice_never_raises(fresh_notices, monkeypatch, stderr):
    class Broken:
        def write(self, text):
            raise OSError("stderr is closed")

    def announced(x):
        return x * 2

    _numba_cffi.announce_compiles(announced, "a test function")
    monkeypatch.setattr(sys, "stderr", None if stderr is None else Broken())
    assert numba.njit(announced)(3) == 6


# Calls every kernel twice, the Emanuel one also with Fortran-ordered arrays
# (a new specialisation): each wrapper may be announced only once.
_CALLS_SCRIPT = """
import sys
import numpy as np
sys.path.insert(0, %r)
import test_parity
for call in ("emanuel", "rrtmg_lw", "rrtmg_lw_mcica", "rrtmg_sw", "rrtmg_sw_mcica"):
    for _ in range(2):
        test_parity.run_call(call, test_parity.reference_inputs(call))
inputs = {name: np.asfortranarray(value) if isinstance(value, np.ndarray) else value
          for name, value in test_parity.reference_inputs("emanuel").items()}
test_parity.run_call("emanuel", inputs)
""" % HERE
EXPECTED_NOTICES = [
    NOTICE + "the Emanuel convection wrapper (this happens once after install, and again for each new "
    "combination of array layouts given to convect; later runs load it from cache)",
    NOTICE + "the RRTMG longwave wrapper (this happens once after install; later runs load it from cache)",
    NOTICE + "the RRTMG longwave McICA wrapper (this happens once after install; later runs load it from cache)",
    NOTICE + "the RRTMG shortwave wrapper (this happens once after install; later runs load it from cache)",
    NOTICE + "the RRTMG shortwave McICA wrapper (this happens once after install; later runs load it from cache)",
]


def _notices_from_a_process(cache_dir):
    """Notice lines a new process writes to stderr making ``_CALLS_SCRIPT``'s calls."""
    import importlib.util

    if importlib.util.find_spec("climt_fortran._rrtmg_lw_cffi") is None:
        pytest.skip("the cffi modules are not built")
    result = subprocess.run([sys.executable, "-c", _CALLS_SCRIPT], capture_output=True, text=True,
                            env=dict(os.environ, NUMBA_CACHE_DIR=str(cache_dir)), check=False)
    assert result.returncode == 0, result.stderr
    return [line for line in result.stderr.splitlines() if line.startswith(NOTICE)]


def test_notice_on_cold_cache_only(tmp_path):
    # A cold numba cache: each wrapper's compile is announced once. A second
    # process loads them from the same cache directory and prints nothing.
    assert _notices_from_a_process(tmp_path) == EXPECTED_NOTICES
    assert _notices_from_a_process(tmp_path) == []

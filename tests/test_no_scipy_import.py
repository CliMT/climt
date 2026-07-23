import subprocess
import sys


def test_import_climt_does_not_import_scipy():
    """``import climt`` and running the CORK demo component must not need scipy.

    Run in a fresh subprocess with ``scipy`` blocked via a ``sys.modules``
    sentinel (``sys.modules["scipy"] = None``), the same technique already
    used in ``tests/test_correlated_k_npz.py`` to force an ``ImportError`` on
    any attempted scipy use. This is required because the *dev* conda
    environment has both numba and scipy installed: numba's own BLAS/LAPACK
    bootstrap (``numba.np.arraymath._check_blas``) and xarray's optional
    scipy netCDF backend probe both opportunistically import scipy
    submodules when scipy happens to be importable, even though neither
    failure would occur in the actual Pyodide/pure-wheel target where scipy
    is never installed at all. Both of those libraries already handle
    ``ImportError`` on scipy gracefully, so blocking it here faithfully
    simulates the real deployment environment rather than merely testing
    "does this dev env happen to have scipy installed".

    The assertion ignores the sentinel entry itself (``sys.modules["scipy"]
    is None``) and only fails on a genuine (non-``None``) scipy module
    object, which would indicate climt (or something it imports) actually
    performed a real scipy import.
    """
    code = (
        "import sys\n"
        "sys.modules['scipy'] = None\n"
        "import climt\n"
        "from climt import get_default_state, get_grid\n"
        "from climt._components.cork import CorkLongwaveRadiation\n"
        "lw = CorkLongwaveRadiation(optics='correlated_k', table='single_band_unit_lw')\n"
        "state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=20))\n"
        "lw(state)\n"
        "bad = [m for m in sys.modules if (m == 'scipy' or m.startswith('scipy.')) "
        "and sys.modules[m] is not None]\n"
        "assert not bad, bad\n"
        "print('no-scipy OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, (
        f"stdout={result.stdout!r}\nstderr={result.stderr!r}"
    )
    assert "no-scipy OK" in result.stdout

"""Parity of the numba DCMIP initial conditions with the Fortran reference.

The reference outputs in ``tests/fixtures/fortran_reference/dcmip.npz`` were
produced by the Fortran-backed ``DcmipInitialConditions.array_call`` (see
``scripts/generate_fortran_reference.py``).
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

from climt import DcmipInitialConditions

FIXTURES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "fixtures", "fortran_reference")

with open(os.path.join(FIXTURES, "cases.json")) as _f:
    CASES = json.load(_f)["dcmip"]
DCMIP_CASES = sorted(CASES)


def _fixture(case, kind):
    prefix = "{}/{}/".format(case, kind)
    with np.load(os.path.join(FIXTURES, "dcmip.npz")) as data:
        return {key[len(prefix):]: data[key] for key in data.files
                if key.startswith(prefix)}


@pytest.mark.parametrize("case", DCMIP_CASES)
def test_component_matches_fortran(case):
    component = DcmipInitialConditions(**CASES[case]["kwargs"])
    out = component.array_call(_fixture(case, "in"))
    expected_outputs = _fixture(case, "out")
    assert expected_outputs
    for name, expected in expected_outputs.items():
        np.testing.assert_allclose(out[name], expected, rtol=1e-12,
                                   atol=1e-12, err_msg=name)


def test_runs_without_climt_fortran():
    code = (
        "import sys\n"
        "sys.modules['climt_fortran'] = None\n"
        "import climt\n"
        "component = climt.DcmipInitialConditions()\n"
        "state = climt.get_default_state([component])\n"
        "component(state)\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().endswith("ok")


_NO_JIT_RUN = """
import numpy as np
import sympl
import climt

sympl.set_backend(climt.UnytBackend())
component = climt.DcmipInitialConditions()
state = climt.get_default_state(
    [component], grid_state=climt.get_grid(nx=4, ny=2, nz=10))
diagnostics = component(state)
print(float(np.asarray(diagnostics['air_temperature']).ravel()[0]))
"""


def test_component_runs_with_the_jit_disabled():
    """Runs the Pyodide-shaped path: no numba JIT, UnytBackend, one call.

    numba reads NUMBA_DISABLE_JIT at import time, so this has to happen in a
    fresh interpreter.
    """
    result = subprocess.run(
        [sys.executable, "-c", _NO_JIT_RUN],
        env={**os.environ, "NUMBA_DISABLE_JIT": "1"},
        capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, (
        "DcmipInitialConditions failed with the JIT disabled -- it will fail "
        "the same way in Pyodide, which has no numba at all:\n"
        + result.stderr)
    assert np.isfinite(float(result.stdout.strip().splitlines()[-1]))

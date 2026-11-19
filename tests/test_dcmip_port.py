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
from climt._components.dcmip import kernels  # noqa: F401,E501

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
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().endswith("ok")

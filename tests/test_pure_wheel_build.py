import os
import subprocess
import sys
import glob
import zipfile
from pathlib import Path


def test_pure_wheel_builds_and_contains_subpackages_and_data(tmp_path):
    repo = Path(__file__).resolve().parent.parent
    env = dict(os.environ, CLIMT_PURE_PYTHON="1")
    out = tmp_path / "wheelhouse"
    subprocess.run(
        [sys.executable, "-m", "pip", "wheel", ".", "--no-deps", "-w", str(out)],
        cwd=repo, env=env, check=True,
    )
    wheels = glob.glob(str(out / "climt-*-py3-none-any.whl"))
    assert wheels, "expected a py3-none-any wheel"
    names = zipfile.ZipFile(wheels[0]).namelist()
    # subpackages present
    assert any(n.endswith("climt/_components/cork/lw/component.py") for n in names)
    assert any(n.endswith("climt/_core/initialization.py") for n in names)
    # data files present
    assert any("climt/_data/cork/correlated_k/" in n and n.endswith(".nc") for n in names)
    assert any(n.endswith("climt/_data/ozone_profile.npy") for n in names)

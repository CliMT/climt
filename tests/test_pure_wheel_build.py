import os
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FORBIDDEN_SUFFIXES = (".so", ".pyd", ".pyx", ".c", ".a", ".f90", ".F90", ".f", ".mod")
STALE_PROBE = REPO / "build" / "lib" / "climt" / "_stale_probe.so"


def _build_wheel(out, extra_env=None):
    """Build the climt wheel the way the release does (sdist, then wheel from
    the sdist), so stale products in a developer checkout's ``build/`` cannot
    leak into it. A stale probe is planted in ``build/lib/climt`` to prove it."""
    env = dict(os.environ)
    env.pop("CLIMT_PURE_PYTHON", None)
    env.update(extra_env or {})
    created_dirs = [d for d in (STALE_PROBE.parents[2], STALE_PROBE.parents[1], STALE_PROBE.parent)
                    if not d.exists()]
    STALE_PROBE.parent.mkdir(parents=True, exist_ok=True)
    STALE_PROBE.touch()
    try:
        subprocess.run(
            [sys.executable, "-m", "build", "--outdir", str(out), str(REPO)],
            cwd=REPO, env=env, check=True,
        )
    finally:
        STALE_PROBE.unlink()
        for d in reversed(created_dirs):
            d.rmdir()
    wheels = sorted(out.glob("climt-*.whl"))
    assert len(wheels) == 1, wheels
    return wheels[0]


@pytest.fixture(scope="module")
def built_wheel(tmp_path_factory):
    return _build_wheel(tmp_path_factory.mktemp("wheelhouse"))


def test_pure_wheel_builds_and_contains_subpackages_and_data(tmp_path):
    wheel = _build_wheel(tmp_path / "wheelhouse", {"CLIMT_PURE_PYTHON": "1"})
    assert wheel.name.endswith("-py3-none-any.whl"), "expected a py3-none-any wheel"
    names = zipfile.ZipFile(wheel).namelist()
    # subpackages present
    assert any(n.endswith("climt/_components/cork/lw/component.py") for n in names)
    assert any(n.endswith("climt/_core/initialization.py") for n in names)
    # data files present
    assert any("climt/_data/cork/correlated_k/" in n and n.endswith(".nc") for n in names)
    assert any(n.endswith("climt/_data/ozone_profile.npy") for n in names)
    assert not [n for n in names if n.endswith(FORBIDDEN_SUFFIXES)]


def test_wheel_is_pure_and_clean(built_wheel):
    assert built_wheel.name.endswith("-py3-none-any.whl")
    names = zipfile.ZipFile(built_wheel).namelist()
    assert not [n for n in names if n.endswith(FORBIDDEN_SUFFIXES)]
    assert not [n for n in names if n.startswith(("fortran/", "climt_fortran/")) or "/_lib/" in n]


def test_wheel_metadata(built_wheel):
    archive = zipfile.ZipFile(built_wheel)
    metadata = archive.read(next(n for n in archive.namelist() if n.endswith(".dist-info/METADATA"))).decode()
    requires = [line.split(":", 1)[1].strip() for line in metadata.splitlines() if line.startswith("Requires-Dist:")]
    fortran = [r for r in requires if r.startswith("climt-fortran")]
    assert len(fortran) == 1
    assert ">=1.0" in fortran[0] and "<2" in fortran[0]
    assert "sys_platform" in fortran[0] and "emscripten" in fortran[0] and "extra" not in fortran[0]
    assert not [r for r in requires if r.lower().startswith("cython")]
    assert [r for r in requires if r.startswith("numba") and 'extra == "numba"' in r.replace("'", '"')]
    assert "Requires-Python: >=3.10" in metadata
    assert "Provides-Extra: numba" in metadata
    assert "Programming Language :: Python :: 2" not in metadata


def test_wheel_keeps_every_data_file(built_wheel):
    names = zipfile.ZipFile(built_wheel).namelist()
    for needle in (
        "climt/_data/ozone_profile.npy",
        "climt/_data/atmospheric_properties/earth.toml",
        "climt/_data/cork/correlated_k/earth_low_res_lw.nc",
        "climt/_data/cork/correlated_k/earth_low_res_lw.npz",
        "climt/_data/land_mask/earth_landmask_2deg.nc",
        "climt/_data/data_ocean/earth_sst_climatology_1deg.nc",
        "climt/_data/topography/earth_topography_2deg.nc",
        "climt/_data/cork/parmentier/freedman2014.npz",
        "climt/_data/cork/parmentier/solar_composition.npz",
        "climt/_data/cork/stellar_spectra/sun.npz",
        "climt/_data/cork/stellar_spectra/trappist1.npz",
        "climt/_core/fortran.py",
    ):
        assert any(n.endswith(needle) for n in names), needle
    # offline build inputs and dev notes stay out of the wheel
    assert not [n for n in names if "climt/_data/" in n and n.endswith((".cia", ".md"))]

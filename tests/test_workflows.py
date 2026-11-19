from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
ROOT = Path(__file__).resolve().parent.parent
WORKFLOWS = ROOT / ".github" / "workflows"


def load(name):
    return yaml.safe_load((WORKFLOWS / name).read_text())


def text(name):
    return (WORKFLOWS / name).read_text()


def test_ci_runs_both_modes():
    matrix = load("python-package.yml")["jobs"]["build"]["strategy"]["matrix"]
    assert matrix["mode"] == ["full", "pure"]
    assert matrix["python-version"] == ["3.10", "3.11", "3.12"]


def test_core_release_never_compiles():
    body = text("release_climt.yml")
    assert "cibuildwheel" not in body and "qemu" not in body.lower()
    assert "startsWith(github.event.release.tag_name, 'v')" in body
    assert "python -m build" in body


def test_fortran_release_publishes_only_complete_sets():
    workflow = load("release_climt_fortran.yml")
    assert set(workflow["jobs"]["publish"]["needs"]) == {"build_wheels", "build_stub_and_sdist"}
    body = text("release_climt_fortran.yml")
    assert "startsWith(github.event.release.tag_name, 'fortran-v')" in body
    assert "CLIMT_FORTRAN_STUB" in body and "CLIMT_REQUIRE_FORTRAN=1" in body
    assert "macos-15" in body and "cp310-* cp311-* cp312-*" in body
    assert "https://pypi.org/p/climt-fortran" in body


def test_tox_has_both_kinds_of_environment():
    tox = (ROOT / "tox.ini").read_text()
    for name in ("py310", "py311", "py312", "py310-pure", "py311-pure", "py312-pure"):
        assert name in tox
    assert "py39" not in tox
    assert "setup.py" not in tox

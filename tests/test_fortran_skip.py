import os, subprocess, sys
from pathlib import Path

CASE = str(Path(__file__).with_name("fortran_skip_case.py"))


def run(require):
    env = dict(os.environ)
    env.pop("CLIMT_REQUIRE_FORTRAN", None)
    if require:
        env["CLIMT_REQUIRE_FORTRAN"] = "1"
    return subprocess.run([sys.executable, "-m", "pytest", CASE, "-q", "-p", "no:cacheprovider"],
                          capture_output=True, text=True, env=env)


def test_unavailable_is_skipped_by_default():
    result = run(require=False)
    assert "2 skipped" in result.stdout and "1 failed" in result.stdout, result.stdout


def test_unavailable_is_an_error_when_required():
    result = run(require=True)
    assert "skipped" not in result.stdout, result.stdout
    assert "FortranUnavailableError" in result.stdout
    assert result.returncode != 0

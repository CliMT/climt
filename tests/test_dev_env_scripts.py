import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BOOTSTRAP = REPO / "scripts" / "bootstrap_dev_env.sh"
VERIFY = REPO / "scripts" / "verify_dev_env.py"


def _fake_conda(tmp_path, existing_envs=""):
    """A stand-in `conda` that lists `existing_envs` and logs every call."""
    log = tmp_path / "conda_calls.log"
    fake = tmp_path / "conda"
    fake.write_text(
        "#!/bin/bash\n"
        'echo "$@" >> "{log}"\n'
        'if [ "$1" = "env" ] && [ "$2" = "list" ]; then\n'
        '  printf "# conda environments:\\nbase  /x/base\\n{envs}"\n'
        "fi\n".format(log=log, envs=existing_envs))
    fake.chmod(0o755)
    return fake, log


def test_requirements_file_pins_the_jax_stack_and_omits_jax_metal():
    text = (REPO / "envs" / "climt-dev.in").read_text()
    for line in ("jax==0.10.1", "jaxlib==0.10.1", "s2fft==1.4.0",
                 "flax==0.12.7", "optax==0.2.8", "numba>=0.64", "numpy>=2.0",
                 "pyzmq", "zarr", "panel", "holoviews"):
        assert line in text.splitlines()
    assert "jax-metal" not in text


def test_bootstrap_refuses_to_touch_an_existing_environment(tmp_path):
    fake, log = _fake_conda(tmp_path, existing_envs="climt-dev  /x/envs/climt-dev\\n")
    env = dict(os.environ, CONDA_EXE=str(fake))
    result = subprocess.run(["bash", str(BOOTSTRAP)], env=env,
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert "already exists; pass --recreate" in result.stderr
    calls = log.read_text()
    assert "remove" not in calls and "create" not in calls


def test_bootstrap_dry_run_runs_nothing_and_names_each_step(tmp_path):
    fake, log = _fake_conda(tmp_path)
    env = dict(os.environ, CONDA_EXE=str(fake), CLIMT_DIR=str(REPO))
    env.pop("GFS_DYCORE_DIR", None)
    result = subprocess.run(["bash", str(BOOTSTRAP), "--dry-run"], env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0
    out = result.stdout
    assert "DRY-RUN: " in out
    assert "python=3.11" in out
    assert "pip install -e" in out
    assert "GFS_DYCORE_DIR not set or not a directory" in out
    calls = log.read_text()
    assert "create" not in calls and "run" not in calls


def test_bootstrap_dry_run_installs_dycore_without_deps(tmp_path):
    fake, _ = _fake_conda(tmp_path)
    dycore = tmp_path / "gfs"
    dycore.mkdir()
    env = dict(os.environ, CONDA_EXE=str(fake), CLIMT_DIR=str(REPO),
               GFS_DYCORE_DIR=str(dycore))
    result = subprocess.run(["bash", str(BOOTSTRAP), "--dry-run"], env=env,
                            capture_output=True, text=True)
    assert "--no-deps -e {}".format(dycore) in result.stdout


def _load_verify():
    spec = importlib.util.spec_from_file_location("verify_dev_env", VERIFY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_verify_reports_every_check_by_name_and_never_raises():
    checks = _load_verify().run_checks()
    assert [c["name"] for c in checks] == [
        "climt import", "numba", "jax cpu x64", "no jax-metal",
        "climt component", "io stack", "fortran extensions",
        "gfs_dynamical_core"]
    assert [c["required"] for c in checks] == [True] * 6 + [False] * 2
    for c in checks:
        assert isinstance(c["ok"], bool) and isinstance(c["detail"], str)


def test_verify_json_output_and_exit_status_follow_required_checks():
    result = subprocess.run([sys.executable, str(VERIFY), "--json"],
                            capture_output=True, text=True)
    checks = json.loads(result.stdout)
    required_ok = all(c["ok"] for c in checks if c["required"])
    assert result.returncode == (0 if required_ok else 1)


def test_verify_marks_a_failing_optional_check_as_missing_optional(monkeypatch, capsys):
    module = _load_verify()
    monkeypatch.setattr(module, "run_checks", lambda: [
        {"name": "climt import", "required": True, "ok": True, "detail": "x"},
        {"name": "fortran extensions", "required": False, "ok": False, "detail": "no"}])
    assert module.main([]) == 0
    assert "MISSING (optional)" in capsys.readouterr().out
    assert module.main(["--strict"]) == 1

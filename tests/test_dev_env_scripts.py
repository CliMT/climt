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


def test_bootstrap_dry_run_links_dycore_with_a_pth_file(tmp_path):
    fake, _ = _fake_conda(tmp_path)
    dycore = tmp_path / "gfs"
    dycore.mkdir()
    env = dict(os.environ, CONDA_EXE=str(fake), CLIMT_DIR=str(REPO),
               GFS_DYCORE_DIR=str(dycore))
    result = subprocess.run(["bash", str(BOOTSTRAP), "--dry-run"], env=env,
                            capture_output=True, text=True)
    line = next(l for l in result.stdout.splitlines()
                if "gfs_dynamical_core.pth" in l)
    assert str(dycore.resolve()) in line
    assert "pip" not in line


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


def test_bootstrap_continues_to_verify_when_the_dycore_link_fails(tmp_path):
    # Run a copy of the script from a scratch repository so the lock file it
    # writes lands in tmp_path, not in this checkout.
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    (repo / "envs").mkdir()
    script = repo / "scripts" / "bootstrap_dev_env.sh"
    script.write_text(BOOTSTRAP.read_text())
    (repo / "envs" / "climt-dev.in").write_text("numpy\n")
    log = tmp_path / "conda_calls.log"
    fake = tmp_path / "conda"
    fake.write_text(
        "#!/bin/bash\n"
        'echo "$@" >> "{log}"\n'
        'case "$*" in *sysconfig*) exit 1 ;; esac\n'
        'if [ "$1" = "env" ] && [ "$2" = "list" ]; then\n'
        '  printf "# conda environments:\\nbase  /x/base\\n"\n'
        "fi\n".format(log=log))
    fake.chmod(0o755)
    dycore = tmp_path / "gfs"
    dycore.mkdir()
    env = dict(os.environ, CONDA_EXE=str(fake), CLIMT_DIR=str(REPO),
               GFS_DYCORE_DIR=str(dycore))
    result = subprocess.run(["bash", str(script)], env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert "gfs_dynamical_core link failed; continuing" in result.stderr
    calls = log.read_text().splitlines()
    dycore_call = next(i for i, c in enumerate(calls) if "sysconfig" in c)
    later = calls[dycore_call + 1:]
    assert any("pip freeze" in c for c in later)
    assert any("verify_dev_env.py" in c for c in later)


def test_bootstrap_links_dycore_when_conda_run_adds_a_trailing_blank_line(tmp_path):
    # Real `conda run` prints a blank line after the command's output; the
    # site-packages path must still be found.
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    (repo / "envs").mkdir()
    script = repo / "scripts" / "bootstrap_dev_env.sh"
    script.write_text(BOOTSTRAP.read_text())
    (repo / "envs" / "climt-dev.in").write_text("numpy\n")
    purelib = tmp_path / "site-packages"
    purelib.mkdir()
    log = tmp_path / "conda_calls.log"
    fake = tmp_path / "conda"
    fake.write_text(
        "#!/bin/bash\n"
        'echo "$@" >> "{log}"\n'
        'case "$*" in *sysconfig*) printf "%s\\n\\n" "{purelib}"; exit 0 ;; esac\n'
        'if [ "$1" = "env" ] && [ "$2" = "list" ]; then\n'
        '  printf "# conda environments:\\nbase  /x/base\\n"\n'
        "fi\n".format(log=log, purelib=purelib))
    fake.chmod(0o755)
    dycore = tmp_path / "gfs"
    dycore.mkdir()
    env = dict(os.environ, CONDA_EXE=str(fake), CLIMT_DIR=str(REPO),
               GFS_DYCORE_DIR=str(dycore))
    result = subprocess.run(["bash", str(script)], env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert "link failed" not in result.stderr
    pth = purelib / "gfs_dynamical_core.pth"
    assert pth.read_text() == str(dycore.resolve()) + "\n"


def test_verify_reports_missing_show_versions_without_raising(monkeypatch, capsys):
    import types
    module = _load_verify()
    monkeypatch.setattr(module, "run_checks", lambda: [])
    fake = types.ModuleType("climt")
    fake.__version__ = "9.9.9"
    fake.__file__ = "/somewhere/climt/__init__.py"
    monkeypatch.setitem(sys.modules, "climt", fake)
    assert module.main([]) == 0
    out = capsys.readouterr().out
    assert ("climt.show_versions unavailable in climt 9.9.9 at "
            "/somewhere/climt/__init__.py") in out
    assert "failed" not in out


def _scratch_bootstrap(tmp_path, fake_body):
    """Copy the bootstrap into a scratch repo and run it against a fake conda.

    `fake_body` is bash run after the call is logged; it sees the call in "$*".
    Returns (result, calls, lock_path).
    """
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    (repo / "envs").mkdir()
    script = repo / "scripts" / "bootstrap_dev_env.sh"
    script.write_text(BOOTSTRAP.read_text())
    (repo / "envs" / "climt-dev.in").write_text("numpy\n")
    log = tmp_path / "conda_calls.log"
    fake = tmp_path / "conda"
    fake.write_text(
        "#!/bin/bash\n"
        'echo "PURE=${{CLIMT_PURE_PYTHON:-}} $*" >> "{log}"\n'
        'if [ "$1" = "env" ] && [ "$2" = "list" ]; then\n'
        '  printf "# conda environments:\\nbase  /x/base\\n"\n'
        "fi\n".format(log=log) + fake_body)
    fake.chmod(0o755)
    dycore = tmp_path / "gfs"
    dycore.mkdir()
    env = dict(os.environ, CONDA_EXE=str(fake), CLIMT_DIR=str(REPO),
               GFS_DYCORE_DIR=str(dycore))
    result = subprocess.run(["bash", str(script)], env=env,
                            capture_output=True, text=True)
    return result, log.read_text().splitlines(), repo / "envs" / "climt-dev.lock.txt"


def test_lock_never_pins_the_path_linked_dycore(tmp_path):
    # The .pth link makes pip see the dycore checkout's egg-info as an installed,
    # non-editable package; the lock must not pin it.
    result, _, lock = _scratch_bootstrap(
        tmp_path,
        'case "$*" in *"pip freeze"*) '
        'printf "numpy==2.0.0\\ngfs_dynamical_core==0.3.0\\ngfs-dynamical-core==0.3.0\\n" ;; esac\n')
    assert result.returncode == 0, result.stderr
    lines = lock.read_text().splitlines()
    assert "numpy==2.0.0" in lines
    assert not [l for l in lines if l.lower().replace("-", "_").startswith("gfs_dynamical_core")]


def test_bootstrap_falls_back_to_pure_python_when_the_compiled_build_fails(tmp_path):
    result, calls, _ = _scratch_bootstrap(
        tmp_path,
        'case "$*" in *"pip install -e"*) '
        '[ "${CLIMT_PURE_PYTHON:-}" = "1" ] || exit 1 ;; esac\n')
    assert result.returncode == 0, result.stderr
    assert "compiled build failed; falling back to pure Python" in result.stdout
    installs = [c for c in calls if "pip install -e" in c]
    assert len(installs) == 2
    assert installs[0].startswith("PURE= ") and installs[1].startswith("PURE=1 ")

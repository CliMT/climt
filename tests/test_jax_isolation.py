import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SHIM = str(Path(__file__).resolve().parent / "no_jax_shim")


def run(code, block_jax=False):
    env = dict(os.environ)
    paths = [str(REPO)] + ([SHIM] if block_jax else [])
    env["PYTHONPATH"] = os.pathsep.join(paths[::-1] + [env.get("PYTHONPATH", "")])
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)


def test_import_climt_does_not_import_jax():
    result = run("import sys, climt; assert 'jax' not in sys.modules, 'jax was imported'")
    assert result.returncode == 0, result.stderr


def test_climt_works_when_jax_cannot_be_imported():
    code = (
        "import climt, numpy as np\n"
        "from climt._core.jax_support import is_jax_array\n"
        "assert is_jax_array(np.zeros(3)) is False\n"
        "climt.EmanuelConvectionPython(); climt.CorkLongwaveRadiation()\n"
    )
    result = run(code, block_jax=True)
    assert result.returncode == 0, result.stderr


def test_no_module_level_jax_import_in_loaded_modules():
    code = (
        "import sys, climt, re, pathlib\n"
        "bad = []\n"
        "for name, module in list(sys.modules.items()):\n"
        "    path = getattr(module, '__file__', None) or ''\n"
        "    if name.startswith('climt') and path.endswith('.py'):\n"
        "        if re.search(r'^(import jax|from jax)', pathlib.Path(path).read_text(), re.M):\n"
        "            bad.append(name)\n"
        "assert not bad, bad\n"
    )
    result = run(code)
    assert result.returncode == 0, result.stderr


def test_climt_never_changes_jax_config():
    import ast

    offenders = []
    for p in (REPO / "climt").rglob("*.py"):
        tree = ast.parse(p.read_text(), filename=str(p))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "update"
                and isinstance(node.func.value, ast.Attribute)
                and node.func.value.attr == "config"
            ):
                offenders.append(str(p))
                break
    assert offenders == []

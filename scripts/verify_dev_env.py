"""Check that the current Python environment is a working climt-dev environment.

Usage: python scripts/verify_dev_env.py [--strict] [--json]

Exits 0 when every required check passes (with --strict, every optional
check too), 1 otherwise.
"""
import argparse
import json
import sys


def _climt_import():
    import climt
    return climt.__version__


def _numba():
    import numba
    import numpy as np

    @numba.njit
    def total(values):
        result = 0.0
        for value in values:
            result += value
        return result

    value = total(np.arange(10.0))
    if value != 45.0:
        raise AssertionError("expected 45.0, got {!r}".format(value))
    return numba.__version__


def _jax_cpu_x64():
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    platform = jax.devices()[0].platform
    if platform != "cpu":
        raise AssertionError("first device platform is {!r}, expected 'cpu'".format(platform))
    value = jnp.arange(3.0).sum()
    if float(value) != 3.0 or str(value.dtype) != "float64":
        raise AssertionError("arange(3.0).sum() gave {!r} ({})".format(value, value.dtype))
    return "jax {}".format(jax.__version__)


def _no_jax_metal():
    from importlib.metadata import PackageNotFoundError, version
    try:
        installed = version("jax-metal")
    except PackageNotFoundError:
        return "jax-metal not installed"
    raise AssertionError("jax-metal {} is installed".format(installed))


def _climt_component():
    import numpy as np
    import climt
    component = climt.GrayLongwaveRadiation()
    state = climt.get_default_state([component])
    tendencies, _ = component(state)
    values = np.asarray(tendencies["air_temperature"])
    if not np.all(np.isfinite(values)):
        raise AssertionError("air_temperature tendency has non-finite values")
    return "GrayLongwaveRadiation air_temperature tendency finite"


def _io_stack():
    import holoviews
    import panel
    import zarr
    import zmq
    return "zmq {}, zarr {}, panel {}, holoviews {}".format(
        zmq.__version__, zarr.__version__, panel.__version__, holoviews.__version__)


def _fortran_extensions():
    import climt
    if not climt.has_fortran_extensions():
        raise AssertionError("climt.has_fortran_extensions() is False")
    return "available"


def _gfs_dynamical_core():
    import gfs_dynamical_core
    import gfs_dynamical_core.jax.stepper  # noqa: F401
    return getattr(gfs_dynamical_core, "__version__", "imported")


CHECKS = [
    ("climt import", True, _climt_import),
    ("numba", True, _numba),
    ("jax cpu x64", True, _jax_cpu_x64),
    ("no jax-metal", True, _no_jax_metal),
    ("climt component", True, _climt_component),
    ("io stack", True, _io_stack),
    ("fortran extensions", False, _fortran_extensions),
    ("gfs_dynamical_core", False, _gfs_dynamical_core),
]


def run_checks():
    """Run every check; return a list of dicts with name, required, ok, detail."""
    results = []
    for name, required, check in CHECKS:
        try:
            detail = str(check())
            ok = True
        except Exception as exc:
            detail = repr(exc)
            ok = False
        results.append({"name": name, "required": required, "ok": ok, "detail": detail})
    return results


def _status(check):
    if check["ok"]:
        return "OK"
    return "FAIL" if check["required"] else "MISSING (optional)"


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--strict", action="store_true",
                        help="also fail when an optional check fails")
    parser.add_argument("--json", action="store_true",
                        help="print the check results as JSON and nothing else")
    args = parser.parse_args(argv)

    checks = run_checks()
    if args.json:
        print(json.dumps(checks))
    else:
        for check in checks:
            print("{:<22} {}  {}".format(check["name"], _status(check), check["detail"]))
        try:
            import climt
            if hasattr(climt, "show_versions"):
                climt.show_versions()
            else:
                print("climt.show_versions unavailable in climt {} at {}".format(
                    getattr(climt, "__version__", "unknown"), climt.__file__))
        except Exception as exc:
            print("climt.show_versions() failed: {!r}".format(exc))

    considered = checks if args.strict else [c for c in checks if c["required"]]
    return 0 if all(c["ok"] for c in considered) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

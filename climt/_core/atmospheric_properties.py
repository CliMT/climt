"""Multi-planet atmospheric properties via TOML profiles.

Loads named or custom atmospheric profiles that set sympl constants
for gravitational acceleration, molar masses, heat capacities, etc.
Supports snapshot/restore so users can switch between planets in a
single session.
"""
import os

from sympl import get_constant, set_constant

from .exceptions import ConstantNotFoundError

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # fallback for 3.10

# Path to built-in profile directory
_BUILTIN_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    "_data",
    "atmospheric_properties",
)

# Snapshot stack (list of dicts). Each entry is a dict of
# {constant_name: (value, units)} captured before a profile was loaded.
_snapshot_stack = []

_active_profile = {"name": None, "path": None}


def _resolve_profile_path(name_or_path):
    """Return the TOML file path for a profile name or file path."""
    if os.path.isfile(name_or_path):
        return name_or_path
    builtin = os.path.join(_BUILTIN_DIR, f"{name_or_path}.toml")
    if os.path.isfile(builtin):
        return builtin
    raise FileNotFoundError(
        f"Atmospheric profile '{name_or_path}' not found. "
        f"Looked for file '{name_or_path}' and built-in '{builtin}'."
    )


def _parse_toml(path):
    """Parse a TOML profile and return a flat dict of {name: (value, units)}."""
    with open(path, "rb") as f:
        data = tomllib.load(f)
    constants = {}
    for section in data.values():
        if isinstance(section, dict):
            for key, val in section.items():
                if isinstance(val, dict) and "value" in val and "units" in val:
                    constants[key] = (val["value"], val["units"])
    return constants


def _snapshot_constants(keys):
    """Capture current values of the given constants from sympl."""
    snap = {}
    for key, (_, units) in keys.items():
        try:
            snap[key] = (get_constant(key, units), units)
        except (KeyError, ValueError):
            # Constant didn't exist before — record as absent
            snap[key] = None
    return snap


def load_atmospheric_properties(name_or_path):
    """Load an atmospheric profile and set sympl constants accordingly.

    Takes a snapshot of the current constants before applying changes,
    so that reset_atmospheric_properties() can restore them.

    Args:
        name_or_path: Built-in profile name (e.g., "earth", "mars",
            "hot_jupiter") or path to a custom .toml file.
    """
    path = _resolve_profile_path(name_or_path)
    constants = _parse_toml(path)

    # Snapshot current state of the constants we're about to overwrite
    _snapshot_stack.append(_snapshot_constants(constants))

    _active_profile["name"] = os.path.splitext(os.path.basename(path))[0]
    _active_profile["path"] = path

    # Apply new constants
    for key, (value, units) in constants.items():
        set_constant(key, value, units)


def reset_atmospheric_properties():
    """Restore constants to the state before the last load_atmospheric_properties call.

    Raises:
        RuntimeError: If no profile has been loaded (nothing to reset).
    """
    if not _snapshot_stack:
        raise RuntimeError(
            "No atmospheric profile snapshot to restore. "
            "Call load_atmospheric_properties() first."
        )
    snap = _snapshot_stack.pop()
    for key, val in snap.items():
        if val is not None:
            value, units = val
            set_constant(key, value, units)
    _active_profile["name"] = None
    _active_profile["path"] = None


def get_constant_checked(name, units):
    """Wrapper around sympl's get_constant that raises ConstantNotFoundError
    with a helpful message when the constant is missing from the active profile."""
    try:
        return get_constant(name, units)
    except (KeyError, ValueError) as exc:
        raise ConstantNotFoundError(
            name,
            profile_name=_active_profile["name"],
            profile_path=_active_profile["path"],
        ) from exc

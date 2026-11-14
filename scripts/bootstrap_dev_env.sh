#!/bin/bash
# Create the climt-dev conda environment.
#
# Usage: scripts/bootstrap_dev_env.sh [--recreate] [--relock] [--dry-run]
#
# Environment variables:
#   CLIMT_DIR           climt checkout to install (default: this repository)
#   GFS_DYCORE_DIR      gfs_dynamical_core checkout (optional; installed --no-deps)
#   CLIMT_DEV_ENV_NAME  environment name (default: climt-dev)
#   CONDA_EXE           conda executable (default: conda on PATH)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REQUIREMENTS="$REPO_DIR/envs/climt-dev.in"
LOCK="$REPO_DIR/envs/climt-dev.lock.txt"

CLIMT_DIR="${CLIMT_DIR:-$REPO_DIR}"
ENV_NAME="${CLIMT_DEV_ENV_NAME:-climt-dev}"
CONDA="${CONDA_EXE:-conda}"

RECREATE=0
RELOCK=0
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --recreate) RECREATE=1 ;;
        --relock) RELOCK=1 ;;
        --dry-run) DRY_RUN=1 ;;
        *) echo "unknown argument: $arg" >&2
           echo "usage: $0 [--recreate] [--relock] [--dry-run]" >&2
           exit 64 ;;
    esac
done

run() {
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "DRY-RUN: $*"
    else
        "$@"
    fi
}

pip_in_env() {
    run "$CONDA" run -n "$ENV_NAME" python -m pip "$@"
}

# 2. Never touch an existing environment without --recreate.
ENV_LIST="$("$CONDA" env list)"
if printf '%s\n' "$ENV_LIST" | awk '{print $1}' | grep -Fxq -- "$ENV_NAME"; then
    if [ "$RECREATE" -ne 1 ]; then
        echo "environment '$ENV_NAME' already exists; pass --recreate to rebuild it" >&2
        exit 2
    fi
    run "$CONDA" env remove -n "$ENV_NAME" -y
fi

# 3. Create the environment.
run "$CONDA" create -n "$ENV_NAME" -y python=3.11 pip

# 4. Install requirements.
WRITE_LOCK=0
if [ "$RELOCK" -eq 1 ] || [ ! -f "$LOCK" ]; then
    WRITE_LOCK=1
    pip_in_env install -r "$REQUIREMENTS"
else
    pip_in_env install -r "$LOCK"
fi

# 5. Install climt, falling back to the pure-Python build.
if ! pip_in_env install -e "$CLIMT_DIR"; then
    echo "compiled build failed; falling back to pure Python"
    run env CLIMT_PURE_PYTHON=1 "$CONDA" run -n "$ENV_NAME" python -m pip install -e "$CLIMT_DIR"
fi

# 6. Install gfs_dynamical_core without its (stale) pinned dependencies.
if [ -n "${GFS_DYCORE_DIR:-}" ] && [ -d "$GFS_DYCORE_DIR" ]; then
    pip_in_env install --no-deps -e "$GFS_DYCORE_DIR"
else
    echo "GFS_DYCORE_DIR not set or not a directory; skipping gfs_dynamical_core"
fi

# 7. Write the lock file.
if [ "$WRITE_LOCK" -eq 1 ]; then
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "DRY-RUN: $CONDA run -n $ENV_NAME python -m pip freeze --exclude-editable > $LOCK"
    else
        "$CONDA" run -n "$ENV_NAME" python -m pip freeze --exclude-editable > "$LOCK.tmp"
        mv "$LOCK.tmp" "$LOCK"
    fi
fi

# 8. Verify the new environment.
status=0
run "$CONDA" run -n "$ENV_NAME" python "$REPO_DIR/scripts/verify_dev_env.py" || status=$?
exit "$status"

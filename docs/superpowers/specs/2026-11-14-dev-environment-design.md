# One development environment for climt, numba and JAX — design

Date: 2026-11-14. Size: small change. Branch: `feature/dev-environment`.

## Why

Two conda environments exist on the maintainer's machine and neither runs
everything. `climt` (Python 3.10) has the editable climt checkout and numba,
but JAX 0.5.3 with `jax-metal` (JAX defaults to a broken Metal device) and a
scipy whose Fortran extensions do not load. `gfs_dycore_develop` (Python
3.11) has a healthy JAX 0.10.1 and scipy, but no numba, no pytest, and climt
0.18.4 from PyPI. Every later project (async I/O, the CCPP generator, the
package split, the JAX merge) needs current climt, numba and the JAX
dynamical core in one process. Users reporting bugs (e.g. issue #174, a
NumPy binary incompatibility) also have no one-line way to report versions.

## Goals

- One conda environment, `climt-dev`, in which current climt (with compiled
  extensions when a Fortran toolchain is present), numba, JAX on CPU, the
  `gfs_dynamical_core` checkout, pytest, and the optional I/O stack
  (pyzmq, zarr, panel, holoviews) all import and run in the same process.
- Reproducible from files committed to the repo: a short list of top-level
  requirements, a lock file with exact versions, one bootstrap script.
- `climt.show_versions()` for bug reports.
- A verification script that says, per capability, whether the environment
  has it.

## Non-goals

- Changing CI, tox, `setup.py` or the published package's dependencies.
- Removing or modifying the existing `climt` and `gfs_dycore_develop`
  environments.
- GPU or Metal JAX. `jax-metal` is deliberately absent.
- Cross-platform hash-pinned lock files. The lock pins versions only.

## Anti-goals

- The environment must not silently lose capability: if the Fortran build or
  the dycore install fails, the bootstrap continues but the verification
  script reports it by name.
- `show_versions()` must not import the packages it reports (importing JAX
  or a broken scipy is exactly what fails on the machines we want reports
  from).

## Approaches

1. **conda for Python only, pip for everything, `pip freeze` as the lock**
   (chosen). No new tool to install; conda, mamba and pip are already on the
   machine and conda-lock, uv and pixi are not. The lock is versions-only
   and platform-agnostic enough for macOS arm64 and Linux x86-64.
2. **conda-lock / pixi.** Proper multi-platform solver locks, but adds a
   tool to the base environment and a second resolver to keep in step with
   `setup.py`. Revisit when the package split introduces conda-forge
   feedstocks.
3. **uv workspace.** The likely end state (legoESM uses one), but it belongs
   with the package split, not before it.

## Decisions

- Python 3.11: the newest version for which both repos already have compiled
  artefacts, and the floor for current JAX.
- JAX stack pinned to the combination already proven with the dycore in
  `gfs_dycore_develop`: `jax==0.10.1`, `jaxlib==0.10.1`, `s2fft==1.4.0`,
  `flax==0.12.7`, `optax==0.2.8`. `gfs_dynamical_core`'s own metadata pins
  `jax==0.5.3`, so it is installed with `--no-deps`.
- climt is installed editable from `CLIMT_DIR` (default: the repo that holds
  the script). On the maintainer's machine it is run with
  `CLIMT_DIR=~/github/climt` so the environment does not point into a
  worktree that will later be deleted.
- The dycore checkout location comes from `GFS_DYCORE_DIR`; when unset or
  missing, that step is skipped and reported.

## Left to the builder

Exact wording of the CONTRIBUTING section; the layout of the verification
script's output table; how the bootstrap script detects an existing
environment (it must not delete one without `--recreate`).

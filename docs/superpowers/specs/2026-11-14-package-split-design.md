# Splitting climt along the build toolchain — design

Date: 2026-11-14. Size: project. Branch: `feature/package-split`.

## Intent

climt is one distribution that carries two very different kinds of code:
pure-Python and numba components that change daily, and five Cython
extensions over Fortran libraries that almost never change. In the year to
2026-11-14, 4 of 539 commits on `develop` touched `climt/_lib` or a `.pyx`
file. Yet the compiled code sets the terms for everything:

- `setup.py` is 300 lines of compiler discovery around 20 lines of metadata,
  and installs Cython and NumPy with pip from inside itself.
- Every release must get `cibuildwheel` through on two platforms. The
  0.30.0 release failed on a macOS `delocate` error unrelated to any change
  in that release.
- A development install tries the Fortran build. On the maintainer's machine
  it currently fails (Homebrew gcc 15 cannot find `stdlib.h`) and falls back
  to pure Python, so the Fortran components are missing from the dev
  environment.
- A user with a binary problem (issue #174, a NumPy ABI mismatch) cannot
  keep climt and drop only the binaries.

`develop` already ships two wheels from one source (compiled, and
`py3-none-any` built with `CLIMT_PURE_PYTHON=1`) and already degrades
quietly when the extensions are absent. That is the right behaviour for
users. This project keeps the behaviour and moves the mechanism: the
binaries go into a companion distribution that is released on its own
schedule, and `climt` itself becomes pure Python always.

There is a precedent in the project: the dynamical core left for
`gfs_dynamical_core` for the same reason.

## Goals

- `climt` is a pure-Python distribution: one `py3-none-any` wheel and an
  sdist that builds with no compiler, Cython or NumPy at build time,
  described by a declarative `pyproject.toml`.
- A companion distribution, `climt-fortran` (import name `climt_fortran`),
  in the same repository under `fortran/`. It holds the five extension
  modules and the Fortran sources, has its own version, and has its own
  release workflow.
- `pip install climt` gives exactly what it gives today on every platform:
  the Fortran-backed components where compiled wheels exist, a working
  pure-Python climt everywhere else, and never a failed install for want of
  a compiler.
- A user or developer can tell what they have:
  `climt.has_fortran_extensions()` (kept), and a new
  `climt.list_components()` that says, per component, whether it can be
  used and what it needs.
- Optional stacks are installable by name: `climt[numba]`, and whichever of
  `jax`, `stream`, `dashboard`, `ccpp` exist when this is built.
- CI runs the suite twice: without the binaries and with them built from
  this repository's sources.

## Non-goals

- Changing how the Fortran is built. It stays make + setuptools + Cython;
  the build logic moves, it is not rewritten.
- One distribution per component family (`climt-rrtmg`, ...), or separate
  repositories.
- Moving the Python component classes (`RRTMGLongwave`, ...) out of
  `climt`. Only the compiled modules move.
- A conda-forge recipe. The `enable-conda-packaging` branch stays unmerged;
  a `noarch` recipe for the core becomes easy afterwards.
- Windows or Linux-ARM binaries (there are none today).
- Making numba a default dependency. It stays optional; see "Open".
- Publishing anything. Uploading to PyPI is the maintainer's action.

## Anti-goals

- A platform that gets RRTMG from `pip install climt` today losing it.
- An install failing because no Fortran compiler is present, unless the
  user explicitly asked for a source build.
- The core and the companion disagreeing about the extension interface
  without saying so.
- A test run that silently skips every Fortran test in the job that is
  meant to exercise them.

## Approaches

1. **Pure core plus one compiled companion that ships binaries where it can
   and an empty stub wheel elsewhere** (chosen). `climt` depends on
   `climt-fortran`. For each companion version PyPI holds platform wheels
   and one `py3-none-any` stub wheel with no binaries. pip prefers the
   platform wheel where one fits and otherwise takes the stub, so the
   install always succeeds and never compiles. This is the two-wheel
   mechanism climt uses for itself today, moved to the package that
   actually needs it. Core releases no longer involve a compiler; the
   companion is released only when the Fortran, the Cython or the NumPy ABI
   changes. A development install of the core pulls the companion's binary
   wheel from PyPI, so a developer gets RRTMG without a local toolchain.
2. **Companion as an opt-in extra** (`pip install "climt[fortran]"`).
   Simplest packaging, and probably the right end state once CORK and the
   Python Emanuel are the defaults. Today it would break every tutorial and
   17 of the examples for anyone who types `pip install climt`.
3. **Companion as a hard dependency under environment markers** (Linux
   x86-64, macOS arm64, CPython 3.10–3.12). Markers cannot express the
   macOS version or musl, so an arm64 Mac on macOS 14, or Alpine Linux,
   would stop being able to install climt at all: pip would find no wheel
   and try to compile.
4. **Leave it as it is and tidy `setup.py`.** No migration, but releases
   stay coupled to the compiled build and the developer's toolchain
   problem stays.
5. **One distribution per component family.** Five build pipelines for code
   that changes four times a year.

The 2026-03 Pyodide design chose "a single `climt` package, not a separate
`climt_lite`". That holds: `import climt` remains the only API, every
component class stays in it, and the pure wheel remains the complete
package. What moves is where the binaries are shipped from.

## Design

### Layout

```
pyproject.toml                 core metadata (setup.py is deleted)
climt/                         pure Python only; no .pyx, no _lib
climt/_core/fortran.py         the only place that imports climt_fortran
climt/_core/registry.py        list_components()
fortran/
  pyproject.toml               build requirements
  setup.py                     the build logic that was in the root setup.py
  MANIFEST.in
  README.rst
  climt_fortran/__init__.py    version, interface version, module table
  climt_fortran/_rrtmg_lw.pyx  (and _rrtmg_sw, _emanuel_convection,
                                _simple_physics, _dcmip)
  _lib/                        Fortran sources and Makefiles (from climt/_lib)
  tests/
```

Files are moved with `git mv` so history follows them. The unused
`climt/_components/_berger_solar_insolation.pyx` is deleted
(`BergerSolarInsolation` has been pure Python since the numba work).

### The companion

`climt_fortran/__init__.py` exposes `__version__`, `INTERFACE_VERSION`
(an integer, 1), `MODULES` (the five module names) and `available()`,
which reports per module whether it can be found. It imports none of them.

It builds two ways from the same directory: normally (make, then Cython
extensions, as today) and with `CLIMT_FORTRAN_STUB=1`, which produces a
`py3-none-any` wheel containing only `__init__.py`.

`INTERFACE_VERSION` changes whenever a function that the core calls
changes its name, arguments or return values. The core requires
`climt-fortran>=1.0,<2` and checks `INTERFACE_VERSION == 1` at import.

### The core

`climt/_core/fortran.py` is the single point of contact. It answers three
questions: give me extension module X or `None`; what is the companion's
status (not installed / stub without binaries / wrong interface version /
binaries present but failing to load / available, with its version); and
what should the error say. The five
component modules ask it for their extension instead of importing a
sibling module. A component whose extension is missing raises
`FortranUnavailableError`, a subclass of `ImportError` (so existing
`except ImportError` code keeps working), whose message names the
component, states which reason applies, and gives the one command that
fixes it.

`climt` depends on `climt-fortran>=1.0,<2` everywhere except Pyodide
(`sys_platform != "emscripten"`), where the in-browser demos should not
need a second download. The core still works if the companion is absent
altogether (`pip install --no-deps`, or after `pip uninstall
climt-fortran`), which is the escape hatch for a broken binary.

### Discoverability

`climt.list_components()` returns one record per component: its name, its
kind (the sympl base class), whether it can be instantiated here, and what
it needs if not. It covers climt's own components and any that another
installed package registers under the entry-point group
`climt.components`, which is how `gfs_dynamical_core` and
CCPP-generated packages can appear in the list later. It imports nothing
optional. `repr` of the result is a readable table.

### Tests

Tests that reach a Fortran-backed component are skipped, not failed, when
the binaries are absent: `tests/conftest.py` turns a
`FortranUnavailableError` raised during a test into a skip. With
`CLIMT_REQUIRE_FORTRAN=1` in the environment it is left as an error, and
the with-binaries CI job sets that, so it cannot pass by skipping.

### CI and release

- CI, per Python version and OS: a **pure** job installs the stub companion
  and the core and runs the suite with no compiler involved; a **full** job
  builds the companion from `fortran/` and runs the suite with
  `CLIMT_REQUIRE_FORTRAN=1`.
- `release_climt.yml` (on a release tagged `v*`): sdist and pure wheel for
  `climt`, attached to the GitHub release and published to PyPI. No
  `cibuildwheel`.
- `release_climt_fortran.yml` (on a release tagged `fortran-v*`): platform
  wheels with `cibuildwheel` (same matrix as today: CPython 3.10–3.12,
  Linux x86-64, macOS arm64 on the pinned `macos-15` runner), the stub
  wheel, and an sdist. It publishes only if every build succeeded, so a
  version never exists with a stub but without its binaries.
- First release order: `climt-fortran` 1.0.0, then the next `climt`
  release. Until the companion is on PyPI, a source checkout installs
  `fortran/` first and the core second.

## Constraints

- Python 3.10 or newer; the classifiers stop claiming Python 2.
- `import climt` stays silent whatever the companion's state, including a
  binary that fails to load.
- The pure wheel must contain no `.so`, `.pyx`, `.c` or Fortran file, and
  must still contain every data file it contains today.
- `CLIMT_PURE_PYTHON=1` remains accepted (and now changes nothing), because
  the docs' in-browser boot script and `scripts/serve_wheel.py` set it.
- Whatever `extras_require` exists in `setup.py` when this is built (the
  async I/O and CCPP plans add some) is carried over unchanged.
- Runtime dependencies are otherwise unchanged, except that `cython` is no
  longer one.

## Open

- **numba as a default dependency.** `pip install climt` today installs no
  numba, so the kernels run as plain Python unless the user knows to add
  it. Making it a default dependency would fix that but ties climt's
  installability on a new Python release to numba's. This design adds only
  the extra `climt[numba]` and reports numba's absence in
  `list_components()`'s footer; the default is the maintainer's call.
- **The name `climt-fortran`** must be free on PyPI and needs a trusted
  publisher configured before the first release.

## Left to the builder

The wording of error messages beyond the three required facts (component,
reason, fix); the table
layout of `list_components()`; how the stub build is expressed inside
`fortran/setup.py`; the exact `tox` environment names.

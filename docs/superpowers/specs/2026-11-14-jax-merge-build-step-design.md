# Differentiable kernels in climt: merging the JAX work and adding a build step — design

Date: 2026-11-14. Size: project. Branch: `feature/jax-merge-build-step`.

## Intent

The `jax-research` branch holds two differentiable kernels that work: the
Emanuel convection scheme and the CORK longwave radiation scheme, each
agreeing with the numba kernel it was derived from. The branch is 23
commits ahead of `develop` and 294 behind, so every week it stays out it
gets harder to bring in, and nobody can build on it.

Three things stand between that branch and something usable:

1. **The route into the kernels.** The branch's `JaxBackend` is a sympl
   state backend that returns arrays without converting units. A state in
   hPa handed to a component that wants Pa would be silently wrong. It is
   not exported and should not be merged.
2. **Compile time.** Measured on JAX 0.10.1, one column, 30 levels, CPU:
   the Emanuel forward kernel compiles in 18–29 s and runs in 2.4 ms;
   `jax.grad` through it compiles in **834 s** and runs in 7.6 ms. At 14
   minutes per compile, differentiating the scheme is not something a
   person iterates on. The likely cause (not yet verified) is that most
   per-level loops are Python loops inside the jitted function, each
   writing with `.at[i].set`, which unroll into a very large graph.
3. **Import hygiene.** On the branch, `import climt` imports JAX when it is
   installed and switches it globally to 64-bit.

sympl's component call does unit conversion, dimension matching and
checking on every call, through pint and xarray. None of that can run under
`jax.jit`, and none of it needs to: units and dimensions do not change
during a run. The **build step** does that work once, from the
components' own property dictionaries and one example state, and returns a
pure function from raw arrays to raw arrays. That function is the route
into the JAX kernels, and it is also a faster way to call the numba ones.

## Goals

- The Emanuel and CORK longwave JAX kernels, with their tests, on
  `develop`, with the branch's history preserved.
- numba stays the reference implementation. A JAX kernel is a derived
  artefact: it is accepted only while a parity test against the numba
  kernel passes.
- `import climt` does not import JAX. JAX is an optional extra.
- A compile-time budget for the Emanuel kernel, enforced by tests.
- `climt.compile_component(component, state)` and
  `climt.compile_tendencies(components, state)`: resolve units, dimensions
  and aliases once and return a pure function over a dictionary of raw
  arrays, which works with NumPy arrays and, for the two components that
  have JAX kernels, with JAX arrays under `jit`, `grad`, `vmap` and
  `lax.scan`.
- One worked example: a radiative–convective column (CORK longwave plus
  Emanuel) stepped inside `lax.scan`, with a gradient of the outgoing
  longwave radiation checked against finite differences.

## Non-goals

- `JaxBackend`, or any way to pass JAX arrays through `component(state)`.
  The compiled function is the JAX entry point.
- JAX kernels for any other component (shortwave, surface, boundary
  layer). Compiling them works with NumPy arrays only.
- Compiling `Stepper` components or sympl's wrappers
  (`UpdateFrequencyWrapper`, `TimeDifferencingWrapper`, `ScalingWrapper`),
  tracers, or sympl's time steppers. The example steps with forward Euler.
- Coupling to the JAX dynamical core in `gfs_dynamical_core`.
- GPU or Metal. CPU only, 64-bit only.
- The research notes on the branch (`JAX_KERNEL_FUSION.md`,
  `PROPOSAL_JAX_BACKEND.md`, `docs/HYBRID_STRATEGY.md`,
  `docs/SENSITIVITY_ANALYSIS_STRATEGY.md`, `docs/predict_parameters*.py`).
  They stay on `jax-research`, which is not deleted.

## Anti-goals

- A compiled function whose numbers differ from `component(state)`. sympl's
  own call is the oracle, for every component the build step accepts.
- A silent fall-back from JAX to NumPy. If JAX arrays reach a component or
  a configuration that has no JAX kernel, that is an error naming the
  component and the reason, because the alternative is a gradient that is
  silently zero or a tracer leaking into NumPy.
- climt changing JAX's global configuration. If 64-bit mode is off, the
  kernels refuse to run and say how to turn it on.
- A budget met by weakening parity. The parity tolerances are fixed first.

## Approaches

For the route into the kernels:

1. **A build step that returns a pure function** (chosen). Units and
   dimensions are resolved once with sympl's own machinery, so they are
   right by construction and tested against sympl. The result is an
   ordinary function, which is what `jit`, `grad` and `scan` want.
2. **Fix `JaxBackend` so it converts units.** Keeps `component(state)`
   working with JAX arrays, but pint cannot run on tracers, so conversion
   would have to be re-implemented inside the backend anyway, and the
   checking and xarray wrapping would still run on every traced call.
3. **Call `array_call` by hand**, as the branch's tests do. No new code,
   but every user re-derives each component's units, dimension order and
   aliases, which is exactly the mistake the unit system exists to prevent.

For bringing the branch in:

1. **A real merge of `jax-research`, then clean up in follow-up commits**
   (chosen). History and authorship are kept; the one conflict
   (`pure_python_v3.py`) is small.
2. **Copy the files across.** Cleaner diff, but loses 23 commits of
   history, including the block-by-block record of how the Emanuel port
   was made differentiable.

## Design

### What is merged

`git merge origin/jax-research`. Kept: the JAX code in
`climt/_components/emanuel/pure_python_v3.py`, the JAX dispatch in
`climt/_components/cork/lw/component.py`,
`climt/_components/cork/lw/kernels_jax.py`,
`climt/_components/cork/optics/correlated_k_jax.py`, the five JAX test
files, and the three plan and design documents under `docs/superpowers`.
Removed in the merge commit: `climt/_core/jax_backend.py` and the research
notes listed under Non-goals.

In the conflict, `develop`'s side wins for the NumPy path (it coerces the
timestep to a plain `float`, which the Pyodide build needs) and the
branch's side supplies the JAX path.

### JAX isolation

- The Emanuel JAX kernel moves out of `pure_python_v3.py` into
  `climt/_components/emanuel/kernel_jax.py`. No module that `import climt`
  loads imports JAX at module level.
- `climt/_core/jax_support.py` holds the two questions components ask:
  *is this a JAX array?* (answered without importing JAX: if JAX has not
  been imported, it cannot be) and *is 64-bit mode on?* (if not, raise
  `JaxPrecisionError` with the one line that fixes it).
- A component with a JAX kernel dispatches on the type of its first array
  input. CORK longwave builds its JAX table on first use, not in
  `__init__`. If JAX arrays arrive and the configuration is one the JAX
  path does not cover (optics other than `correlated_k`, a table that is
  not the seven-dimensional pre-mixed kind, `diagnostics_level` above 0),
  it raises `NotImplementedError` saying which.
- A component without a JAX kernel that is handed JAX arrays through a
  compiled function raises `TypeError` naming the component.
- `jax` is an optional extra, `climt[jax]`, pinned at the versions the
  measurements were made with (`jax>=0.10.1`, `jaxlib>=0.10.1`). Tests that
  need JAX are marked `jax` and skipped when it is absent; CI runs them in
  one extra job (Linux, Python 3.12).

### Parity policy

Each JAX kernel has a parity test against the numba kernel on the same
inputs: Emanuel at `rtol = atol = 1e-10` over moist, dry and unstable
profiles and one and four columns; CORK longwave at `rtol = 1e-6`. A change
to a numba kernel that breaks parity fails CI until the JAX kernel is
brought back in step or explicitly retired.

### Compile-time budget

Two measures, because wall time is the thing that matters and graph size
is the thing that can be tested deterministically:

- **Graph size**: the number of equations in the traced program
  (`jax.make_jaxpr`), counted through nested sub-programs. It does not
  depend on the machine. The test suite pins it.
- **Wall time** of the first call, for the forward kernel and for
  `jax.grad`, at 30 levels and one column on CPU.

Targets for the Emanuel kernel on the maintainer's machine with JAX
0.10.1: **gradient compile at most 120 s** (from 834 s) and **forward
compile at most 30 s** (no worse than today). The work is: measure first;
then replace per-level Python loops with whole-array masked expressions
where each level is independent, with `lax.scan` or a cumulative operation
where a level depends on the one before, and with broadcasting where two
level loops are nested. Parity is re-run after each block.

If the targets cannot be reached this way, the numbers achieved and the
largest remaining contributors are reported and the target is not quietly
moved.

### The build step

```python
compiled = climt.compile_component(component, state)
raw = compiled.extract(state)                 # {name: ndarray}, canonical form
tendencies, diagnostics = compiled(raw, dt)   # pure; dt in seconds, where the component takes a timestep
```

**Canonical form** is the example state's own: each quantity keeps the
units and dimension order it has in the state passed to
`compile_component`. The compiled function accepts and returns arrays in
that form. So a quantity has one representation however many components
use it, and the conversion to and from each component's preferred units
and dimension order lives inside the compiled function as a multiply, an
add, a transpose and a reshape.

At build time, for each input the build step records: the scale and offset
from the state's units to the component's (obtained from sympl's unit
registry, and verified to be affine); the transpose, broadcast and
wildcard-flattening that take the state's dimensions to the component's;
and the alias, if the component uses one. For each output it records the
inverse. Tendencies come back in the state's units per second; a
diagnostic that is not in the state comes back in the units and dimensions
the component declares.

`compile_tendencies(components, state)` does the same for a list of
`TendencyComponent`, `ImplicitTendencyComponent` and `DiagnosticComponent`
objects: every component sees the same input arrays, tendencies of the
same quantity are summed, and diagnostics are merged. Two components
producing the same diagnostic is an error at build time.

`compiled.restore(tendencies, diagnostics)` wraps raw outputs back into
sympl quantities, for getting results out at the end of a run.

Errors at build time, not at call time: a component that is a `Stepper` or
a wrapper; an input missing from the state; units that cannot be converted;
a state quantity in an offset unit such as `degC` whose tendency is
requested.

### Example and documentation

`examples/differentiable_column.py` and a user-guide page,
`docs/user-guide/differentiable.qmd`: what has a JAX kernel, how to compile
a model, the 64-bit requirement, the compile times measured, and the
parity policy.

## Constraints

- `import climt` must not import `jax`, with or without JAX installed.
- The whole existing suite passes without JAX installed.
- Works with both state containers (`sympl.DataArray`,
  `UnytStateContainer`): `extract` goes through the active sympl backend.
- The compiled function itself touches no pint, xarray or unyt object.
- Python 3.10 remains supported for everything except the JAX extra.

## Left to the builder

How the Emanuel loops are regrouped, block by block, within the rules
above; the layout of the benchmark's output; the names of private helpers;
the plots, if any, in the example.

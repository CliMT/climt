# A compiled single-column run: the build step — design

Date: 2026-11-19. Size: project. Branch: `feature/scm-build`.

Two smaller projects come first and are designed separately:
`2026-11-19-numba-ports-design.md` (simple_physics and dcmip to numba) and
`2026-11-19-cffi-fortran-design.md` (climt-fortran 2.0: Cython replaced by
cffi). This design assumes both are done: every component in a model is
either a numba kernel or a cffi call into Fortran.

## Intent

Most climt users run single-column models, and many want to run many of
them: a sweep over CO₂, a set of initial profiles, a range of a convection
parameter. The aim is to push that use case to its performance limit while
the sympl components stay the readable, unit-checked way to write and
understand a model.

Moving from pint to unyt removed most of sympl's overhead, but not all of it.
Measured on 2026-11-19 (climt-dev env, `develop` at 065c5e94, the GMD
radiative–convective column of `examples/gmd_radiative_convective_unyt.py`:
Fortran Emanuel, RRTMG longwave and shortwave, slab surface, SimplePhysics;
28 levels, one column, Adams–Bashforth, 5-minute step, 500 steps after
warm-up):

| Where the time goes | ms per step |
|---|---|
| sympl and unyt bookkeeping (wildcard matching, `get_array`, unit conversion) | 2.2 |
| SimplePhysics re-parsing its 16 constants with pint on every call | 1.2 |
| Kernels: RRTMG SW 0.23, LW 0.12, SimplePhysics Fortran ~0.15, Emanuel 0.05, slab 0.04 | ~0.6 |
| **Total** | **4.0** |

The DataArray backend takes 14.4 ms for the same step. The kernels are about
15% of the unyt step, so removing everything else bounds the gain on one
column at about 7×; beyond that only faster kernels help.

A second measurement shapes the design: the numba Emanuel kernel
(`_numpy_vectorized_convect`, `parallel=True`) takes 0.74 ms per call on one
column with numba's default 11 threads and 0.16 ms with one thread. Starting
a parallel region costs more than the column's physics. Thread count has to
follow the number of columns.

## The approach

Three ways to build the compiled run were considered:

- **A. Generate one numba function** for the whole run (all components, the
  time stepper and the time loop), parallel over columns once.
  Reaches the limit; needs a code generator.
- **B. One fixed kernel signature** over a packed `(variable, column, level)`
  array, called from a generic numba loop. No code generation, but every
  kernel is rewritten into index lookups, against the readability the
  components exist for.
- **C. A Python loop over compiled kernels** (chosen). Units, dimensions and
  aliases are resolved once; each step calls the components' existing
  `array_call` on raw arrays from a short Python loop.

The user chose **C** to move conservatively: it reuses each component's
`array_call` and its kernels unchanged, and its result is easy to check
against sympl. The Python that remains per step is estimated at 10–20% of a
0.6 ms step. **A** is deferred, not rejected: this project ships a benchmark
so that A's gain over C can be measured later and decided on numbers.

This build step was planned as Tasks 5–6 of the jax-merge design
(`feature/jax-merge-build-step`). The user moved it here so the performance
work does not wait on jax-merge's open-ended Task 4. jax-merge keeps Tasks
1–4; its Tasks 5–7 are removed from that plan, and its JAX kernels plug into
this build step afterwards.

## Goals

- `climt.build_column_model(...)`: from components, an example state and a
  timestep, build once, then run many steps with no sympl, unyt or pint
  object touched inside the loop.
- Tendency components (`TendencyComponent`, `ImplicitTendencyComponent`),
  `DiagnosticComponent`s and `Stepper`s, in the order a sympl loop applies  them; forward Euler and Adams–Bashforth (orders 1–4, default 3, starting at
  lower order exactly as sympl's `AdamsBashforth` does).
- Many scenarios in one run, one per column. A scenario can change state
  inputs (any quantity in the state) and component parameters.
- Thread count chosen from the number of columns, so one column does not pay
  for a parallel region.
- Numbers that match the equivalent sympl loop.
- A benchmark, committed, that times the sympl loop and the compiled run on
  the GMD column at several column counts and thread counts, and writes its
  results in a form a later comparison with approach A can reuse.

## Non-goals

- Approach A (one generated numba function). Decided later from the
  benchmark.
- JAX arrays and gradients. They arrive through jax-merge later; nothing here
  may make that harder (the call path stays array-namespace-agnostic where it
  costs nothing), but nothing here tests JAX.
- sympl's wrappers (`UpdateFrequencyWrapper`, `TimeDifferencingWrapper`,
  `ScalingWrapper`), leapfrog and Runge–Kutta steppers, tracers. Each is a
  build-time error naming what is unsupported.
- Models with horizontal coupling (dynamical cores, GCMs). Columns are
  independent here.
- Making kernels themselves faster. The measured kernel times are the floor
  this project works down to.

## Anti-goals

- A compiled run whose numbers differ from the sympl loop it replaces. The
  sympl loop is the oracle for every model the build step accepts.
- A silent fallback. Anything the build step cannot compile (a wrapper, a
  per-scenario parameter a component cannot take per column, an unsupported
  stepper) is an error at build time naming the component and the reason.
- A parameter sweep that silently becomes one call per scenario. If a
  configuration would make the per-step cost grow with the number of
  scenarios, the build step says so.
- Physics duplicated between the readable and the compiled path. Both call
  the same kernels; only the glue around them differs.

## Design

### Calling it

```python
model = climt.build_column_model(    state,                                   # example state: one column, or already n columns
    diagnostic_components=[],                # run first, on state n
    tendency_components=[convection, radiation_lw, radiation_sw, slab],
    steppers=[simple_physics],               # applied in order after the time step
    timestep=timedelta(minutes=5),
    time_stepping="adams_bashforth",         # or "euler"; order=3 by default, as sympl
    scenarios=climt.Scenarios(
        n=256,
        inputs={"surface_temperature": sst_values},           # (n,) or (n, ...) arrays, state units
        parameters={(convection, "alpha"): alpha_values},     # (n,) arrays
    ),
    prescribed=["eastward_wind"],            # reset to the scenario's value every step
    threads=None,                            # None: chosen from the number of columns
)
result = model.run(n_steps=20000, output_every=20,
                   outputs=["air_temperature", "surface_temperature"])
result.to_xarray()                           # dims (time, scenario, ...)
model.to_state(scenario=3)                   # final state of one scenario, as sympl quantities
```

One step does what the sympl loop in
`examples/gmd_radiative_convective_unyt.py` does: diagnostic components on
state *n*, their outputs merged into it; tendency components on state *n*
(each sees the same input; tendencies of one quantity add;
diagnostics merge, a name produced twice is a build-time error), the time
stepper's update, then each stepper in order on the updated state, then the
prescribed quantities reset and time advanced. Diagnostics are written into
the state as the sympl loop's `state.update(diagnostics)` does. `time` is
kept as a Python `datetime` outside the arrays and passed to any component
whose inputs include it.

### Building

At build time, from each component's own property dictionaries and the
example state, the build step records once what sympl redoes on every call:

- **Canonical form.** Each quantity keeps the units and dimension order it
  has in the example state. The run holds arrays in that form.
- **Per input:** the affine map from state units to the component's units
  (from the unit registry, checked to be affine), the transpose and reshape
  to the component's dimension order with the wildcard flattened, and the
  alias the component uses.
- **Per output:** the inverse map. Tendencies come back in state units per
  second; diagnostics not in the state in the units and dimensions the
  component declares.
- Unit factors of 1 and offsets of 0 are skipped at call time, so most
  quantities pass through untouched.

The jax-merge plan's Task 5 (`compile_component`, `ArraySpec`,
`CompileError`, the 0/1/2 affine check, the dimension plan) is the
specification for this part and is taken over here, without the JAX
dispatch.

### Columns and scenarios

The run holds `n` columns. If the example state has one column, it is
repeated `n` times; per-scenario inputs then overwrite their column. A
per-scenario input may be one value per scenario (broadcast over its other
dimensions) or a full array per scenario.

### Per-scenario parameters

A component parameter can vary by scenario only if the component's kernel
takes it per column. A component declares which parameters it can take per
column; for those, its kernel receives an `(ncol,)` array and indexes it by
column. The readable path passes the scalar, broadcast. The build step works
on a shallow copy of each component whose parameter record (`_params`, a
`NamedTuple`) holds the arrays; the user's component objects are not
modified. Asking for a
per-scenario value of a parameter the component does not declare is a
build-time error naming the component, the parameter and the ones it does
take.

In this project the components of the benchmark model are covered: the
numba SimplePhysics written per column from the start (numba-ports design),
the numba Emanuel's `EmanuelParams` (`EmanuelConvectionPython`; the Fortran
`EmanuelConvection` sets its parameters as Fortran globals), and the slab
surface (its mixed-layer depth is already a state input). Changing the numba
Emanuel kernel's parameters to arrays touches the kernel jax-merge derives
its JAX Emanuel from; jax-merge's parity test must keep passing after this
lands, with its JAX kernel taking scalars. Fortran components' parameters are Fortran
globals or one value per call, so they are fixed for a build, except RRTMG's
per-call settings (below).

**RRTMG.** Solar constant, day of year, the cloud and ice property options
and the overlap method are one value per call; physical constants and the
dry-air heat capacity in `rrtmg_lw_ini` are Fortran globals set once. If
scenarios vary a per-call setting, the build step sorts columns into groups
that share it and calls RRTMG once per group per step, and warns at build time with the
number of groups, since each group is a separate call. Varying a global is a
build-time error.

### Threads

Kernels that use `prange` run with a thread count the build step sets with
`numba.set_num_threads` for the duration of `run` and restores afterwards:
one thread below a column threshold, more above it. The threshold is taken
from the benchmark, not guessed. The user can pass `threads` to override.
RRTMG is called serially (one call per group per step), so its
thread-safety is not relied on.

### Output

`run` keeps the named quantities every `output_every` steps in preallocated
arrays (canonical form) and returns them with their specs. `to_xarray()`
wraps them with a `scenario` dimension and units as attributes. The model
keeps its state between `run` calls, so a run can be continued.

### Errors at build time

A wrapper or an unsupported stepper; an input missing from the state; units
that cannot be converted affinely; a tendency requested for a quantity in an
offset unit (`degC`); a diagnostic produced twice; a per-scenario parameter
the component cannot take per column; a per-scenario RRTMG global.

### Testing

- **Against sympl.** For every compilable component, the compiled call
  matches `component(state)` (taken from jax-merge Task 5's test). For the
  GMD column and a dry variant, a compiled run of 50 steps matches the sympl
  loop step by step, at one column and at several columns with different
  inputs and parameters (each column compared to a sympl run of that
  scenario alone).
- **Per-column parameters.** A kernel given an `(ncol,)` parameter array
  matches the same kernel run column by column with scalars.
- **Errors.** Each build-time error above, with its message.

### Benchmark

`benchmarks/scm_build.py`: the GMD column (with the Fortran Emanuel as in the
example, and with the numba Emanuel), sympl+unyt loop against the
compiled run, at 1, 16, 256 and 4096 columns and one thread versus all
cores, reporting milliseconds per step per column and the share spent
outside kernels. Results are written as JSON with the machine, versions and
commit, so approach A can later be measured with the same script and
compared.

## Constraints

- `import climt` stays light: no new import-time dependencies.
- The whole existing suite passes; the sympl path of every component is
  unchanged except where a kernel gains per-column parameters, where the
  scalar case gives identical numbers.
- Works with both state containers (`sympl.DataArray`,
  `UnytStateContainer`).
- Python 3.10–3.12, as climt.

## Left to the builder

Private helper names; the exact column threshold for threads (from the
benchmark); the tolerance of the step-by-step comparison (tight, and the
same for all models); the layout of the benchmark's JSON and printed table;
how `Scenarios` validates shapes.

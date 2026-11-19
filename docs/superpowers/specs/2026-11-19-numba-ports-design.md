# simple_physics and dcmip in numba — short design

Date: 2026-11-19. Size: small change (two ports). Branch: `feature/scm-build`.
First of three projects; see `2026-11-19-scm-build-design.md` for why.

## Intent and why

The build step (`2026-11-19-scm-build-design.md`) works on numba kernels and
cffi calls. Two of climt-fortran's five libraries are small enough to port:
`simple_physics` (~1,000 lines, Reed and Jablonowski 2012, runs every step)
and `dcmip` (~2,100 lines, initial conditions for DCMIP tests 4 and 5, runs
once). Once ported they leave climt-fortran, so less Fortran goes through the
cffi switch. Porting SimplePhysics also removes the pint re-parsing of its
constants on every call, measured at 1.2 ms of a 4.0 ms single-column step.

## Goals

- `SimplePhysics` and `DcmipInitialConditions` run on numba kernels, with the
  same inputs, outputs, properties and constructor arguments as today.
- The SimplePhysics kernel takes its parameters (drag coefficients, boundary
  layer top and decay, the physical constants it uses) as per-column arrays,
  so the build step can vary them by scenario. Scalars broadcast.
- Constants are read once, at construction, never per call.
- After both ports land in climt, a later climt-fortran release drops
  `_simple_physics` and `_dcmip` (in the cffi project).

## Non-goals and anti-goals

- No change to the physics. Numbers match the Fortran within the tolerances
  below; a discrepancy is a bug in the port, not a chance to fix the scheme.
- Emanuel is not ported here. A numba Emanuel exists and matches the Fortran,
  but at one column and one thread it takes 0.16 ms per call against the
  Fortran's 0.05 ms, so the Fortran stays (and moves to cffi) until a numba
  version is as fast. Speeding it up is not part of this work.
- Not a fallback: once ported, the components no longer need climt-fortran.

## Approach

Line-by-line ports, one numba kernel per scheme, looping over columns with
`prange` and over levels inside. Alternatives considered: wrapping the
Fortran through cffi like RRTMG (keeps Fortran the user wants gone);
vectorised NumPy (slower than numba for per-level loops at this size).

## Checks

- **Parity.** Before the Fortran is removed, reference outputs are generated
  from the Fortran over a set of inputs (moist and dry profiles, tropical
  cyclone and baroclinic-wave cases for dcmip, several columns, each switch of
  SimplePhysics on and off) and committed as test fixtures. The numba ports
  must match them at `rtol = 1e-10` for SimplePhysics and `1e-12` for dcmip
  (pure arithmetic). The fixtures keep the comparison alive after the Fortran
  is gone.
- **Speed (SimplePhysics only).** At one column and one thread the numba
  SimplePhysics call is no slower than 1.5× the Fortran call, measured
  through `array_call` with the constant re-parsing excluded from the Fortran
  number. Recorded, with the machine, in the benchmark output.
- **Per-column parameters.** An `(ncol,)` parameter array gives the same
  numbers as running each column alone with its scalar.

## Left to the builder

Kernel structure and helper names; the exact fixture inputs within the cases
above; how the fixtures are generated (a script committed beside them).

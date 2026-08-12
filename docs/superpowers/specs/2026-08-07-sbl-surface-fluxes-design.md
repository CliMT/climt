# SimpleBoundaryLayer surface fluxes — design

**Date:** 2026-08-07
**Status:** approved, ready for planning

## Problem

`SimpleBoundaryLayer` (SBL) diffuses heat, moisture and momentum but never
injects anything at the surface. Its `_diffuse_profile` uses no-flux boundaries
(`diag_m[0] = diag_p[n-1] = 0`), so `surface_temperature` and
`surface_specific_humidity` influence only the *diffusivity* (via the
surface-layer Richardson number `Ri_a` → exchange coefficient `C`) and the
wind-stress *diagnostic*. The stress is computed and discarded.

Verified empirically. With a 320 K saturated surface under a 250 K, bone-dry
atmosphere, run 20 × 1200 s:

```
column enthalpy  2.591876e+09 -> 2.591876e+09  J/m2   change +0.0000e+00
column water     1.032620e-02 -> 1.032620e-02  kg/m2  change +0.0000e+00
lowest-level T   250.000 K   (surface is 320 K)
wind stress diag 0.11410 Pa
lowest-level u   10.000 m/s (started at 10.0)
```

The nonzero stress confirms `C` and the diffusivity are active; nothing reaches
the atmosphere.

This is not a regression. The behaviour was deliberate and documented in the
original spec (`2026-07-22-boundary-layer-and-jit-tridiagonal-design.md`, lines
198–203) and the user guide, and the upstream reference implementation
(`~/github/other_libs/component_abel.py`, byte-identical to `component.py`) has
the same no-flux solve and the same "assumes a surface flux component has been
already run" caveat. Frierson et al. (2006) the *paper* has both halves; the
thesis code split them and only the diffusion half was ever implemented.

### Consequences

- A grey model built from SBL + `GrayLongwaveRadiation` + `EmanuelConvection`
  has **no energy source at all** — radiation only cools, convection only
  redistributes — and cools monotonically.
- No surface drag reaches the flow, so there is no Ekman friction.
- `BucketHydrology`, `SlabSurface` and `IceSheet` emit
  `surface_upward_sensible_heat_flux` / `surface_upward_latent_heat_flux` as
  **diagnostics only** (`output_properties` carries just `surface_temperature`
  and, for the bucket, soil moisture). Those fluxes reach nothing.
  `examples/full_gcm_land_ocean_ice.py` pairs `BucketHydrology` with SBL and no
  `SimplePhysics`, so it has the same gap.

## Goals

1. Apply Frierson (2006) bulk surface fluxes of heat, moisture and momentum,
   reusing the `C`, `wind_int[0]` and `col_rho[0]` already in the kernel.
2. Alternatively apply **externally supplied** surface heat/moisture fluxes, so
   surface components' existing flux diagnostics reach the atmosphere.
3. Preserve today's no-flux behaviour behind an explicit opt-out.
4. Make the wind-stress diagnostics report the momentum flux the scheme
   actually applies, rather than a separately-derived number.

## Non-goals

- Switching heat diffusion from `T` to dry static energy or potential
  temperature. Frierson uses DSE; diffusing `T` ignores the adiabatic gradient.
  The existing scheme already diffuses `T`, so a `T`-based surface flux is at
  least self-consistent. Changing the diffused variable affects the whole
  scheme and is separate work.
- Evaporative resistance / β for land. `surface_specific_humidity` is already an
  input; land dryness stays the caller's responsibility.
- Separate exchange coefficients for momentum, heat and moisture. Frierson uses
  one `C`; keep one.
- Changing `_richardson_diffusivity`. Its use of `Ri_a` (rather than the
  reference implementation's local `Rich[i]`) in the stability denominator is
  **settled as correct** and stays: `Ri_a` is what makes `K_b` continuous at
  `Ri_a == 0`. See Resolved questions.

## Design

### Boundary condition in the tridiagonal

Layer 0 has mass `Δp₀/g` per unit area, `Δp₀ = p_int[0] − p_int[1] > 0`. With
`F_surf` the upward flux of `X` into the layer:

```
(Δp₀/g) dX₀/dt = F_surf − F_top
```

**Bulk.** `F_surf = ρ_s C |v| (X_s − X₀)` depends on the unknown, so it is taken
implicitly. With `β = g·ρ_s·C·|v|·dt / Δp₀`:

```
diag[0] += β
rhs[0]  += β · X_s
```

`X_s` is `T_s` for temperature, `q_s` for humidity, and `0` for both wind
components (no-slip), which turns the existing stress calculation into actual
drag on the flow.

**External.** `F` is known, so it is a pure Neumann source:

```
rhs[0] += g·dt·F / Δp₀        (diag unchanged)
```

with, for `surface_upward_sensible_heat_flux` `SH` and
`surface_upward_latent_heat_flux` `LH`, both W m⁻² positive upward (into the
atmosphere, hence positive sign):

```
rhs_T[0] += g·dt·SH / (cp·Δp₀)
rhs_q[0] += g·dt·LH / (L_v·Δp₀)      L_v = latent_heat_of_condensation = 2.5e6
```

Momentum stays bulk-internal in external mode: surface components do not
compute wind stress, and `northward_wind_stress` / `eastward_wind_stress` are
already SBL *diagnostics*, so accepting them as inputs would collide.

`ρ_s` is `col_rho[0]` — the density at the first interface — matching what the
existing stress diagnostic already uses.

### Surface-flux diagnostics

**One rule for all of them: a surface-flux diagnostic reports the flux the
solve actually delivered, evaluated with the post-solve layer-0 value.**
Because the bulk flux is implicit, that value is `X₀ⁿᵉʷ`, not `X₀ᵒˡᵈ`. With
`κ = ρ_s·C·|v|` (`col_rho[0] * C * wind_int[0]`, one quantity shared by all
four):

```
surface_upward_sensible_heat_flux = cp · κ · (T_s − T₀ⁿᵉʷ)
surface_upward_latent_heat_flux   = L_v · κ · (q_s − q₀ⁿᵉʷ)
eastward_wind_stress              = κ · u₀ⁿᵉʷ
northward_wind_stress             = κ · v₀ⁿᵉʷ
```

The heat and moisture diagnostics are new, and belong to `'bulk'` mode only —
in `'external'` mode those two names are *inputs*, so they must not also be
diagnostics. They are what makes the budget and bulk-vs-external consistency
tests exact rather than approximate: substituting the diagnosed `F` into the
external mode's Neumann term reproduces the bulk system algebraically.

**This changes the wind-stress diagnostics.** Today they are
`κ · u_int[0]` — the *pre-step interface* wind, `0.5·(u₀ᵒˡᵈ + u₁ᵒˡᵈ)`. The drag
the solve applies is `κ · u₀ⁿᵉʷ`, so the old diagnostic never matched the
momentum the column actually lost. The reconciled form closes the momentum
budget the same way the heat and moisture ones close theirs:

```
(1/g) Σ Δp_k (u_k ⁿᵉʷ − u_k ᵒˡᵈ) = −eastward_wind_stress · dt
```

The formula is uniform across all three modes. In `surface_fluxes=None` no drag
is applied, so the stress there is advisory only — the bulk stress a separate
surface component could apply — and its *value* still changes relative to
today, because `u₀ⁿᵉʷ ≠ u_int[0]`. That is accepted: a single definition beats
a mode-dependent one, and no in-tree consumer reads these diagnostics.

### Implementation

`_diffuse_profile` gains two optional scalars, `surface_exchange` (the `β`
term, zero in external mode) and `surface_source` (the Neumann term, zero in
bulk mode). One code path serves all three modes; today's behaviour is both
scalars zero.

### API

```python
SimpleBoundaryLayer(surface_fluxes='bulk')     # default
SimpleBoundaryLayer(surface_fluxes='external')
SimpleBoundaryLayer(surface_fluxes=None)       # today's no-flux behaviour
```

`'bulk'` is the default so the component matches its own name and docstring.
This is a behaviour change for existing callers, accepted deliberately: the
component is one commit old and the current default is a footgun.

`input_properties` is built in `__init__` — `'external'` adds
`surface_upward_sensible_heat_flux` and `surface_upward_latent_heat_flux` as
required inputs; the other modes must not declare them.

Double-counting is the caller's responsibility and must be documented: pairing
`surface_fluxes='bulk'` with `SimplePhysics(surface_fluxes=True)` applies
surface fluxes twice.

In `'external'` mode the surface component computed its fluxes with its own
drag coefficient (`BucketHydrology` uses `bulk_coefficient=0.0011`), which
differs from SBL's Monin–Obukhov `C`. That inconsistency is inherent to the
modular split and must be documented, not silently reconciled.

## Tests

- **`surface_fluxes=None` conserves exactly.** Keep the existing
  `test_simple_boundary_layer_conserves_column_integrals` body, pinned to this
  mode.
- **Budget test (bulk).** Change in column enthalpy and column water equals the
  time-integrated surface flux to tight tolerance. Strictly stronger than the
  current conservation test.
- **Budget test (external).** Prescribed `SH`/`LH` produce exactly the implied
  column enthalpy and water change.
- **Budget test (momentum).** Change in column momentum equals
  `−wind_stress · dt` for both components, in bulk and external mode. This is
  the test the old stress diagnostic could not have passed.
- **Sign tests.** A warm, moist surface warms and moistens the column; drag
  decelerates the flow. The empirical case in Problem above becomes a
  regression test with the assertions inverted.
- **Property test.** `input_properties` gains the two flux fields in
  `'external'` mode only, and `diagnostic_properties` gains them in `'bulk'`
  mode only; `test_simple_boundary_layer_properties` gets a per-mode
  parametrisation.
- **Consistency.** Bulk mode and external mode fed the bulk-diagnosed fluxes
  agree to solver tolerance.

## Impact

- `tests/test_components.py:809` must be rewritten (see Tests).
- `northward_wind_stress` / `eastward_wind_stress` change value in **every**
  mode, including `surface_fluxes=None` (see Surface-flux diagnostics). No
  in-tree code reads them, so this is a diagnostic-only change, but it must be
  called out in the user guide and the release notes.
- `docs/user-guide/simple_boundary_layer.qmd` — the "assumes a surface-flux
  component has already run" framing (lines 15–16, 78–79) and the conservation
  claim (78–79) become mode-dependent. The class docstring's contradiction
  ("This is the surface-flux / boundary-layer formulation…" followed by "assumes
  a surface-flux component has already run") is resolved by this work.
- `examples/full_gcm_land_ocean_ice.py` should move to
  `surface_fluxes='external'` so `BucketHydrology`'s fluxes drive the
  atmosphere, rather than SBL computing a second, inconsistent set.

## Component ordering and the external-mode lag

Resolved by inspection of `sympl.ImplicitTendencyComponentComposite.__call__`
(`sympl/_core/composite.py`): it calls **every component with the same input
`state`** and merges diagnostics only into its *return* dict, never back into
`state`. Diagnostics therefore do not propagate between components within a
step. Two consequences for external mode:

- **SBL inside the dycore** (wrapped in `TimeDifferencingWrapper`, since
  `GFSDynamicalCore` accepts only `TendencyComponent` /
  `ImplicitTendencyComponent` and SBL is a `Stepper`): SBL reads the flux
  fields as they stood at the start of the step, i.e. the surface component's
  *previous*-step diagnostics, provided the caller does
  `state.update(diagnostics)` in the time loop. A one-step lag, and zero flux
  on step 0 from `get_default_state` defaults.
- **SBL called directly in the time loop**, after the surface component and
  with `state.update(...)` between: no lag.

Both patterns are legitimate; the docs must state the lag explicitly so callers
can choose.

## Resolved questions

1. **`Ri_a` vs local `Ri`.** climt's `_richardson_diffusivity` uses `Ri_a` in
   both the branch and the stability denominator; the reference `K_b` branches
   on `Ri_a` but uses local `Rich[i]` in the denominator. **Resolved: keep
   `Ri_a`.** It is the choice that makes `K_b` continuous at `Ri_a == 0`, which
   is the physically desirable behaviour and what the climt docstring already
   claims as the deliberate "thesis Eqn 2.8" modification. No code change, and
   no need to check the Abel thesis PDF to settle it.

## Open questions / risks

1. **Stability of implicit bulk drag at large `dt`.** `β` grows with `dt`; the
   implicit form should stay stable, but the GCM timestep (1200 s) should be
   exercised in tests.
2. **`TimeDifferencingWrapper` and implicitness.** Wrapping SBL converts its
   solve to `(new − old)/dt` tendencies that the dycore then sums with others.
   The surface flux BC stays implicit *within* SBL's own solve; the operator
   splitting against the other physics is unchanged from today. No action
   expected — noted so it is not mistaken for a regression.

## Adjacent bug (not part of this work)

`examples/full_gcm_land_ocean_ice.py:59-63` passes `BucketHydrology` and
`SimpleBoundaryLayer` to `GFSDynamicalCore` **unwrapped**. Both are `Stepper`s,
which the composite rejects:

```
TypeError: Component should be of type (TendencyComponent, ImplicitTendencyComponent)
           but is type BucketHydrology
```

`build_model()` therefore raises, and the example has never run despite its
docstring advertising a smoke test. (`UpdateFrequencyWrapper`, used on the
radiation components, preserves component type; `TimeDifferencingWrapper` is
what converts a `Stepper`.) The fix is to wrap both — worth doing alongside the
`surface_fluxes='external'` migration under Impact, but tracked separately.

# Land Surface Components Design

**Date:** 2026-07-20
**Branch:** feature/land-ocean-ice-components
**Goal:** Make climt's land surface physics more capable and complete: (1) a fully-featured **BucketHydrology** with an optional two-layer (deep + shallow) mode to capture both subseasonal and seasonal soil memory, and (2) a modular, intermediate-complexity **SecondBEST** land surface model implementing the BEST scheme (Pitman et al.), structured so individual process representations can be swapped or extended without forking the component.

This is one of two companion design docs. The other, `2026-07-20-ocean-ice-design.md`, covers the slab ocean, data ocean, and the sea-ice / land-ice split.

---

## Problem

climt's land surface today is thin:

| Component | What it does | Gaps |
|---|---|---|
| `BucketHydrology` (`_components/bucket_hydrology/component.py`) | Single-bucket Manabe hydrology: one moisture store (`lwe_thickness_of_soil_moisture_content`), β-limited evaporation, one thermal slab. | One store cannot separate a fast (subseasonal, precip/evap-driven) reservoir from a slow (seasonal) one. `soil_moisture_max` appears both as a config arg and a hardcoded `0.15` clamp (`component.py:216`). No sub-surface heat memory beyond a single slab. |
| `SecondBEST` (`_components/second_best/`) | Stub — `__init__.py` imports a `component.py` that does not exist; only `BEST-Description.pdf` is present. | No implementation at all. A full description exists in `docs/Description_of_SecondBEST.ipynb` and the PDF. |

We want land surface schemes that can support realistic multi-timescale behaviour and, together with the ocean/ice components, enable AMIP-style GCM runs.

---

## Design overview

Two independent workstreams, each shippable on its own:

1. **BucketHydrology** gains a backward-compatible `num_layers` option (`1` = today, `2` = shallow + deep). Single code path; the 1-layer default reproduces current results bit-for-bit.
2. **SecondBEST** is built new as a thin **orchestrator** `Stepper` that composes a set of **pluggable process objects** (soil properties, albedo, surface layer, fluxes, subsurface transport). Each process has a small documented interface with a BEST-default implementation, so changing or adding a process representation is a constructor argument, not a code fork.

Both follow climt/sympl conventions: `Stepper` subclasses, `input_properties` / `output_properties` / `diagnostic_properties`, `array_call`, state quantities registered in climt's state initialisation, and tests using the cached-result pattern.

---

# Part 1 — Two-layer BucketHydrology

## Section 1.1 — Configuration

```python
BucketHydrology(
    num_layers=1,                    # 1 (default, current behaviour) or 2
    soil_moisture_max=0.15,          # per-layer field capacity (m)
    beta_parameter=0.75,
    specific_latent_heat_of_water=2260000,
    bulk_coefficient=0.0011,
    # two-layer-only params (ignored when num_layers == 1):
    deep_soil_moisture_max=0.50,     # deep reservoir capacity (m)
    moisture_diffusion_timescale=None,   # s; shallow<->deep exchange (default derived below)
    deep_layer_thickness_ratio=10.0, # deep thermal/moisture layer thickness / shallow
    deep_drainage_timescale=None,    # s; loss from deep layer to groundwater (None = no drainage)
)
```

The hardcoded `0.15` clamp at `component.py:216` is replaced by `self._smax` (bug fix — the clamp currently ignores the configured `soil_moisture_max`).

## Section 1.2 — State representation

The single-layer path keeps the current scalar-per-column field for backward compatibility. The two-layer path adds a soil-layer axis.

| Quantity | 1-layer dims | 2-layer dims | units |
|---|---|---|---|
| `lwe_thickness_of_soil_moisture_content` | `["*"]` | `["soil_levels", "*"]` (shallow = index 0, deep = index 1) | m |
| `soil_temperature` (new, prognostic) | `["*"]` | `["soil_levels", "*"]` | degK |

`surface_temperature` remains the diagnosed skin/shallow temperature exposed to the atmosphere and radiation. In 2-layer mode it equals the shallow soil temperature; the deep soil temperature is an additional prognostic store coupled by conduction.

Registration: add `soil_levels` as a vertical coordinate for land columns in climt's state init, defaulting to the layer count from the component.

## Section 1.3 — Water balance (2-layer)

Let `w_s`, `w_d` be shallow/deep water paths (m), with capacities `S_s = soil_moisture_max`, `S_d = deep_soil_moisture_max`.

- **Precipitation** enters the shallow layer: `P = convective + stratiform`.
- **Evaporation** is drawn from the shallow layer only, with the existing β limiter keyed on shallow wetness:
  `β = min(1, w_s / (g · S_s))`, `E = β · E_pot`, `E_pot = c · |U| · (q*_surf − q_air)`.
- **Shallow → deep exchange** relaxes the two toward equal relative saturation on `moisture_diffusion_timescale` τ_m:
  `F_sd = (w_s/S_s − w_d/S_d) · S_eff / τ_m`, where `S_eff` normalises to a water depth.
  Default `τ_m` derived from layer thicknesses so it is O(days) — documented, tunable.
- **Deep drainage** (optional): `D = w_d / τ_drain` if `deep_drainage_timescale` set, else runoff-only at capacity.
- **Overflow**: water above capacity in either layer becomes runoff (diagnostic `runoff_rate`).

Update (explicit Euler, matching current component):
```
w_s ← clamp(w_s + (P − E − F_sd) · dt, 0, S_s)
w_d ← clamp(w_d + (F_sd − D) · dt, 0, S_d)
```

Single-layer mode collapses to the current update (no deep layer, no exchange).

## Section 1.4 — Thermal balance (2-layer)

Two heat stores of thickness `Δz_s` (shallow) and `Δz_d = deep_layer_thickness_ratio · Δz_s` (deep), coupled by conduction:

```
C_s dT_s/dt = Q_net − G_sd
C_d dT_d/dt = G_sd
G_sd = k_soil (T_s − T_d) / (0.5(Δz_s + Δz_d))
```

where `Q_net` is the current net surface heat flux (SW + LW − SHF − LHF), `C_x = ρ · c · Δz_x`. `k_soil` is a new (constant) soil thermal conductivity. This gives the surface a fast diurnal response while the deep store integrates the seasonal cycle — the stated motivation.

Single-layer mode keeps today's single-slab update exactly.

## Section 1.5 — Diagnostics & tests

New diagnostics: `runoff_rate` (m s⁻¹), `deep_soil_moisture_fraction`, `deep_soil_temperature` (2-layer only).

Tests:
- 1-layer regression: assert bit-for-bit agreement with current cached results (guards backward compatibility).
- 2-layer conservation: total water `w_s + w_d + ∫(E−P)dt + runoff + drainage` closes to machine precision over an integration.
- 2-layer timescale behaviour: forced with a seasonal + subseasonal signal, deep store variance is dominated by the seasonal band, shallow by the subseasonal band.

---

# Part 2 — SecondBEST (modular BEST land surface model)

## Section 2.1 — Architecture: orchestrator + pluggable processes

`SecondBEST` is a `Stepper` that owns state I/O and time-stepping, but delegates every physical calculation to a **process object**. Each process is a plain Python class with a narrow, documented `__call__` contract and a BEST-default implementation. The orchestrator wires them together each step.

```python
class SecondBEST(Stepper):
    def __init__(
        self,
        soil_type="clay",                 # "clay" | "sand"
        num_soil_layers=3,                # BEST default: upper/lower/base
        minimum_wind_speed=1.0,           # m/s, floor on |U| for drag/fluxes
        soil_properties=None,             # SoilProperties instance  (default: BestSoilProperties)
        albedo=None,                      # SurfaceAlbedo            (default: BestSurfaceAlbedo)
        surface_layer=None,               # SurfaceLayer             (default: BestSurfaceLayer)
        fluxes=None,                      # SurfaceFluxes            (default: BestSurfaceFluxes)
        subsurface=None,                  # SubsurfaceTransport      (default: BestSubsurfaceTransport)
        **kwargs,
    ):
        ...
```

Passing `None` selects the BEST default. Passing an instance swaps the representation. This mirrors how climt lets users recombine components — but at the sub-process level within one land model.

### Data flow per step

```
                   state (atmos lowest level, radiation, soil/snow state)
                                     │
              ┌──────────────────────┼───────────────────────────┐
              ▼                      ▼                            ▼
      SoilProperties          SurfaceLayer                  SurfaceAlbedo
   (porosity, X_FC, X_wilt)  (C_DN, ζ, Ri, C_Dm, C_Dh)   (α_SW, α_LW  -> diagnostics)
              │                      │
              └──────────┬───────────┘
                         ▼
                   SurfaceFluxes
             (SHF, LHF via β, E_usmax; momentum flux)
                         │
                         ▼
               SubsurfaceTransport
   (implicit heat+liquid+ice diffusion, freeze/melt Γ, Neumann bottom BC)
                         │
                         ▼
      new soil/snow temperature, X_w, X_i, snow depth  ->  output state
```

## Section 2.2 — Process interfaces

Each contract is a dataclass-style input/output. Interfaces are documented in the module and in a dedicated "Extending SecondBEST" doc section.

**`SoilProperties`** — pure function of soil type / area type, computed once:
```
__call__(soil_type, area_type) -> {
    colour, texture, porosity X_v,
    field_capacity X_FC, wilting_point X_wilt,
    B, K_H0, psi_0            # exfiltration / conductivity params
}
```
BEST default (`BestSoilProperties`): Eqs 4.10–4.12 from the notebook (porosity `0.6 − 0.03·tex`, `X_FC = (0.95 − 0.086·tex)·X_v`, `X_wilt = X_v − 0.03`, land-ice overrides).

**`SurfaceLayer`** — stability-dependent drag:
```
__call__(z_mid, z0, wind_speed, T_surf, T_air, area_type) -> {C_Dm, C_Dh, C_DN, Ri}
```
BEST default (`BestSurfaceLayer`): roughness length by surface, `C_DN = (κ/(ln z_m − ln z0))²`, `ζ = exp(−κ/√C_DN)`, `Ri = −(g z_m/(T_s U²))(T_u − T_s)`, and the unstable/stable `C_Dm`, `C_Dh` branches (Section 6 of the notebook), with `ε = 0.01` land / `1.0` ocean.

**`SurfaceAlbedo`**:
```
__call__(soil_props, wetness W_Lu, area_type) -> {alpha_sw, alpha_lw}
```
BEST default (`BestSurfaceAlbedo`): Eqs 5.5–5.8 (bare soil `α_SW = 0.10 + 0.1·clr + 0.06(1−W_Lu)`, `α_LW = 2α_SW`; land ice variants). Exposed as `surface_albedo_for_direct_shortwave` / `..._diffuse_shortwave` diagnostics matching RRTMG inputs.

**`SurfaceFluxes`**:
```
__call__(drag, atmos_state, soil_state, soil_props, timestep) -> {SHF, LHF, momentum_flux, evaporation}
```
BEST default (`BestSurfaceFluxes`): bulk `SHF = ρ_s C_pd U C_Dh (T_u − T_s)`, `LHF = L_v ρ_s U C_Dh β_u (q*_u − q_s)`, wetness factor `β_u` with the exfiltration limit `E_usmax = K_HD Θ_u^{0.5B+2} − K_H0 Θ_u^{2B+3}` and `Θ_u = (W_Lu − 0.01)/(1 − W_Fu)` (Section 8). Momentum flux `U_s = −ρ_s C_Dm U u`.

**`SubsurfaceTransport`** — the numerically heaviest piece, isolated behind the smallest interface:
```
__call__(profiles{T, X_w, X_i}, surface_flux_bc, timestep, dz) -> new profiles{T, X_w, X_i}
```
BEST default (`BestSubsurfaceTransport`): implicit solve of the three conservation laws from the notebook —
```
c_v T_t = G_z + L_f Γ + c_pw R T_z            (heat)
ρ_w (X_w)_t = (R − (1−f_i)E)_z − Γ            (liquid; Π transpiration = 0, no vegetation)
ρ_i (X_i)_t = −(f_i E)_z + Γ                  (ice)
Γ = (c_v/L_f)(T_f − T)/Δt                      (potential ice production)
```
Implicit (Crank–Nicolson-style, reusing the tridiagonal solver pattern already in `surface_ice.py`) for the diffusive density terms; explicit Euler for the freeze/melt and advective sources. Neumann bottom BC `()_z = 0`. The advective term `c_pw R T_z` sets a CFL limit `Δz/(c_pw R) > Δt/2`, checked at init and documented.

## Section 2.3 — climt component surface

**Inputs** (lowest atmospheric level + radiation + soil/snow state): `air_temperature`, `specific_humidity`, `northward_wind`, `eastward_wind`, `air_pressure`, downwelling/upwelling SW & LW fluxes, `area_type`, and prognostic soil/snow fields below.

**Prognostic state** (`output_properties`, all on a `soil_levels` dim except snow):
`soil_temperature`, `soil_liquid_water_content` (X_w), `soil_ice_content` (X_i), `surface_snow_thickness`, `surface_temperature`.

**Diagnostics**: `surface_upward_sensible_heat_flux`, `surface_upward_latent_heat_flux`, `evaporation_rate`, `surface_albedo_for_direct_shortwave`, `surface_albedo_for_diffuse_shortwave`, `surface_drag_coefficient_for_heat`, `surface_drag_coefficient_for_momentum`.

These flux/albedo diagnostics are deliberately the same names BucketHydrology/RRTMG use, so SecondBEST is a drop-in swap for BucketHydrology on land columns.

## Section 2.4 — "Extending SecondBEST" (documentation deliverable)

A dedicated docstring section + notebook cell showing:
- The five process interfaces and their contracts.
- A worked example: replacing `BestSurfaceFluxes` with a simple bucket-β flux, or `BestSurfaceLayer` with a constant-drag scheme, by passing an instance to the constructor.
- How to add a *new* process (e.g. a vegetation/transpiration term) by extending the orchestrator's wiring and the `SurfaceFluxes` contract.

This is what "modular like the rest of climt" buys us: process representations are first-class, swappable objects.

## Section 2.5 — Tests

- **Per-process unit tests**: each BEST-default process checked against hand-computed values from the notebook equations (soil properties for clay/sand, albedo, drag coefficients across stable/unstable Ri, exfiltration limit).
- **Substitution test**: swap in a trivial stub process and assert the orchestrator calls it and threads its outputs (proves the interface, not the physics).
- **Conservation**: energy and water closure over an integration with the full BEST stack.
- **Integration**: SecondBEST driven by a fixed atmospheric forcing reaches a sensible equilibrium surface temperature and soil moisture; cached-result regression once validated.

---

## Build order

1. BucketHydrology `num_layers` (small, self-contained, immediate value; fixes the `0.15` clamp bug).
2. SecondBEST scaffolding: orchestrator + interfaces + `SoilProperties`, `SurfaceAlbedo`, `SurfaceLayer` (all algebraic, easy to unit-test).
3. `SurfaceFluxes` (depends on 2).
4. `SubsurfaceTransport` (heaviest; reuses the `surface_ice.py` tridiagonal solver — coordinate with the ocean/ice doc's shared `_snow_ice_column` core so the implicit solver is written once).
5. Integration, conservation tests, notebook "Extending" section.

## Open questions / decisions deferred to implementation

- Exact default `moisture_diffusion_timescale` and `k_soil` values — pick from BEST/literature during implementation; expose as config regardless.
- Whether `SubsurfaceTransport` and the sea/land-ice `_snow_ice_column` solver share one tridiagonal core module in `_core/` (recommended) — settle when Part 4 of the ocean/ice doc is planned.

# Picket-Fence Radiation Scheme for climt

## Design Specification

**Date:** 2026-04-12
**Status:** Draft
**Components:** `PicketFenceLongwave`, `PicketFenceShortwave`

---

## 1. Motivation

climt currently has two radiation options:

- **GrayLongwaveRadiation** — single broadband, no gas-specific absorption, useful only for idealized experiments
- **RRTMG LW/SW** — full 16/14-band correlated-k with Fortran backend, accurate but heavy and not pure Python

There is no intermediate option: a scheme that captures real gas absorption (CO2, H2O, O3, etc.) with realistic spectral structure, runs in pure Python with numba acceleration, and is extensible to arbitrary atmospheres and absorbers.

The picket-fence radiation scheme fills this gap. It provides two optics modes:

1. **Parmentier mode** — bulk mixture opacities from Rosseland mean fits and dimensionless ratio coefficients (Parmentier & Guillot 2014, 2015; Lee et al. 2021). Fast, proven in 3D GCMs for exoplanet atmospheres. 5 bands total (2 LW + 3 SW).

2. **Correlated-k mode** — per-gas k-distribution tables generated from HITRAN line data via `linepyline`. Enables gas-specific forcing experiments (e.g., doubling CO2). Configurable band/g-point resolution.

Both modes share the same solver kernels, cloud hooks, and component interface.

### Design goals

- Pure Python + numba, with pure Python fallback (Pyodide-compatible at low resolution)
- Low-res correlated-k: ~400 iterations/sec in pure Python for single-column RCE
- High-res correlated-k: similar performance with numba acceleration
- Pedagogically transparent: students can read the solver and understand the physics
- Extensible: new gases via new tables, clouds via optical depth input, new planets via new data files
- Consistent with existing climt conventions (sympl components, state dict, naming)

---

## 2. Architecture

### 2.1 Three-layer design

```
Layer 3: climt Components (user-facing API)
  PicketFenceLongwave(TendencyComponent)
  PicketFenceShortwave(TendencyComponent)

Layer 2: Solver kernels (shared, numba-jitted)
  _lw_transport(tau, source, weights)
  _sw_two_stream(tau, ssa, g, zenith, albedo, solar_flux, weights)

Layer 1: Gas optics functions (two flavors)
  _parmentier_gas_optics(T, p, T_eff, kappa_R_coeffs, ratio_coeffs)
  _correlated_k_gas_optics(T, p, gas_amounts, table)
```

**Layer 1 (Gas Optics)** computes optical depths per band from atmospheric state. The two flavors differ in how they do this but produce the same output shape: optical depth arrays and source functions.

**Layer 2 (Solver Kernels)** consumes optical depths + source functions, performs radiative transfer sweeps, and returns per-band and broadband fluxes. These are pure numba functions with `prange` parallelization over columns. They are agnostic to the optics mode.

**Layer 3 (Components)** wraps everything in climt's TendencyComponent interface. Handles state dict extraction, reshaping, unit conversion, and mode dispatch.

### 2.2 Data flow

```
state dict
  → extract arrays (T, p, gas concentrations, cloud tau, surface properties)
  → reshape to 2D (nlev, ncol)
  → gas optics function (Parmentier or correlated-k)
      → tau[nband, ngpt, nlev, ncol], source[nband, ngpt, nlev, ncol]
  → add cloud optical depth
  → solver kernel
      → fluxes_per_band[nband, nlev+1, ncol]
      → fluxes_broadband[nlev+1, ncol]
  → heating rate from broadband flux divergence
  → reshape back to original dims
  → return (tendencies, diagnostics)
```

---

## 3. Component Interface

### 3.1 `PicketFenceLongwave(TendencyComponent)`

**Initialization:**

```python
# Parmentier mode (exoplanet work)
lw = PicketFenceLongwave(
    optics="parmentier",
    rosseland_mean_fit="freedman2014",   # or path to custom fit
    coefficients="solar_composition",     # or path to custom coefficient file
)

# Correlated-k mode (terrestrial climate)
lw = PicketFenceLongwave(
    optics="correlated_k",
    table="earth_low_res",    # or "earth_high_res", "mars", "venus", "trappist1e", or path
)
```

**Input properties (both modes):**

| Quantity | Dims | Units | Alias |
|---|---|---|---|
| `air_temperature` | mid_levels, * | degK | T |
| `air_pressure` | mid_levels, * | Pa | p |
| `air_pressure_on_interface_levels` | interface_levels, * | Pa | p_int |
| `surface_temperature` | * | degK | T_surf |
| `surface_longwave_emissivity` | num_longwave_bands, * | dimensionless | emissivity |
| `longwave_optical_thickness_due_to_cloud` | mid_levels, *, num_longwave_bands | dimensionless | tau_cloud |

**Additional inputs (Parmentier mode only):**

| Quantity | Dims | Units |
|---|---|---|
| `irradiation_temperature` | * | degK |
| `internal_temperature` | * | degK |

Note: these quantities are new and must be added to `default_values` in `climt/_core/initialization.py` with sensible defaults (e.g., `irradiation_temperature`: 0 K, `internal_temperature`: 0 K).

**Additional inputs (correlated-k mode only):**

| Quantity | Dims | Units |
|---|---|---|
| `specific_humidity` | mid_levels, * | kg/kg |
| `mole_fraction_of_carbon_dioxide_in_air` | mid_levels, * | mole/mole |
| (+ other gases as defined by the loaded table) | | |

**Tendency properties:**

| Quantity | Units |
|---|---|
| `air_temperature` | degK s^-1 |

**Diagnostic properties:**

| Quantity | Dims | Units |
|---|---|---|
| `upwelling_longwave_flux_in_air` | interface_levels, * | W m^-2 |
| `downwelling_longwave_flux_in_air` | interface_levels, * | W m^-2 |
| `upwelling_longwave_flux_in_air_per_band` | interface_levels, *, num_longwave_bands | W m^-2 |
| `downwelling_longwave_flux_in_air_per_band` | interface_levels, *, num_longwave_bands | W m^-2 |
| `longwave_heating_rate` | mid_levels, * | degK day^-1 |

### 3.2 `PicketFenceShortwave(TendencyComponent)`

**Initialization:**

```python
# Parmentier mode
sw = PicketFenceShortwave(
    optics="parmentier",
    stellar_spectrum="sun",           # or "hd209458", "trappist1", or path
    rosseland_mean_fit="freedman2014",
    coefficients="solar_composition",
)

# Correlated-k mode
sw = PicketFenceShortwave(
    optics="correlated_k",
    table="earth_low_res",
    stellar_spectrum="sun",
)
```

**Additional input properties (beyond the LW set, both modes):**

| Quantity | Dims | Units |
|---|---|---|
| `zenith_angle` | * | radians |
| `surface_albedo_for_direct_shortwave` | * | dimensionless |
| `flux_adjustment_for_earth_sun_distance` | * | dimensionless |

**Additional diagnostic properties (beyond the LW pattern):**

| Quantity | Dims | Units |
|---|---|---|
| `upwelling_shortwave_flux_in_air` | interface_levels, * | W m^-2 |
| `downwelling_shortwave_flux_in_air` | interface_levels, * | W m^-2 |
| `upwelling_shortwave_flux_in_air_per_band` | interface_levels, *, num_shortwave_bands | W m^-2 |
| `downwelling_shortwave_flux_in_air_per_band` | interface_levels, *, num_shortwave_bands | W m^-2 |
| `shortwave_heating_rate` | mid_levels, * | degK day^-1 |

---

## 4. Gas Optics: Parmentier Mode

### 4.1 Overview

Based on Parmentier & Guillot (2014, 2015) and as implemented in a GCM by Lee et al. (2021). The opacity structure is described by a small number of dimensionless parameters derived from the Rosseland mean opacity of the atmospheric mixture.

### 4.2 Thermal (LW) bands

Two bands with opacities kappa_1 and kappa_2, where kappa_1 > kappa_2. The Planck function is split between bands by the parameter beta:

- Band 1 (high opacity): fraction beta of the Planck function, opacity kappa_1
- Band 2 (low opacity): fraction (1 - beta) of the Planck function, opacity kappa_2

The opacities are derived from:

```
kappa_R(p, T)          — Rosseland mean opacity from pre-computed fit (e.g., Freedman et al. 2014)
gamma_1 = kappa_1 / kappa_R
gamma_2 = kappa_2 / kappa_R
gamma_P = kappa_P / kappa_R   — Planck-to-Rosseland mean ratio
R = kappa_1 / kappa_2 = gamma_1 / gamma_2
```

The parameters gamma_P, beta, and R are functions of T_eff, provided by fitting tables (Parmentier et al. 2015, Table 1). Given these, gamma_1 and gamma_2 are computed from the relations in Parmentier & Guillot (2014) Eqs. 87-96.

### 4.3 Visible (SW) bands

Three bands of equal spectral width (beta_vi = 1/3 each), with opacities kappa_v1 > kappa_v2 > kappa_v3. The ratios gamma_vi = kappa_vi / kappa_R are also functions of T_eff from the same fitting tables.

### 4.4 Procedure per column (following Lee et al. 2021, Section 2.4)

1. Calculate T_eff for the column from irradiation and internal temperatures
2. Look up gamma_v1, gamma_v2, gamma_v3, beta, gamma_P, R from coefficient tables
3. Compute kappa_R(p, T) in each layer from the Rosseland mean opacity fit
4. Compute the 3 visible band opacities: kappa_vi = gamma_vi * kappa_R
5. Compute the 2 thermal band opacities: kappa_1, kappa_2 from gamma_1, gamma_2, kappa_R
6. Build optical depth structure, pass to solver

### 4.5 Data files

- **Coefficient tables**: Parmentier et al. (2015) Table 1 — piecewise linear fits of log10(gamma_vi), beta, log10(gamma_P) as functions of log10(T_eff). Small (~1 KB).
- **Rosseland mean opacity fits**: Freedman et al. (2014) or Valencia et al. (2013) polynomial coefficients. Small (~1 KB).
- **Custom coefficient files**: same format, for users who calibrate for other atmospheric compositions (see student project brief).

---

## 5. Gas Optics: Correlated-k Mode

### 5.1 Overview

Uses pre-computed k-distribution tables derived from HITRAN line-by-line calculations. Each gas has absorption coefficients tabulated as a function of (temperature, pressure) for each spectral band and g-point. This enables gas-specific forcing experiments.

### 5.2 Table format (netCDF)

**Dimensions:**

| Name | Description |
|---|---|
| `band` | Number of spectral bands |
| `gpoint` | Maximum number of g-points per band |
| `temperature` | T grid points for interpolation |
| `pressure` | p grid points for interpolation (stored as log(p)) |
| `gas` | Number of absorbing gases |

**Variables:**

| Name | Shape | Description |
|---|---|---|
| `k_coefficients` | (gas, band, gpoint, temperature, pressure) | Absorption coefficients, m^2/kg |
| `gpoint_weights` | (band, gpoint) | Quadrature weights, sum to 1 per band |
| `planck_fraction` | (band, gpoint, temperature) | Fraction of Planck function per g-point (LW) |
| `band_wavenumber_limits` | (band, 2) | Lower/upper wavenumber bounds, cm^-1 |
| `temperature_grid` | (temperature) | K |
| `pressure_grid_log` | (pressure) | ln(Pa) |
| `solar_source_per_gpoint` | (band, gpoint) | W/m^2 per g-point (SW only) |
| `rayleigh_coefficient` | (band) | Rayleigh scattering cross-section per band (SW) |

**Attributes:**

| Name | Description |
|---|---|
| `gas_names` | Ordered list matching gas dimension, e.g., ["h2o", "co2"] |
| `overlap_method` | "additive" or "esft" |
| `resolution` | "low" or "high" (informational) |
| `source` | Provenance string |

### 5.3 Table interpolation

Absorption coefficients are looked up by bilinear interpolation in (log(p), T) space for each gas, band, and g-point. This is a standard approach shared with RRTMGP.

Water vapor is accepted as `specific_humidity` (kg/kg) in the state dict and converted to mole fraction internally, consistent with existing climt conventions.

### 5.4 Gas overlap

**Additive (low-res mode):** Optical depths from different gases are summed directly:

```
tau_total[band, gpt, layer] = sum_over_gases(k[gas, band, gpt] * column_amount[gas, layer])
```

This is equivalent to multiplying transmittances and assumes the gases have correlated spectral structure within each band. Simple and pedagogically transparent.

**ESFT (high-res mode):** The equivalent sum of exponentials with full terms. For N_gas gases each with G g-points, the solver evaluates G^N_gas combined g-points per band. For combined g-point (i, j, ...):

```
tau_combined = k_gas1[i] * amount1 + k_gas2[j] * amount2 + ...
weight_combined = w_gas1[i] * w_gas2[j] * ...
```

The combined weights sum to 1. This correctly handles spectrally uncorrelated gases at the cost of multiplicative g-point count. Practical for 2-3 major gases per band.

The overlap method is determined by the `overlap_method` attribute in the table file.

### 5.5 Example table sizes

| Configuration | Bands | g-points | Gases | T pts | p pts | Size |
|---|---|---|---|---|---|---|
| Low-res LW | 4 | 2 | 2 | 14 | 20 | ~36 KB |
| High-res LW | 16 | 8 | 4 | 14 | 20 | ~11 MB |

### 5.6 Pre-built tables shipped with climt

| Name | Atmosphere | Bands (LW+SW) | Gases |
|---|---|---|---|
| `earth_low_res` | Earth | 4+4 | H2O, CO2 |
| `earth_high_res` | Earth | 16+14 | H2O, CO2, O3, CH4 |
| `mars` | Mars (95% CO2) | TBD | CO2, H2O |
| `venus` | Venus (96% CO2) | TBD | CO2, H2O, SO2 |
| `trappist1e` | TRAPPIST-1e | TBD | H2O, CO2, N2 |

Table generation is performed externally using `linepyline` (GPL-3.0). climt (BSD) only loads the resulting netCDF files. See Section 10 for the table generation workflow.

---

## 6. Solver Kernels

### 6.1 Longwave transport

Absorption-only (no scattering) sweep, identical to the existing gray radiation kernel but executed per g-point. A diffusivity factor D=1.66 is applied to convert vertical optical depth to effective diffuse optical depth, consistent with the Eddington approximation (Elsasser 1942):

```python
@njit
def _lw_transport(T, T_surface, tau, planck_source, emissivity, weights, sigma):
    # D = 1.66  (diffusivity factor)
    # For each g-point:
    #   Upward sweep: surface emission + layer-by-layer absorption/emission
    #     trans = exp(-D * tau)
    #   Downward sweep: top-of-atmosphere (zero) + layer-by-layer
    #     trans = exp(-D * tau)
    # Weight by g-point weights
    # Sum over bands for broadband
    # Return per-band fluxes and broadband fluxes
```

The Planck source at each layer is computed from the `planck_fraction` table (correlated-k mode) or from beta weighting (Parmentier mode), multiplied by sigma * T^4.

### 6.2 Shortwave two-stream

Meador & Weaver (1980) two-stream approximation, following the same formulation used by RRTMGP:

```python
@njit
def _sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights):
    # For each g-point:
    #   Direct beam: exponential attenuation by total optical depth / cos(zenith)
    #   Diffuse: two-stream solution for layer reflectance and transmittance
    #   Adding method (Shonk & Hogan 2008) for inter-layer coupling
    # Weight by g-point weights
    # Sum over bands for broadband
```

Rayleigh scattering is added to the optical depth and single-scattering albedo before the solve. For the Parmentier mode, scattering is handled implicitly through the Bond albedo parameterization (Parmentier et al. 2015, Table 2).

### 6.3 Heating rate

Computed from broadband flux divergence, same as existing gray radiation:

```
heating_rate = (g / cp) * d(F_net) / dp
```

where F_net = F_up - F_down.

### 6.4 Numba acceleration pattern

All kernels use `@njit` with `prange` over columns, following the existing climt pattern:

```python
from climt._core.backend import prange

try:
    from numba import njit
except ImportError:
    def njit(x, **kwargs):
        return x
```

The outermost loop is over columns (parallelized), with inner loops over bands and g-points. All inputs are raw numpy arrays — no Python objects inside kernels.

### 6.5 Pedagogical introspection and entry points

The solver kernels implement the same Meador & Weaver / PIFM two-stream algorithm used in production codes (RRTMG, RRTMGP). Their pedagogical value comes from being pure Python: readable, modifiable, and debuggable interactively. To preserve and strengthen that value, the following features should be provided.

#### 6.5.1 Optional intermediate diagnostics

Both `sw_two_stream` and `lw_transport` should accept an optional `diagnostics_level` parameter (default 0). When set to a higher level, the solver returns additional per-layer quantities alongside the standard fluxes:

- **Level 1**: Per-layer diffuse reflectance and transmittance (`Rdif`, `Tdif`), direct beam transmittance (`Tnoscat`), and the direct beam flux profile at interfaces — enough to see what the adding method is building from.
- **Level 2**: Additionally, the direct-beam source terms (`Rdir`, `Tdir` per layer), the combined albedo array from the upward sweep, and the delta-scaled optical properties — the full internal state of the solver.

These are returned as a dictionary of arrays keyed by name. When `diagnostics_level=0`, no additional work is done and no dictionary is allocated. The component's `array_call` passes the flag through and, if non-zero, includes the extra arrays in the diagnostics output with names like `sw_layer_diffuse_reflectance`, `sw_direct_beam_profile`, etc.

This is the single most valuable pedagogical feature: it lets a student see *inside* the two-stream solve, which compiled Fortran codes cannot easily expose.

#### 6.5.2 Component-level pedagogical diagnostics

The primary pedagogical value of the picket-fence scheme is that it exposes the spectral physics that compiled codes hide: per-band optical depths, per-band transmittances, and per-band flux profiles. These let students see *why* the stratosphere cools (the window band is nearly transparent while the absorption band is opaque, so the stratosphere radiates more than it absorbs), and they can vary gas abundances and directly observe how optical depth changes — something rarely possible in full radiation codes.

Both `PicketFenceLongwave` and `PicketFenceShortwave` expose these quantities as standard diagnostics in their output dict, with no special flags required:

**LW additional diagnostics:**

| Quantity | Dims | Units | Description |
|---|---|---|---|
| `longwave_optical_depth_per_band` | mid_levels, \*, num_longwave_bands | dimensionless | Per-layer, per-band optical depth (tau) as computed by the optics function from gas abundances |
| `longwave_transmittance_per_band` | mid_levels, \*, num_longwave_bands | dimensionless | Per-layer, per-band transmittance exp(-D·tau), where D=1.66 is the diffusivity factor |
| `longwave_heating_rate_per_band` | mid_levels, \*, num_longwave_bands | degK day^-1 | Per-band contribution to the heating rate, so students can see which band drives cooling/warming at each level |

**SW additional diagnostics:**

| Quantity | Dims | Units | Description |
|---|---|---|---|
| `shortwave_optical_depth_per_band` | mid_levels, \*, num_shortwave_bands | dimensionless | Per-layer, per-band extinction optical depth |
| `shortwave_heating_rate_per_band` | mid_levels, \*, num_shortwave_bands | degK day^-1 | Per-band contribution to the SW heating rate |

These are always computed and returned — they add negligible cost since the tau arrays are already available inside `array_call` before being passed to the kernel. The per-band fluxes (`upwelling_longwave_flux_in_air_per_band`, etc.) are already in the standard output.

Together with the existing per-band fluxes and the `diagnostics_level` kernel introspection, this gives students a complete spectral decomposition of the radiation budget at every level — the key ingredient for understanding why non-grey radiation matters and why stratospheric cooling is a spectral phenomenon.

#### 6.5.3 Kernel function decomposition

The current decomposition into small, single-purpose functions (`_delta_scale`, `_sw_dif_and_source`, `_adding`, `sw_two_stream`) must be preserved as the code evolves. Each function represents one conceptual step of the algorithm and can be understood in isolation. Resist merging them for performance — the numba JIT already inlines them, and the readability cost is not worth the negligible gain.

#### 6.5.4 Worked Jupyter notebooks

Notebooks (in `examples/`) use the components directly to explore spectral radiation physics:

- **`examples/spectral_radiation_anatomy.ipynb`**: Run `PicketFenceLongwave` (correlated-k mode) on a realistic Earth temperature profile. Plot the per-band optical depth, per-band transmittance, and per-band flux profiles. Show that the window band is nearly transparent (transmittance close to 1 everywhere) while the absorption band is opaque (transmittance near zero in the lower troposphere). Use the per-band heating rates to demonstrate why the stratosphere cools: it emits in all bands but only absorbs in the opaque bands. Then double CO₂ and show how optical depth changes — making the "doubling CO₂" experiment concrete in optical depth space.
- **`examples/k_distribution_demo.ipynb`**: Start from a line-by-line absorption spectrum, construct a k-distribution, and show that a 2-gpoint quadrature reproduces the broadband transmission. Connect this to the correlated-k tables used by the component.
- **`examples/picket_fence_vs_rrtmg.ipynb`**: Compare picket-fence (Parmentier mode and correlated-k mode) against RRTMG for an Earth standard atmosphere, showing where the simplified scheme agrees and where it diverges.

---

## 7. Shortwave: Stellar and Orbital Parameters

### 7.1 Stellar spectrum

The component accepts a stellar spectrum identifier or file path at initialization. The spectrum (irradiance vs. wavenumber) is integrated over each band's wavenumber limits to produce per-band top-of-atmosphere flux. For Parmentier mode, the 3 visible bands have equal spectral width (1/3 of the stellar flux each, adjusted for the actual spectrum).

Pre-built spectra shipped with climt:

- `sun` — Lean & DeLand (2012) solar spectrum
- `trappist1` — TRAPPIST-1 spectrum
- Additional spectra can be provided as text files (wavenumber, irradiance columns)

### 7.2 Orbital parameters

The component does not compute orbital geometry internally. Instead, it accepts:

- `zenith_angle` — solar/stellar zenith angle (already in climt defaults)
- `flux_adjustment_for_earth_sun_distance` — scales the top-of-atmosphere flux (already in climt defaults)

Users compute these from orbital parameters (eccentricity, obliquity, longitude of perihelion, latitude, time of day) before calling the component. This keeps the radiation component simple and reusable across different planetary configurations. A utility function or separate diagnostic component can handle the orbital mechanics.

---

## 8. Atmospheric Properties: Multi-Planet Support

### 8.1 Motivation

Different planetary atmospheres require different physical constants: the molar mass of the bulk atmosphere, gas-specific molar masses for the absorbers that matter, gravitational acceleration, and so on. Hardcoding Earth values in component code prevents reuse for Mars, Titan, TRAPPIST-1 planets, or user-defined atmospheres.

climt solves this by routing all atmosphere-dependent constants through sympl's constants dictionary. Components call `get_constant("molar_mass_of_dry_air", "g/mol")` rather than hardcoding 28.97. A dedicated API loads the full set of relevant constants for a named atmosphere in one call, replacing sympl's defaults.

### 8.2 API

```python
import climt

# Load Earth constants (default at import)
climt.load_atmospheric_properties("earth")

# Switch to Mars
climt.load_atmospheric_properties("mars")

# Restore previous state
climt.reset_atmospheric_properties()

# Load from a custom TOML file
climt.load_atmospheric_properties("/path/to/my_atmosphere.toml")
```

`load_atmospheric_properties` takes a snapshot of the current sympl constants dictionary (deepcopy) before applying the new profile. `reset_atmospheric_properties` restores from that snapshot, reverting all changes made by the last load call.

### 8.3 Profile format (TOML)

Each profile is a TOML file containing every constant relevant to that atmosphere. Built-in profiles ship in `climt/_data/atmospheric_properties/`. The file name (without `.toml`) is the profile identifier.

```toml
# climt/_data/atmospheric_properties/earth.toml

[planetary]
gravitational_acceleration = { value = 9.80665, units = "m/s^2" }
planetary_radius            = { value = 6371000.0, units = "m" }
planetary_rotation_rate     = { value = 7.292e-05, units = "s^-1" }

[bulk_atmosphere]
molar_mass_of_dry_air                          = { value = 28.970, units = "g/mol" }
gas_constant_of_dry_air                        = { value = 287.0,  units = "J/kg/K" }
heat_capacity_of_dry_air_at_constant_pressure  = { value = 1004.64, units = "J/kg/K" }

[gas_species]
molar_mass_of_water_vapor      = { value = 18.015, units = "g/mol" }
molar_mass_of_carbon_dioxide   = { value = 44.010, units = "g/mol" }
molar_mass_of_ozone            = { value = 47.998, units = "g/mol" }
molar_mass_of_methane          = { value = 16.043, units = "g/mol" }
molar_mass_of_nitrous_oxide    = { value = 44.013, units = "g/mol" }
```

A Mars profile would set `gravitational_acceleration = 3.721 m/s^2`, `molar_mass_of_dry_air = 43.34 g/mol` (CO2-dominated), and omit water vapour or include it at lower priority. There is no requirement that profiles cover the same set of keys — each profile declares only the constants relevant to its atmosphere.

### 8.4 Fail-fast with helpful errors

If a component calls `get_constant` for a key not present in the loaded profile (e.g., requesting `molar_mass_of_water_vapor` on a dry Mars run), sympl raises a `KeyError`. climt wraps this in a more informative error:

```
ConstantNotFoundError: 'molar_mass_of_water_vapor' is not set in the current
atmospheric profile. To add it, either:
  1. Add it to your profile TOML under [gas_species]:
       molar_mass_of_water_vapor = { value = 18.015, units = "g/mol" }
  2. Set it directly: climt.set_constant("molar_mass_of_water_vapor", 18.015, "g/mol")
Current profile: mars (climt/_data/atmospheric_properties/mars.toml)
```

### 8.5 Built-in profiles

| Profile name | Description |
|---|---|
| `earth` | Earth standard atmosphere (loaded by default at `import climt`) |
| `mars` | Mars CO2-dominated atmosphere |
| `titan` | Titan N2/CH4 atmosphere |
| `hot_jupiter` | Generic H2/He-dominated hot Jupiter |
| `trappist1e` | TRAPPIST-1e estimated atmosphere |

### 8.6 File organization addition

```
climt/
  _data/
    atmospheric_properties/
      earth.toml
      mars.toml
      titan.toml
      hot_jupiter.toml
      trappist1e.toml
  _core/
    atmospheric_properties.py   # load_atmospheric_properties, reset_atmospheric_properties
```

`atmospheric_properties.py` exposes `load_atmospheric_properties` and `reset_atmospheric_properties` at the top-level `climt` namespace alongside the existing component exports.

---

## 9. Cloud Extensibility

Clouds contribute additional optical depth per layer per band. The input property `longwave_optical_thickness_due_to_cloud` (and the SW equivalent) defaults to zero and is added to the gas optical depth before the solver:

```
tau_total = tau_gas + tau_cloud
```

For the SW two-stream, clouds also contribute to single-scattering albedo and asymmetry parameter. These inputs (`shortwave_optical_thickness_due_to_cloud`, `single_scattering_albedo_due_to_cloud`, `cloud_asymmetry_parameter`) follow the existing RRTMG naming convention and default to zero/neutral values.

A future cloud optics component would fill these arrays. No changes to the radiation components are needed.

---

## 10. File Organization

```
climt/
  _components/
    picket_fence/
      __init__.py                  # exports PicketFenceLongwave, PicketFenceShortwave
      lw/
        component.py               # PicketFenceLongwave(TendencyComponent)
        kernels.py                 # _lw_transport and related numba functions
      sw/
        component.py               # PicketFenceShortwave(TendencyComponent)
        kernels.py                 # _sw_two_stream and related numba functions
      optics/
        parmentier.py              # _parmentier_gas_optics, coefficient loading
        correlated_k.py            # _correlated_k_gas_optics, table loading/interpolation
      common.py                    # shared utilities (heating rate, flux reduction, etc.)
  _data/
    picket_fence/
      parmentier/
        solar_composition.nc       # Parmentier et al. 2015 coefficients
        freedman2014.nc            # Rosseland mean opacity fit
      correlated_k/
        earth_low_res_lw.nc
        earth_low_res_sw.nc
        earth_high_res_lw.nc
        earth_high_res_sw.nc
        mars_lw.nc
        mars_sw.nc
        venus_lw.nc
        venus_sw.nc
        trappist1e_lw.nc
        trappist1e_sw.nc
      stellar_spectra/
        sun.nc
        trappist1.nc
```

---

## 11. Table Generation Workflow (External)

Table generation is not part of climt. It uses `linepyline` (GPL-3.0) and lives in a separate repository or in `linepyline`'s examples.

### Workflow:

1. User specifies: band boundaries, number of g-points per band, T/p grid, list of gases
2. For each gas, (T, p) point, and wavenumber in a band:
   - Call `linepyline.rtm()` to get high-resolution absorption cross-sections
3. For each band and (T, p) point:
   - Sort absorption coefficients by magnitude within the band
   - Compute the cumulative distribution function (CDF)
   - Fit the CDF with the specified number of g-points (Gaussian quadrature)
   - Record k-coefficients and weights
4. For LW: compute Planck fractions per g-point from band-integrated Planck function
5. For SW: integrate stellar spectrum over each band and distribute across g-points
6. Write netCDF file in the format specified in Section 5.2

### Coordination with linepyline

`linepyline` (maintained by Rodrigo Caballero) provides the line-by-line engine. The table generation scripts will be developed in coordination with the `linepyline` project. Since Rodrigo is a close collaborator, we can request updates or new features in `linepyline` as needed.

---

## 12. Documentation Requirements

The documentation for this component should be detailed and instructive, serving as a teaching resource for students learning radiative transfer. Topics to cover:

1. **Why non-grey?** The exponential-of-the-mean vs. mean-of-the-exponential problem. Why band-mean opacities fail and how k-distributions fix this.

2. **The k-distribution concept.** Sorting absorption coefficients by magnitude within a band. The CDF representation. Quadrature with g-points.

3. **Correlated-k and why it works.** The assumption that spectral ranking is preserved across (T, p). When it breaks down. Connection to the full line-by-line calculation.

4. **ESFT for gas overlap.** TODO: detailed explanation with worked example of the outer product of g-points, why the weights multiply, and the physical interpretation as random spectral overlap. Include a comparison showing when additive overlap is sufficient vs. when ESFT matters.

5. **The picket-fence model.** The Parmentier & Guillot formulation: Rosseland and Planck means, the ratio parameters (gamma_P, beta, R), what they mean physically, how they shape the temperature profile.

6. **Two-stream approximation.** The Meador & Weaver (1980) formulation, what the Eddington approximation is, when it's accurate (~2% for grey, ~4% for non-grey per Parmentier et al. 2015).

7. **Worked examples.** Jupyter notebooks as described in Section 6.5.4: spectral radiation anatomy with per-band diagnostics (optical depth, transmittance, heating rates), k-distribution construction from line-by-line, and picket-fence vs. RRTMG comparison. Additional scenarios: Earth RCE with CO2 doubling, Hot Jupiter T-p profile comparison between Parmentier and correlated-k modes, Mars greenhouse effect.

8. **Switching planetary atmospheres.** How to use `load_atmospheric_properties` and `reset_atmospheric_properties`. Cover: (a) the built-in profiles and what each sets, (b) how to write a custom TOML profile — document every valid key, the `{ value, units }` format, and which keys are required vs. optional for each component, (c) what happens when a required constant is missing (the error message and how to fix it), (d) worked examples: running the same RCE model for Earth then Mars in the same session, and defining a novel exoplanet atmosphere from scratch.

---

## 13. Validation Strategy

### 12.1 Against existing climt components

- Compare broadband fluxes and heating rates against RRTMG for Earth standard atmosphere profiles
- Compare gray limit (single band, uniform opacity) against GrayLongwaveRadiation

### 12.2 Against analytical solutions

- Compare Parmentier mode equilibrium T-p profiles against the analytical expressions in Parmentier & Guillot (2014) Eq. 76
- Verify grey limit recovers Milne's solution

### 12.3 Against published results

- Reproduce Lee et al. (2021) HD 209458b T-p profiles with Parmentier mode
- Compare Earth correlated-k results against RFMIP benchmark profiles (Pincus et al. 2016)

### 12.4 Internal consistency

- Verify that low-res and high-res correlated-k modes converge as resolution increases
- Verify that per-band fluxes sum to broadband fluxes
- Verify energy conservation: net flux at TOA = surface emission - OLR (for LW, no internal heat source)

---

## 14. Performance Targets

| Configuration | Mode | Target |
|---|---|---|
| Low-res correlated-k, pure Python, single column | No numba | ~400 iter/sec |
| High-res correlated-k, numba, single column | Numba | ~400 iter/sec |
| Parmentier mode, pure Python, single column | No numba | >1000 iter/sec |
| Low-res correlated-k, Pyodide | No numba | ~400 iter/sec |

---

## 15. References

- Parmentier, V. & Guillot, T. (2014). A non-grey analytical model for irradiated atmospheres. I. Derivation. A&A 562, A133.
- Parmentier, V., Guillot, T., Fortney, J. J., & Marley, M. S. (2015). A non-grey analytical model for irradiated atmospheres. II. Analytical vs. numerical solutions. A&A 574, A35.
- Lee, E. K. H., Parmentier, V., Hammond, M., et al. (2021). Simulating gas giant exoplanet atmospheres with EXO-FMS. MNRAS 506, 2695-2711.
- Pincus, R., Mlawer, E. J., & Delamere, J. S. (2019). Balancing accuracy, efficiency, and flexibility in radiation calculations for dynamical models. JAMES 11, 3074-3089.
- Freedman, R. S., Lustig-Yaeger, J., Fortney, J. J., et al. (2014). Gaseous mean opacities for giant planet and ultracool dwarf atmospheres. ApJS 214, 25.
- Meador, W. E. & Weaver, W. R. (1980). Two-stream approximations to radiative transfer in planetary atmospheres. JAS 37, 630-643.
- Koll, D. D. B. & Cronin, T. W. (2018). Earth's outgoing longwave radiation linear due to H2O greenhouse effect. PNAS 115, 10293-10298.

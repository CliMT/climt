# Notes on the Cork Radiation Scheme

Welcome to the review for the **Cork Radiation Scheme**! This feature introduces a dual-mode (Parmentier + correlated-k) intermediate radiation option for `climt`, written in pure Python and accelerated with `numba`.

## What This PR Does

Currently, `climt` forces users to choose between the ultra-simplified `GrayLongwaveRadiation` and the highly accurate but heavy, Fortran-compiled `RRTMG`. The cork scheme fills the gap by providing realistic spectral structure without compiled C/Fortran extensions.

### Two Modes of Operation
1. **Parmentier Mode:** A 5-band (2 LW + 3 SW) bulk mixture model based on Parmentier & Guillot (2014, 2015). It's incredibly fast and relies on dimensionless ratio coefficients derived from Rosseland mean opacities. This is primarily used for exoplanet (gas giant) atmospheres.
2. **Correlated-k Mode:** A classic $k$-distribution approach supporting per-gas absorption tables (generated via HITRAN/`linepyline`). This allows for dynamic gas concentrations (e.g. CO2 doubling experiments) and uses an additive overlap assumption in its low-res configuration.

Both modes output optical depth arrays and source functions that are consumed by shared solver kernels.

## Key Files to Review

### 1. The Core Kernels (`climt/_components/cork/*/kernels.py`)
These are the mathematical engines of the scheme, optimized with `@njit` and parallelized across columns using `prange`.
- **`lw/kernels.py`**: Implements the longwave transport sweep. It uses a standard diffusivity approximation ($D=1.66$) to convert vertical to effective diffuse optical depth.
- **`sw/kernels.py`**: Implements the shortwave two-stream solver. Currently, this is a simplified direct-beam-only solver (no diffuse scattering yet, just surface reflection).

### 2. The Components (`climt/_components/cork/*/component.py`)
These wrap the kernels in the `sympl` `TendencyComponent` interface.
- **`lw/component.py`**: Contains `CorkLongwaveRadiation`. Handles both Parmentier and correlated-k modes. Notice how it dynamically adjusts its `input_properties` based on the gases present in the correlated-k table!
- **`sw/component.py`**: Contains `CorkShortwaveRadiation`. Currently only supports Parmentier mode (3 visible bands). Notice the sign convention for `shortwave_heating_rate` (it's driven by `up - down` net flux absorption).

### 3. The Optics Layer (`climt/_components/cork/optics/`)
This translates atmospheric state into optical depths.
- **`parmentier.py`**: Evaluates piecewise linear fits for $\gamma$ coefficients and derives $\kappa_1$ and $\kappa_2$ from the Rosseland mean.
- **`correlated_k.py`**: Performs bilinear interpolation of $k$-coefficients in $(\log p, T)$ space and sums optical depths according to the additive overlap method.

### 4. Integration & Testing (`tests/test_cork_*.py`)
- We have comprehensive unit tests for the kernels, optics math, and individual components.
- **`test_cork_integration.py`** verifies that combining LW and SW produces physically sensible results (net cooling at TOA, net heating from SW absorption, etc.).
- You can find the synthetic 2-band table generator in the git history if you want to see how `test_2band_lw.npz` was built for the correlated-k tests!

## Notes for Reviewers

- **Initialization Defaults:** To support the Parmentier mode without crashing standard Earth-like tests, `irradiation_temperature` and `internal_temperature` were added to `climt/_core/initialization.py` with default values of $0$ K.
- **Placeholder Values:** You will notice a placeholder `kappa_R = 1e-4` in the Parmentier optics blocks. This is a stopgap until the full Freedman (2014) polynomial fit is integrated, but it allows the integration tests to run with Earth-like fluxes.
- **SW Diffuse:** The SW solver is currently "v1" (direct beam + reflection). Full Meador & Weaver two-stream is slated for a future follow-up, but the interface will remain identical.

Feel free to run `pytest tests/ -v` to see all 304 tests passing cleanly!

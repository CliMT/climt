# Student Project: Calibrating Cork Coefficients for Terrestrial Atmospheres

**Date:** 2026-04-12
**Prerequisites:** Familiarity with atmospheric radiation, Python, basic spectroscopy
**Tools:** `linepyline`, climt, HITRAN database
**Estimated scope:** 1 semester undergraduate research project or MSc thesis chapter

---

## 1. Background

The Parmentier & Guillot (2014, 2015) cork radiation model describes atmospheric opacity using a small number of dimensionless parameters:

- **gamma_P** = kappa_P / kappa_R — ratio of Planck mean to Rosseland mean opacity, quantifying non-greyness
- **beta** — relative spectral width of the high-opacity thermal channel (fraction of Planck function in the strong-absorption band)
- **R** = kappa_1 / kappa_2 — ratio of the two thermal opacities
- **gamma_v1, gamma_v2, gamma_v3** — ratios of visible band opacities to Rosseland mean

These parameters, combined with a Rosseland mean opacity fit kappa_R(T, p), fully determine the radiative transfer. The original calibration (Parmentier et al. 2015, Table 1; Freedman et al. 2014) is for **solar-composition gas giant atmospheres** (H2/He dominated). This project extends the calibration to terrestrial atmospheres: Earth, Mars, Venus, and potentially other rocky planet compositions.

### Why this matters

The Parmentier mode in climt's cork radiation scheme currently only works for gas giant atmospheres. Extending it to terrestrial compositions would enable fast, physically motivated radiation calculations for Earth-like planets without requiring k-distribution tables — useful for parameter sweeps, idealized GCM experiments, and building physical intuition about how atmospheric composition shapes the temperature profile.

---

## 2. Objectives

1. **Compute Rosseland and Planck mean opacities** for terrestrial atmospheric compositions (Earth, Mars, Venus) as functions of temperature and pressure, using line-by-line calculations from `linepyline`.

2. **Derive the cork parameters** (gamma_P, beta, R, gamma_v1, gamma_v2, gamma_v3) from the opacity cumulative distribution functions (CDFs) at representative atmospheric conditions.

3. **Fit the parameters** as functions of an appropriate control variable (T_eff for irradiated cases, or surface temperature for non-irradiated cases) and produce coefficient tables in the same format as Parmentier et al. (2015) Table 1.

4. **Fit the Rosseland mean opacity** kappa_R(T, p) as a polynomial or piecewise function for each atmospheric composition.

5. **Validate** the resulting analytical temperature profiles against numerical solutions from climt's correlated-k mode and/or RRTMG.

---

## 3. Method

### 3.1 Computing line-by-line opacities

Use `linepyline` to compute mass absorption coefficients kappa(nu, T, p) for each gas in the target atmosphere:

**Earth atmosphere:**
- Major absorbers: H2O, CO2, O3
- Minor absorbers: CH4, N2O
- Background gas: N2 + O2 (transparent, set molecular weight)
- T range: 180-320 K
- p range: 1-101325 Pa (1 Pa to 1 atm)
- Wavenumber range: 0-50000 cm^-1 (covers thermal + solar)

**Mars atmosphere:**
- Major absorbers: CO2 (95%), H2O (variable)
- T range: 140-300 K
- p range: 1-700 Pa
- Note: CO2 is not dilute — the pressure broadening is by CO2 itself

**Venus atmosphere:**
- Major absorbers: CO2 (96%), SO2, H2O
- T range: 200-740 K
- p range: 1-9.2e6 Pa (surface pressure ~92 atm)
- Note: extreme T/p range, CO2 continuum absorption important

### 3.2 Computing mean opacities

For each (T, p) point, compute:

**Rosseland mean:**

```
1/kappa_R = integral[ (1/kappa_nu) * (dB_nu/dT) dnu ] / integral[ (dB_nu/dT) dnu ]
```

**Planck mean:**

```
kappa_P = integral[ kappa_nu * B_nu dnu ] / integral[ B_nu dnu ]
```

where B_nu is the Planck function at temperature T.

The ratio gamma_P = kappa_P / kappa_R quantifies the non-greyness. For a grey atmosphere, gamma_P = 1. For real atmospheres, gamma_P >> 1 because the Planck mean is dominated by the strongest absorption lines while the Rosseland mean is dominated by the weakest.

### 3.3 Deriving cork parameters from the opacity CDF

Following Parmentier et al. (2015) Section 5.1.2:

1. At a representative atmospheric level (tau ~ 2/3, where the temperature is close to the effective temperature), compute the **cumulative distribution function** of the monochromatic opacities weighted by the local Planck function.

2. The CDF shows what fraction of the thermal emission encounters opacities below a given value. The Rosseland mean corresponds roughly to the 25th percentile, and the Planck mean to the 75th percentile (these are rough guidelines — see their Fig. 10).

3. **beta** is the fraction of the Planck-weighted spectral range where opacities exceed the Rosseland mean. Read this directly from the CDF.

4. **kappa_1** and **kappa_2** are representative opacities for the high and low opacity channels. One approach: kappa_1 is the Planck-weighted mean of opacities above kappa_R, and kappa_2 is the Planck-weighted mean of opacities below kappa_R. Then R = kappa_1/kappa_2.

5. For the **visible bands**, repeat the analysis but weight by the stellar spectrum instead of the Planck function, and split into 3 bands (highest, middle, lowest opacity thirds).

### 3.4 Fitting as functions of temperature

The parameters will vary with atmospheric conditions. For each composition:

1. Compute the parameters across a grid of temperatures and pressures
2. Identify the appropriate control variable:
   - For irradiated planets: T_eff (as in Parmentier et al. 2015)
   - For non-irradiated or weakly irradiated planets: may need to parameterize differently (e.g., as a function of surface temperature or pressure)
3. Fit piecewise linear or polynomial functions, following the format of Parmentier et al. (2015) Table 1
4. Pay attention to regime transitions (e.g., for Earth: the transition from water-vapor-dominated to CO2-dominated opacity in the upper troposphere/stratosphere)

### 3.5 Fitting Rosseland mean opacity

Fit kappa_R(T, p) as a function of temperature and log-pressure. For gas giants, Freedman et al. (2014) use a polynomial in (T, log p). For terrestrial atmospheres, the functional form may need to differ because:

- Water vapor concentration varies enormously (unlike gas giants where all species are well-mixed)
- The dominant absorber changes with altitude (H2O in troposphere, CO2 + O3 in stratosphere)

One approach: fit kappa_R as a function of (T, p, q) where q is specific humidity, treating the variable gas separately from the well-mixed ones.

---

## 4. Validation

### 4.1 Against correlated-k calculations

Run climt's cork radiation in correlated-k mode for standard atmospheric profiles (e.g., US Standard Atmosphere for Earth, Mars Climate Database profiles). Compare the equilibrium temperature profiles against those obtained using the newly calibrated Parmentier coefficients.

Target accuracy: temperature profiles within 10% of the correlated-k solution at pressures above the tropopause, following the accuracy achieved by Parmentier et al. (2015) for gas giants.

### 4.2 Against RRTMG

For Earth, compare broadband fluxes (OLR, surface downwelling LW) against RRTMG for a set of standard profiles. Document the bias and discuss which spectral features are captured vs. lost.

### 4.3 Sensitivity tests

- Does the calibrated model correctly capture the sign and approximate magnitude of CO2 doubling forcing for Earth?
- Does it reproduce the strong greenhouse effect of Venus?
- Does it capture the thin-atmosphere limit for Mars?

These tests will reveal which climate phenomena the 5-band cork model can and cannot represent for terrestrial atmospheres.

---

## 5. Deliverables

1. **Coefficient files** for Earth, Mars, and Venus in the format expected by climt's `CorkLongwaveRadiation` and `CorkShortwaveRadiation` components (netCDF).

2. **Rosseland mean opacity fits** for each atmosphere.

3. **A Jupyter notebook** documenting the full calibration procedure: line-by-line calculations, CDF analysis, parameter extraction, fitting, and validation plots. This notebook serves as both documentation and a reproducible recipe for calibrating new atmospheric compositions.

4. **A short report** (or paper draft) discussing:
   - How the terrestrial opacity structure differs from gas giants
   - Which parameters are most sensitive to composition
   - Where the 5-band cork approximation breaks down for terrestrial atmospheres
   - Recommendations for whether more thermal or visible bands are needed

---

## 6. Extensions (if time permits)

- **Variable gas handling**: For Earth, water vapor varies by orders of magnitude. Can the Parmentier framework be extended to handle a variable absorber, or does this fundamentally require per-gas tables (i.e., the correlated-k approach)?

- **Additional compositions**: Titan (N2 + CH4), early Earth (high CO2, no O2), Archean haze atmospheres.

- **Comparison with analytical solutions**: Use the Parmentier & Guillot (2014) Eq. 76 analytical T(tau) profile with the newly derived coefficients and compare against the numerical solution. This tests whether the Eddington approximation errors are similar for terrestrial and gas giant atmospheres.

- **Optimizing band structure**: Is the fixed 3 SW + 2 LW band structure optimal for terrestrial atmospheres, or would a different split (e.g., 2 SW + 3 LW, or 4 SW + 2 LW) give better results? This connects to the broader question of how many bands are "enough" for a given accuracy target.

---

## 7. Key References

- Parmentier, V. & Guillot, T. (2014). A non-grey analytical model for irradiated atmospheres. I. Derivation. A&A 562, A133.
- Parmentier, V., Guillot, T., Fortney, J. J., & Marley, M. S. (2015). A non-grey analytical model for irradiated atmospheres. II. Analytical vs. numerical solutions. A&A 574, A35.
- Lee, E. K. H. et al. (2021). Simulating gas giant exoplanet atmospheres with EXO-FMS. MNRAS 506, 2695-2711.
- Freedman, R. S. et al. (2014). Gaseous mean opacities for giant planet and ultracool dwarf atmospheres. ApJS 214, 25.
- Pierrehumbert, R. T. (2010). Principles of Planetary Climate. Cambridge University Press. (Chapters 4-5 on grey and non-grey radiation)
- Koll, D. D. B. & Cronin, T. W. (2018). Earth's outgoing longwave radiation linear due to H2O greenhouse effect. PNAS 115, 10293-10298.

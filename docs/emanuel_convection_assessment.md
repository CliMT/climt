# Assessment: Emanuel Convection Feature Parity

This document assesses the parity between the pure Python implementation of the Emanuel Convection scheme (`climt/_components/emanuel/pure_python.py`) and the reference Fortran 4.3c implementation (`climt/_lib/emanuel/convect43c.f90`), as described in Emanuel and Živković-Rothman (1999).

## Feature Parity Matrix

| Feature | Reference Fortran (`convect43c.f90`) | Pure Python (`pure_python.py`) | Status |
| :--- | :--- | :--- | :--- |
| **Relaxation Mechanism** | Relaxed CBMF based on LCL buoyancy (Eq. 2) | Implemented identically | **Parity** |
| **Parcel Lifting** | Iterative Newton-Raphson in `TLIFT` | Implemented in `_tlift` | **Parity** |
| **Buoyancy Sorting** | Full ensemble of mixtures (`SIJ`, `MENT`) | Matrix implemented; simplified weighting | **Partial** |
| **Momentum Transport** | Convective momentum transport (`FU`, `FV`) | Implemented (`UP`, `VP` logic needed) | **Partial** |
| **Passive Tracers** | Transport for `NTRA` tracers | Not implemented | **Missing** |
| **Dry Adjustment** | Dry adiabatic adjustment if `IPBL != 0` | Not implemented | **Missing** |
| **Surface Diagnostics** | Outputs `wd`, `tprime`, `qprime` (Eq. 4-6) | Implemented | **Parity** |
| **LNB / CAPE** | Finds Neutral Buoyancy Level and CAPE | Implemented with matrix search | **Parity** |
| **Downdraft Physics** | Hydrostatic `MP` with surface friction | Implemented; tuning likely needed | **Partial** |

## Implementation Status (Feb 2026 Update)

A parity test run on Feb 12, 2026, using an isothermal moist profile confirms that the Python implementation correctly **triggers** convection and calculates key thermodynamic levels identically to Fortran.

### Verified Parity
The following core variables match between Python and Fortran (within numerical precision):
*   **Parcel Origin Level (NK):** Matches (Fortran 5 vs Python index 4).
*   **LCL Pressure (PLCL):** Matches exactly (~914.4 hPa).
*   **Cloud Base Level (ICB):** Matches (Fortran 7 vs Python index 6).
*   **Cloud Base Mass Flux (CBMF):** Matches closely (~0.0075).
*   **Buoyancy/Trigger (DTMA):** Matches closely (~0.75).

### Critical Discrepancies

#### 1. Low-Level Thermodynamic Tendencies
**Issue:** In the lowest model levels (below cloud base), the Python implementation produces significant **heating (~16 K/day)**, whereas the Fortran implementation produces **zero tendency**.
*   **Hypothesis:** This error is likely in the vertical integration of mass fluxes (`AMP1`, `AD`) or the downdraft mass flux (`MP`) calculation in the tendency loop. The Fortran code has specific conditions (e.g., `IF(I.GE.NK)`) for `AMP1` accumulation that might be subtly different in the vectorized Python version.
*   **Impact:** This results in spurious warming in the sub-cloud layer, destabilizing the column artificially.

#### 2. Cloud Layer Tendency Profile
**Issue:**
*   **Fortran:** Shows strong cooling (~-15 K/day) at the top of the cloud (Level 7) and warming/drying at Level 4.
*   **Python:** Shows warming (~22 K/day) at Level 6/7 and weak cooling at Level 7.
*   **Cause:** The simplified "equal weighting" of mixtures in the `MENT`/`SIJ` matrix (vs the probability-mass weighting in Fortran) leads to a different vertical distribution of detrainment.

## Detailed Component Analysis

### Buoyancy Sorting (Mixing)
The Python code now implements the `SIJ` (mixing fraction) and `MENT` (entrainment mass flux) matrices.
*   **Gap:** The re-normalization step (Fortran lines 240-300) which weights mixtures based on their probability of occurrence (`SMIN`, `SMAX` logic) is simplified in Python. This affects *where* mass is detrained, leading to the observed profile differences.

### Momentum Transport
The Python code includes `fu` and `fv` arrays and updates them based on mass flux divergence.
*   **Gap:** The specific transport of momentum within the updraft/downdraft (`UP`, `VP` arrays in Fortran) is partially implemented but needs verification against the Fortran `UENT`/`VENT` logic.

### Conservation
*   **Enthalpy:** Python implements a column enthalpy correction.
*   **Strict Conservation:** Fortran includes a final adjustment loop (lines 640-700) to force exact enthalpy, momentum, and tracer conservation. This is partially implemented in Python but may lack the tracer component.

## Next Steps
1.  **Fix Low-Level Tendencies:** Debug the `FT` (temperature tendency) loop for levels `i < icb`. Specifically, check the `AMP1` and `AD` summation limits and the `MP` (downdraft) contribution terms.
2.  **Refine Mixing Weights:** Implement the probability-mass weighting for `MENT` normalization to match Fortran's detrainment profile.
3.  **Dry Adjustment:** Add the dry adiabatic adjustment block (optional for this specific test case but needed for full completeness).

"""Reproduce the HD 209458b purely-radiative equilibrium T-p profile.

Three sub-checks (A/B/C) run in one test after equilibrating with AdamsBashforth:
  A. OLR energy-balance: LW OLR ≈ σ*(T_int⁴+T_irr⁴) within 10%
  B. T-profile regression: converged profile matches hd209458b_model_equilibrium.csv within 2%
  C. Deep-atmosphere P&G analytic: at p>0.3 bar, T within 15% of P&G (2014)

Integration uses:
  - AdamsBashforth([lw, sw, slab]) — radiation + prognostic surface temperature
  - SimplePhysics(boundary_layer=False, use_external_surface_specific_humidity=True)
    called outside the AB stepper to provide sensible-heat coupling; latent heat
    is zeroed because surface_specific_humidity defaults to 0 for a dry atmosphere.

References:
  Torres et al. (2008), ApJ 677, 1324         -- planetary parameters
  Parmentier & Guillot (2014), A&A 562, A133  -- opacity / Teff parameterisation
  Guillot (2010), A&A 520, A27                -- analytic T(tau) formula
  Freedman et al. (2014), ApJS 214, 25        -- Rosseland mean opacity fit
"""

import os

import numpy as np
import pandas as pd
import pytest
from sympl import AdamsBashforth, set_backend

# ── Physical constants / reference values ────────────────────────────────────
STEFAN_BOLTZMANN = 5.670374419e-8  # W m⁻² K⁻⁴
T_IRR = 1450.0  # K  (Torres et al. 2008)
T_INT = 500.0  # K
OLR_TARGET = STEFAN_BOLTZMANN * (T_IRR**4 + T_INT**4)  # ≈ 254 203 W/m²

# ── File paths ────────────────────────────────────────────────────────────────
_FIGS = os.path.join(
    os.path.dirname(__file__), "..", "docs", "radiative_transfer", "figures"
)
REF_PATH = os.path.join(_FIGS, "hd209458b_reference.csv")
EQUIL_PATH = os.path.join(_FIGS, "hd209458b_model_equilibrium.csv")

# ── Integration parameters ────────────────────────────────────────────────────
N_STEPS = 150_000  # 60s × 1M ≈ 1.9 years of simulation time


@pytest.mark.slow
def test_hd209458b_equilibrium_profile():
    """HD 209458b radiative equilibrium: OLR balance, regression, and deep P&G check.

    Long RCE integration (~150k steps) — excluded from routine CI via the ``slow``
    marker; run on demand with ``pytest -m slow``.
    """
    from climt import (
        DryConvectiveAdjustment,
        SimplePhysics,
        SlabSurface,
        UnytBackend,
        UnytTimeDelta,
        get_default_state,
        load_atmospheric_properties,
        reset_atmospheric_properties,
    )

    set_backend(UnytBackend())
    from climt._components.cork import (
        CorkLongwaveRadiation,
        CorkShortwaveRadiation,
    )

    try:
        from tqdm import tqdm

        _iter = tqdm(range(N_STEPS), desc="HD 209458b equilibration", unit="step")
    except ImportError:
        _iter = range(N_STEPS)

    ref = pd.read_csv(REF_PATH, comment="#")

    load_atmospheric_properties("hot_jupiter")
    try:
        lw = CorkLongwaveRadiation(
            optics="parmentier",
            rosseland_mean_fit="freedman2014",
            coefficients="solar_composition",
        )
        sw = CorkShortwaveRadiation(
            optics="parmentier",
            stellar_spectrum="sun",
            rosseland_mean_fit="freedman2014",
            coefficients="solar_composition",
            bond_albedo_feedback=True,
        )
        slab = SlabSurface()
        dry_convection = DryConvectiveAdjustment()
        # boundary_layer=False: no BL mixing — pure radiative equilibrium.
        # use_external_surface_specific_humidity=True: reads q_surf from state,
        # which defaults to 0 for a dry hot-Jupiter atmosphere → LHF = 0.
        simple_physics = SimplePhysics(
            boundary_layer=False,
            large_scale_condensation=False,
            use_external_surface_specific_humidity=True,
        )

        state = get_default_state([lw, sw, slab, simple_physics, dry_convection])
        state["irradiation_temperature"].data[:] = T_IRR
        state["internal_temperature"].data[:] = T_INT
        state["zenith_angle"].data[:] = 0.0
        state["ocean_mixed_layer_thickness"].data[:] = 5.0  # thin slab
        state["area_type"].data[:] = "sea"
        # Constant 5 m/s surface forcing — Reed–Jablonowski SHF scales with |U|,
        # so without this the surface and atmosphere are radiatively decoupled.
        # Reset every step (see loop below) since SimplePhysics evolves the wind.
        U_FORCING = 5.0
        state["eastward_wind"].data[:] = U_FORCING
        state["northward_wind"].data[:] = 0.0

        # Warm-start: interpolate P&G analytic profile onto model pressure grid.
        # Level 0 is the bottom level (p≈1 bar) in the hot_jupiter default state.
        p_bar = np.asarray(state["air_pressure"].data)[:, 0, 0] / 1e5
        T_init = np.interp(
            np.log10(p_bar),
            np.log10(ref["pressure_bar"].values),
            ref["temperature_K"].values,
        )
        state["air_temperature"].data[:, :, :] = T_init[:, np.newaxis, np.newaxis]
        state["surface_temperature"].data[:] = T_init[0]

        TIMESTEP = UnytTimeDelta(seconds=120)
        stepper = AdamsBashforth([lw, sw, slab])

        print(
            f"\n{'step':>6}  {'T_surf':>10}  {'net_rad':>12}  {'SHF':>10}  {'OLR':>10}"
        )
        for step in _iter:
            diagnostics, state = stepper(state, TIMESTEP)
            state.update(diagnostics)
            diag_sp, new_state_sp = simple_physics(state, TIMESTEP)
            state.update(diag_sp)
            state.update(new_state_sp)
            diag_dc, new_state_dc = dry_convection(state, TIMESTEP)
            state.update(diag_dc)
            state.update(new_state_dc)
            # Reset wind to keep SHF forcing constant — SimplePhysics applies
            # surface drag that would otherwise spin the column down to zero.
            state["eastward_wind"].data[:] = U_FORCING
            state["northward_wind"].data[:] = 0.0

            if step % 5000 == 0:
                T_s = float(np.asarray(state["surface_temperature"].data).flat[0])
                lw_up = float(
                    np.asarray(diagnostics["upwelling_longwave_flux_in_air"].data).flat[
                        0
                    ]
                )
                lw_dn = float(
                    np.asarray(
                        diagnostics["downwelling_longwave_flux_in_air"].data
                    ).flat[0]
                )
                sw_dn = float(
                    np.asarray(
                        diagnostics["downwelling_shortwave_flux_in_air"].data
                    ).flat[0]
                )
                sw_up = float(
                    np.asarray(
                        diagnostics["upwelling_shortwave_flux_in_air"].data
                    ).flat[0]
                )
                shf = float(
                    np.asarray(state["surface_upward_sensible_heat_flux"].data).flat[0]
                )
                net_rad = lw_dn + sw_dn - lw_up - sw_up
                olr_now = float(
                    np.asarray(diagnostics["upwelling_longwave_flux_in_air"].data).flat[
                        -1
                    ]
                )
                print(
                    f"{step:>6}  {T_s:>10.2f}  {net_rad:>12.2f}  {shf:>10.2f}  {olr_now:>10.2f}"
                )

        # Final LW call for OLR diagnostic
        _, diag_lw = lw(state)
        # Upwelling flux shape is (n_lev+1, ncol, ncol); flat index -1 is TOA.
        olr = float(np.asarray(diag_lw["upwelling_longwave_flux_in_air"].data).flat[-1])
        T_model = np.asarray(state["air_temperature"].data)[:, 0, 0].copy()

    finally:
        reset_atmospheric_properties()

    # ── Check A: OLR energy balance ───────────────────────────────────────────
    olr_err = abs(olr - OLR_TARGET) / OLR_TARGET
    assert olr_err < 0.10, (
        f"OLR {olr:.1f} W/m² deviates {olr_err:.1%} from target "
        f"{OLR_TARGET:.1f} W/m² (tolerance 10%)"
    )

    # ── Check B: T-profile regression ─────────────────────────────────────────
    if not os.path.exists(EQUIL_PATH):
        df = pd.DataFrame({"pressure_bar": p_bar, "temperature_K": T_model})
        df.to_csv(EQUIL_PATH, index=False)
        pytest.skip(
            "Generated new equilibrium reference — re-run tests to validate regression"
        )

    ref_eq = pd.read_csv(EQUIL_PATH)
    T_ref_eq = np.interp(
        np.log10(p_bar),
        np.log10(ref_eq["pressure_bar"].values),
        ref_eq["temperature_K"].values,
    )
    reg_err = np.abs(T_model - T_ref_eq) / T_ref_eq
    worst_reg = int(np.nanargmax(reg_err))
    assert np.nanmax(reg_err) < 0.02, (
        f"Regression error {np.nanmax(reg_err):.2%} > 2% at "
        f"p={p_bar[worst_reg]:.4f} bar  "
        f"T_model={T_model[worst_reg]:.1f} K  T_ref={T_ref_eq[worst_reg]:.1f} K"
    )

    # ── Check C: deep-atmosphere P&G analytic (p > 0.3 bar only) ─────────────
    T_ref_pg = np.interp(
        np.log10(p_bar),
        np.log10(ref["pressure_bar"].values),
        ref["temperature_K"].values,
    )
    mask = p_bar > 0.3
    T_m = T_model[mask]
    T_r = T_ref_pg[mask]
    p_m = p_bar[mask]

    pg_err = np.abs(T_m - T_r) / T_r
    worst_pg = int(np.nanargmax(pg_err))
    assert np.nanmax(pg_err) < 0.15, (
        f"Deep P&G deviation {np.nanmax(pg_err):.2%} > 15% at "
        f"p={p_m[worst_pg]:.4f} bar  "
        f"T_model={T_m[worst_pg]:.1f} K  T_P&G={T_r[worst_pg]:.1f} K"
    )


if __name__ == "__main__":
    test_hd209458b_equilibrium_profile()

# Doc seed — explaining correlated-k (with a real failure case)

**Status:** note / raw material, not a finished doc. When the cork
correlated-k mode gets a proper write-up (a docs page alongside
`docs/PLUG_AND_PLAY_ARCHITECTURE.md`, or a user guide for `CorkLongwaveRadiation`),
fold this in. The value here is a *worked example, grounded in a real debugging
session*, of both how correlated-k works and how it fails — far more memorable
than an abstract description.

## Why this is worth writing up

While diagnosing a radiative-convective-equilibrium warm bias (CORK ran ~13 K too
warm vs RRTMG), we traced ~38 W/m² of a 51 W/m² OLR over-trapping to a single
cause: the 800–1800 cm⁻¹ "window" band of the Earth LW table **lumps two
spectrally distinct populations** — the transparent 8–12 µm window and the
opaque H₂O ν₂ band — into one 8-g-point k-distribution. That is the classic
correlated-k failure mode, and we have the `k(g)` curves to show it. (Full
story: experiment #16 in
`docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md`.)

## Narrative beats to use

1. **The problem.** Band-mean transmission `T(u) = ⟨exp(−k(ν)u)⟩` can't use a
   mean k (average of exp ≠ exp of average). Line-by-line is exact but costly.
2. **The k-distribution trick.** The integral only depends on the *distribution*
   of k values in the band, not where in ν they sit. Sort k(ν) weakest→strongest
   → a smooth monotonic `k(g)`, g ∈ [0,1] = cumulative band fraction. Then a few
   quadrature points (**g-points**) integrate it: low-g = transparent gaps
   (escape to space), high-g = saturated line cores.
3. **The "correlated" assumption.** Exact for one homogeneous layer; stacking
   layers requires that the rank ordering of k(ν) is preserved with height (a
   given g-point = the same wavenumbers at all levels). Holds when a band is
   dominated by one absorber scaling uniformly with p, T.
4. **How it fails — band too wide.** If a band straddles two regions with
   different opacity *and* different emission height, one shared Planck weight +
   shared g-points can't represent both → the transparent part can't open, gets
   over-absorbed. The rank assumption also decorrelates (the two regions scale
   differently with p, T).
5. **The fix.** Split the band at the opacity boundary so each piece is a single
   smooth population with its own g-points and Planck weight. (Same surgery the
   plan applied to the 500–800 cm⁻¹ CO₂ band, where a strong 667 cm⁻¹ core
   diluted the weak wings.)

## The figure

`scripts/experiments/plot_kg_curves.py` → `debug_data/kg_curves_window_band.png`
(regenerate with the climt env; move/copy into docs assets when the page is
written). Panel A — the 800–1800 band dry vs moist: dry is a low, smooth,
transparent curve; **moist develops a low-g "shelf" (window) + a high-g "cliff"
(H₂O ν₂ lines) with a ~30× step near g≈0.5 — the fingerprint of two populations
in one band.** Panel B — the lumped band's two-piece curve vs the smooth,
single-population curves of spectrally coherent bands (CO₂ core, H₂O rotation).

The dry-vs-moist contrast is the punchline: the error is *invisible when dry*
(both halves transparent) and switches on with H₂O (the cliff appears while the
shelf stays low) — which is exactly the dry/moist signature measured in the
forward-flux comparison (dry deficit −13 W/m², moist −51 W/m²).

## Caveat to carry into the doc

Keep the *other* ~13 W/m² (dry, CO₂-only) separate: that residual is a
solver-level issue (isothermal-layer Planck source vs RRTMG's linear-in-τ), not
a correlated-k band issue. Don't conflate the two when teaching.

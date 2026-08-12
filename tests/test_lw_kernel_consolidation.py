import numpy as np


def _reference_lw_transport(tau, planck_source, surface_source, emissivity, weights):
    """Frozen copy of the ORIGINAL python-loop transport (pre-consolidation),
    used as the bit-for-bit oracle."""
    D = 1.66
    nband, ngpt, nlev, ncol = tau.shape
    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))
    for b in range(nband):
        for g in range(ngpt):
            up = np.zeros((nlev + 1, ncol))
            dn = np.zeros((nlev + 1, ncol))
            for i in range(ncol):
                up[0, i] = emissivity[b, i] * surface_source[b, g, i]
                for k in range(nlev):
                    t = np.exp(-D * tau[b, g, k, i])
                    up[k + 1, i] = up[k, i] * t + planck_source[b, g, k, i] * (1.0 - t)
                dn[nlev, i] = 0.0
                for k in range(nlev - 1, -1, -1):
                    t = np.exp(-D * tau[b, g, k, i])
                    dn[k, i] = dn[k + 1, i] * t + planck_source[b, g, k, i] * (1.0 - t)
            w = weights[b, g]
            for k in range(nlev + 1):
                for i in range(ncol):
                    up_band[b, k, i] += w * up[k, i]
                    down_band[b, k, i] += w * dn[k, i]
    up_broad = up_band.sum(axis=0)
    down_broad = down_band.sum(axis=0)
    return up_band, down_band, up_broad, down_broad


def _make_inputs(seed=0):
    rng = np.random.RandomState(seed)
    nband, ngpt, nlev, ncol = 14, 8, 20, 3
    tau = rng.uniform(0.0, 2.0, size=(nband, ngpt, nlev, ncol))
    planck = rng.uniform(50.0, 400.0, size=(nband, ngpt, nlev, ncol))
    surf = rng.uniform(50.0, 400.0, size=(nband, ngpt, ncol))
    emis = rng.uniform(0.8, 1.0, size=(nband, ncol))
    weights = rng.uniform(0.1, 0.5, size=(nband, ngpt))
    return tau, planck, surf, emis, weights


def test_consolidated_kernel_matches_reference_bitwise():
    """The refactored kernel must reproduce current fluxes exactly."""
    from climt._components.cork.lw.kernels import lw_transport

    tau, planck, surf, emis, weights = _make_inputs()
    T = np.zeros((20, 3))
    T_surf = np.zeros(3)
    sigma = 5.67e-8

    ub_r, db_r, ubr_r, dbr_r = _reference_lw_transport(tau, planck, surf, emis, weights)
    ub, db, ubr, dbr = lw_transport(
        T, T_surf, tau, planck, surf, emis, weights, sigma, diagnostics_level=0
    )

    np.testing.assert_allclose(ub, ub_r, rtol=1e-15)
    np.testing.assert_allclose(db, db_r, rtol=1e-15)
    np.testing.assert_allclose(ubr, ubr_r, rtol=1e-15)
    np.testing.assert_allclose(dbr, dbr_r, rtol=1e-15)


def test_diagnostics_level_one_still_returns_diag():
    """diagnostics_level=1 must still return diag dict with correct shapes."""
    from climt._components.cork.lw.kernels import lw_transport

    tau, planck, surf, emis, weights = _make_inputs(seed=1)
    T = np.zeros((20, 3))
    T_surf = np.zeros(3)
    sigma = 5.67e-8

    out = lw_transport(
        T, T_surf, tau, planck, surf, emis, weights, sigma, diagnostics_level=1
    )
    assert len(out) == 5
    diag = out[4]
    assert diag["transmittance"].shape == (14, 8, 20, 3)
    assert diag["up_per_gpoint"].shape == (14, 8, 21, 3)
    assert diag["down_per_gpoint"].shape == (14, 8, 21, 3)
    # Diagnostic weighted up-flux sums over g to up_band.
    np.testing.assert_allclose(diag["up_per_gpoint"].sum(axis=1), out[0], rtol=1e-12)

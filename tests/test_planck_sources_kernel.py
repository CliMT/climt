import numpy as np


def _reference(planck_frac, T_grid, T, T_surf, sigma, nband, ngpt,
               is_esft, ngpt_orig, nband_orig):
    """Frozen copy of the original pure-Python Planck-source loop (oracle)."""
    nT = len(T_grid)
    nlev, ncol = T.shape
    planck_src = np.zeros((nband, ngpt, nlev, ncol))
    surf_src = np.zeros((nband, ngpt, ncol))
    for icol in range(ncol):
        T_s = T_surf[icol]
        iTs = np.searchsorted(T_grid, T_s) - 1
        iTs = max(0, min(iTs, nT - 2))
        fTs = (T_s - T_grid[iTs]) / (T_grid[iTs + 1] - T_grid[iTs])
        fTs = max(0.0, min(1.0, fTs))
        surf_planck = sigma * T_s ** 4
        for ib in range(nband):
            for igp in range(ngpt):
                g_idx_orig = igp % ngpt_orig if (is_esft) else igp
                ib_orig = min(ib, nband_orig - 1)
                frac_s = (planck_frac[ib_orig, g_idx_orig, iTs] * (1 - fTs)
                          + planck_frac[ib_orig, g_idx_orig, iTs + 1] * fTs)
                surf_src[ib, igp, icol] = frac_s * surf_planck
        for k in range(nlev):
            T_l = T[k, icol]
            iTl = np.searchsorted(T_grid, T_l) - 1
            iTl = max(0, min(iTl, nT - 2))
            fTl = (T_l - T_grid[iTl]) / (T_grid[iTl + 1] - T_grid[iTl])
            fTl = max(0.0, min(1.0, fTl))
            layer_planck = sigma * T_l ** 4
            for ib in range(nband):
                for igp in range(ngpt):
                    g_idx_orig = igp % ngpt_orig if (is_esft) else igp
                    ib_orig = min(ib, nband_orig - 1)
                    frac_l = (planck_frac[ib_orig, g_idx_orig, iTl] * (1 - fTl)
                              + planck_frac[ib_orig, g_idx_orig, iTl + 1] * fTl)
                    planck_src[ib, igp, k, icol] = frac_l * layer_planck
    return planck_src, surf_src


def test_planck_sources_kernel_matches_reference():
    from climt._components.picket_fence.lw.kernels import planck_sources_kernel

    rng = np.random.RandomState(0)
    nband, ngpt, nT, nlev, ncol = 14, 8, 12, 30, 5
    planck_frac = rng.uniform(0.0, 1.0, size=(nband, ngpt, nT)).astype(np.float32)
    T_grid = np.linspace(50.0, 1000.0, nT).astype(np.float64)
    T = rng.uniform(150.0, 320.0, size=(nlev, ncol))
    T_surf = rng.uniform(250.0, 320.0, size=ncol)
    sigma = 5.670374419e-8

    ps_ref, ss_ref = _reference(
        planck_frac, T_grid, T, T_surf, sigma, nband, ngpt,
        False, ngpt, nband)

    ps = np.zeros((nband, ngpt, nlev, ncol))
    ss = np.zeros((nband, ngpt, ncol))
    planck_sources_kernel(
        planck_frac, T_grid, T, T_surf, sigma,
        nband, ngpt, False, ngpt, nband, ps, ss)

    np.testing.assert_allclose(ps, ps_ref, rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(ss, ss_ref, rtol=1e-6, atol=1e-12)


def test_planck_sources_kernel_clamps_out_of_range_T():
    """T below/above the grid clamps the fraction, matching the oracle."""
    from climt._components.picket_fence.lw.kernels import planck_sources_kernel

    rng = np.random.RandomState(1)
    nband, ngpt, nT, nlev, ncol = 4, 2, 6, 8, 3
    planck_frac = rng.uniform(0.0, 1.0, size=(nband, ngpt, nT)).astype(np.float32)
    T_grid = np.linspace(200.0, 300.0, nT).astype(np.float64)
    # Deliberately out of [200,300] at both ends.
    T = rng.uniform(120.0, 360.0, size=(nlev, ncol))
    T_surf = rng.uniform(120.0, 360.0, size=ncol)
    sigma = 5.670374419e-8

    ps_ref, ss_ref = _reference(
        planck_frac, T_grid, T, T_surf, sigma, nband, ngpt,
        False, ngpt, nband)
    ps = np.zeros((nband, ngpt, nlev, ncol))
    ss = np.zeros((nband, ngpt, ncol))
    planck_sources_kernel(
        planck_frac, T_grid, T, T_surf, sigma,
        nband, ngpt, False, ngpt, nband, ps, ss)
    np.testing.assert_allclose(ps, ps_ref, rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(ss, ss_ref, rtol=1e-6, atol=1e-12)

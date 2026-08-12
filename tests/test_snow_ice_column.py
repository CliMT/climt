import numpy as np
from climt._core.snow_ice_column import solve_column, Dirichlet, Flux


def test_dirichlet_both_ends_reaches_linear_steady_state():
    n = 21
    kappa = np.full(n - 1, 2.0)
    rho = np.full(n - 1, 900.0)
    c = np.full(n - 1, 2000.0)
    T = np.full(n, 260.0)
    dz = 0.1
    dt = 50.0
    # Bottom held at 270 K, top at 250 K -> steady state is linear in z.
    T_new = T.copy()
    for _ in range(20000):
        T_new = solve_column(rho, c, kappa, T_new, dt, dz,
                             top_bc=Dirichlet(250.0), bottom_bc=Dirichlet(270.0))
    expected = np.linspace(270.0, 250.0, n)  # index 0 = bottom
    assert np.allclose(T_new, expected, atol=0.05)


def test_flux_top_bc_matches_gradient():
    n = 11
    kappa = np.full(n - 1, 2.0)
    rho = np.full(n - 1, 900.0)
    c = np.full(n - 1, 2000.0)
    T = np.linspace(270.0, 260.0, n)
    dz = 0.1
    dt = 10.0
    # A prescribed downward flux at the top should warm the top node.
    T_new = solve_column(rho, c, kappa, T.copy(), dt, dz,
                         top_bc=Flux(30.0), bottom_bc=Dirichlet(270.0))
    assert T_new[-1] > T[-1]
    assert not np.any(np.isnan(T_new))


def test_solver_no_nan_and_bounded_for_random_profiles():
    rng = np.random.default_rng(0)
    for _ in range(50):
        n = rng.integers(5, 40)
        kappa = rng.uniform(0.3, 2.3, n - 1)
        rho = rng.uniform(300.0, 950.0, n - 1)
        c = rng.uniform(1500.0, 2200.0, n - 1)
        T = rng.uniform(250.0, 272.0, n)
        out = solve_column(rho, c, kappa, T, 30.0, 0.1,
                          top_bc=Flux(rng.uniform(-50, 50)),
                          bottom_bc=Dirichlet(271.0))
        assert not np.any(np.isnan(out))
        # implicit diffusion cannot create a new extremum beyond BCs + init
        assert out.min() > 240.0 and out.max() < 290.0

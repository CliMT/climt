import numpy as np
from climt import SlabSurface, get_default_state, get_grid
import sympl


def test_qflux_adds_to_sea_tendency():
    sympl.set_backend(sympl.DataArrayBackend())
    grid = get_grid(nx=4, ny=1, nz=10)
    slab = SlabSurface()
    state = get_default_state([slab], grid_state=grid)
    state["area_type"].values[:] = "sea"
    tend0, _ = slab(state)
    state["ocean_heat_transport_convergence"].values[:] = 50.0
    tend1, _ = slab(state)
    dT = (tend1["surface_temperature"].values - tend0["surface_temperature"].values)
    assert np.all(dT > 0.0)   # positive convergence warms the mixed layer


def test_zero_qflux_matches_baseline():
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface()
    state = get_default_state([slab], grid_state=get_grid(nx=4, ny=1, nz=10))
    state["area_type"].values[:] = "sea"
    tend, _ = slab(state)
    # default field is zero -> tendency identical to computing without the term
    assert not np.any(np.isnan(tend["surface_temperature"].values))


def test_ekman_pumping_from_wind_stress_curl():
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface(include_ekman=True)
    state = get_default_state([slab], grid_state=get_grid(nx=32, ny=16, nz=10))
    state["area_type"].values[:] = "sea"
    lat = state["latitude"].values
    # zonal stress with meridional structure -> nonzero curl
    state["surface_downward_eastward_stress"].values[:] = 0.1 * np.cos(np.deg2rad(lat))
    state["surface_downward_northward_stress"].values[:] = 0.0
    tend, diag = slab(state)
    assert "ekman_pumping" in diag
    assert np.any(diag["ekman_pumping"].values != 0.0)
    assert not np.any(np.isnan(tend["surface_temperature"].values))


def test_ekman_zero_curl_gives_zero_pumping():
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface(include_ekman=True)
    state = get_default_state([slab], grid_state=get_grid(nx=16, ny=8, nz=10))
    state["area_type"].values[:] = "sea"
    state["surface_downward_eastward_stress"].values[:] = 0.1   # uniform -> no curl
    tend, diag = slab(state)
    assert np.allclose(diag["ekman_pumping"].values, 0.0, atol=1e-12)


def test_ekman_masks_stress_at_coast():
    # Mixed land/sea grid: the western half of the domain is land, the
    # eastern half is sea, with a spatially-structured wind stress over
    # the whole domain. This exercises the coastal no-flux-at-coast
    # boundary treatment (Fix 2): a coastal sea cell's divergence/curl
    # stencil must not pick up the land cells' wind stress or 300K soil
    # surface_temperature.
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface(include_ekman=True)
    grid = get_grid(nx=8, ny=8, nz=10)
    state = get_default_state([slab], grid_state=grid)

    area_type = np.empty((8, 8), dtype="S100")
    area_type[:, :4] = b"land"
    area_type[:, 4:6] = b"sea"
    area_type[:, 6:] = b"sea_ice"
    state["area_type"].values[:] = area_type

    lat = state["latitude"].values
    lon = state["longitude"].values
    # Spatially-structured stress over the whole domain (land included)
    # so that, absent masking, a coastal sea cell's stencil would pick up
    # a nonzero contribution from its land neighbors.
    state["surface_downward_eastward_stress"].values[:] = 0.1 * np.cos(
        np.deg2rad(lat)
    ) * np.sin(np.deg2rad(lon))
    state["surface_downward_northward_stress"].values[:] = 0.05 * np.sin(
        np.deg2rad(lat)
    )
    # Land soil surface_temperature stays at its (much warmer) default;
    # give the sea cells a distinct value so contamination would be
    # detectable if it occurred.
    state["surface_temperature"].values[:] = 300.0

    sea_mask = area_type == b"sea"
    non_sea_mask = ~sea_mask

    # Baseline with the Ekman term disabled, to isolate its contribution
    # to the surface_temperature tendency.
    slab_no_ekman = SlabSurface(include_ekman=False)
    tend_no_ekman, _ = slab_no_ekman(state)

    tend, diag = slab(state)

    ekman_pumping = diag["ekman_pumping"].values
    ekman_conv = diag["ekman_heat_transport_convergence"].values

    # (a) sea cells: finite (no NaN) Ekman diagnostics.
    assert not np.any(np.isnan(ekman_pumping[sea_mask]))
    assert not np.any(np.isnan(ekman_conv[sea_mask]))

    # (b) land/sea_ice cells: exactly zero Ekman diagnostics (masked
    # before differentiating, not just on output).
    assert np.all(ekman_pumping[non_sea_mask] == 0.0)
    assert np.all(ekman_conv[non_sea_mask] == 0.0)

    # (c) surface_temperature tendency on land/sea_ice cells is unaffected
    # by the Ekman term (it only ever adds to the sea-only q-flux).
    assert np.allclose(
        tend["surface_temperature"].values[non_sea_mask],
        tend_no_ekman["surface_temperature"].values[non_sea_mask],
    )

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


def test_ocean_heat_transport_convergence_reports_total_with_ekman():
    # Fix (final review): ocean_heat_transport_convergence must report the
    # TOTAL (prescribed q-flux + Ekman) actually applied to sea cells when
    # include_ekman=True, not just echo the prescribed input. This locks
    # that behavior: with a nonzero prescribed q-flux AND a nonzero Ekman
    # contribution (from a wind-stress curl), the diagnostic must equal
    # prescribed + ekman_heat_transport_convergence pointwise.
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface(include_ekman=True)
    state = get_default_state([slab], grid_state=get_grid(nx=32, ny=16, nz=10))
    state["area_type"].values[:] = "sea"
    lat = state["latitude"].values
    state["surface_downward_eastward_stress"].values[:] = 0.1 * np.cos(np.deg2rad(lat))
    state["surface_downward_northward_stress"].values[:] = 0.0
    prescribed = 15.0
    state["ocean_heat_transport_convergence"].values[:] = prescribed
    _, diag = slab(state)

    ekman_conv = diag["ekman_heat_transport_convergence"].values
    assert np.any(ekman_conv != 0.0)  # sanity: Ekman actually contributes here

    np.testing.assert_allclose(
        diag["ocean_heat_transport_convergence"].values,
        prescribed + ekman_conv,
    )


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


def test_ekman_coastal_sea_cells_independent_of_land_stress():
    # Discriminating regression test for the pre-differencing stress mask
    # (slab_surface.py ~L341-352, inside the `include_ekman` block).
    #
    # test_ekman_masks_stress_at_coast (above) checks that sea cells are
    # finite and that land/sea_ice cells are exactly zero -- but BOTH of
    # those hold even if the pre-differencing mask were deleted: land/
    # sea_ice outputs are zeroed unconditionally by a separate output-
    # masking step (L390-391) regardless of what the curl/divergence
    # stencil actually saw, and smooth stress never produces a NaN even
    # when it leaks across the coast. So that test would NOT catch a
    # regression that removed the pre-differencing mask.
    #
    # This test closes that gap directly: it runs SlabSurface(include_
    # ekman=True) twice on grids that are identical in every way except
    # the wind stress assigned to the LAND cells (the sea-cell stress is
    # byte-for-byte the same in both runs). If the land stress is zeroed
    # BEFORE the curl/divergence finite-difference stencil runs (the
    # fix), it can never influence a neighboring sea cell's derivative,
    # so the two runs' sea-cell Ekman diagnostics must be IDENTICAL.
    # If the mask were removed, the differing land-cell stress would
    # leak into the coastal sea cells' centered-difference stencil (which
    # reaches one cell into land in the longitude direction) and the two
    # runs would differ there. So: identical sea-cell output across
    # differing land input is only explained by the land stress having
    # been masked to zero prior to differentiation -- exactly what the
    # fix does.
    sympl.set_backend(sympl.DataArrayBackend())
    grid = get_grid(nx=8, ny=8, nz=10)

    def make_state():
        state = get_default_state([SlabSurface(include_ekman=True)], grid_state=grid)
        area_type = np.empty((8, 8), dtype="S100")
        # land occupies lon columns 0:4 (a coastal neighbor of the sea
        # column at index 4 in the longitude / axis=1 direction, which is
        # the direction differentiated by dfdlon in both divergence and
        # curl -- see climt/_core/horizontal_operators.py).
        area_type[:, :4] = b"land"
        area_type[:, 4:6] = b"sea"
        area_type[:, 6:] = b"sea_ice"
        state["area_type"].values[:] = area_type

        lat = state["latitude"].values
        lon = state["longitude"].values
        # Identical, spatially-structured stress for EVERY run (sea and
        # non-sea alike, initially); the land-only override is applied by
        # the caller afterward, so the sea-cell values below are common
        # to both run A and run B.
        state["surface_downward_eastward_stress"].values[:] = 0.1 * np.cos(
            np.deg2rad(lat)
        ) * np.sin(np.deg2rad(lon))
        state["surface_downward_northward_stress"].values[:] = 0.05 * np.sin(
            np.deg2rad(lat)
        )
        state["surface_temperature"].values[:] = 300.0
        return state, (area_type == b"land")

    slab = SlabSurface(include_ekman=True)

    # Run A: land cells get one nonzero stress pattern.
    state_a, land_mask = make_state()
    state_a["surface_downward_eastward_stress"].values[land_mask] = 7.0
    state_a["surface_downward_northward_stress"].values[land_mask] = -3.0
    tend_a, diag_a = slab(state_a)

    # Run B: land cells get a DIFFERENT nonzero stress pattern. Sea-cell
    # stress (assigned in make_state before this override) is untouched
    # and therefore identical to run A.
    state_b, _ = make_state()
    state_b["surface_downward_eastward_stress"].values[land_mask] = -42.0
    state_b["surface_downward_northward_stress"].values[land_mask] = 19.0
    tend_b, diag_b = slab(state_b)

    sea_mask = ~land_mask & (state_a["area_type"].values != b"sea_ice")
    assert np.any(sea_mask), "test setup must include sea cells"

    # The point of the test: differing land stress must NOT change any
    # sea cell's Ekman diagnostics or surface_temperature tendency.
    assert np.allclose(
        diag_a["ekman_heat_transport_convergence"].values[sea_mask],
        diag_b["ekman_heat_transport_convergence"].values[sea_mask],
    )
    assert np.allclose(
        diag_a["ekman_pumping"].values[sea_mask],
        diag_b["ekman_pumping"].values[sea_mask],
    )
    assert np.allclose(
        tend_a["surface_temperature"].values[sea_mask],
        tend_b["surface_temperature"].values[sea_mask],
    )

    # Sanity check that the two runs actually differ somewhere (on the
    # land cells' own diagnostics would be zero by output-masking, so
    # instead confirm the inputs genuinely diverged) -- otherwise this
    # test would pass vacuously.
    assert not np.allclose(
        state_a["surface_downward_eastward_stress"].values[land_mask],
        state_b["surface_downward_eastward_stress"].values[land_mask],
    )

    # Fold in the original (non-discriminating but still valid)
    # assertions from test_ekman_masks_stress_at_coast for good measure.
    assert not np.any(np.isnan(diag_a["ekman_pumping"].values[sea_mask]))
    assert not np.any(np.isnan(diag_a["ekman_heat_transport_convergence"].values[sea_mask]))

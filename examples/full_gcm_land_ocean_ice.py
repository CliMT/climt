"""Full moist GCM with realistic land / ocean / ice boundary conditions.

Couples the GFS spectral dynamical core to a realistic Earth surface built from
climt's bundled boundary data:

* ``LandMask``      -> area_type + surface_geopotential (orography the dycore
                       reads) + land_ice_thickness, all from the bundled 2 deg
                       ETOPO/Natural-Earth files.
* ``DataOcean``     -> prescribed, month-looping climatological SST on sea cells.
* ``BucketHydrology`` -> land surface energy + soil-moisture balance.
* Static ice        -> land_ice cells are *not* evolved; they only raise the
                       surface albedo and (through surface_geopotential) the
                       orography. Their surface temperature is held fixed.
* ``SimpleBoundaryLayer`` -> boundary-layer diffusion (Frierson 2006), run with
                       ``surface_fluxes='external'`` so it applies the surface
                       heat and moisture fluxes ``BucketHydrology`` computes and
                       adds its own bulk momentum drag.
* ``RRTMGLongwave`` / ``RRTMGShortwave`` -> radiation.
* ``EmanuelConvection`` + ``GridScaleCondensation`` -> moist processes.
* ``Instellation``  -> zenith angle / insolation.

By default this runs a short *smoke* integration to confirm everything couples
and stays finite. Pass ``--long`` (or ``--steps N``) for a real integration
toward a spun-up climate.

    python examples/full_gcm_land_ocean_ice.py            # short smoke test
    python examples/full_gcm_land_ocean_ice.py --long     # long climate run
    python examples/full_gcm_land_ocean_ice.py --steps 500 --plot
"""
import argparse
import os
from datetime import timedelta

import numpy as np

import climt
from gfs_dynamical_core import GFSDynamicalCore
from sympl import TimeDifferencingWrapper, UpdateFrequencyWrapper


# Per-area-type surface shortwave albedo for the static surface.
_ALBEDO = {"sea": 0.07, "land": 0.20, "land_ice": 0.60}


def build_model(nx, ny, timestep, radiation_interval):
    """Construct the dycore (with its physics list) and the setup-only pieces."""
    mask = climt.LandMask()                      # geography + topography (once)
    ocean = climt.DataOcean(_bundled_sst())      # prescribed SST (each step)
    insolation = climt.Instellation()            # zenith angle (each step)

    # BucketHydrology, SimpleBoundaryLayer and GridScaleCondensation are
    # Steppers; GFSDynamicalCore's component list accepts only
    # TendencyComponent / ImplicitTendencyComponent, so each needs
    # TimeDifferencingWrapper (UpdateFrequencyWrapper, used on the radiation
    # components below, preserves component type and would not do).
    land = TimeDifferencingWrapper(climt.BucketHydrology())
    # 'external' so the bucket's own surface fluxes drive the atmosphere rather
    # than the boundary layer computing a second, inconsistent set. Inside the
    # dycore this costs a one-step lag: sympl's composite hands every component
    # the same input state, so the boundary layer reads the previous step's
    # fluxes (zero on step 0).
    boundary_layer = TimeDifferencingWrapper(
        climt.SimpleBoundaryLayer(surface_fluxes="external")
    )
    convection = climt.EmanuelConvection()
    condensation = TimeDifferencingWrapper(climt.GridScaleCondensation())
    radiation_lw = UpdateFrequencyWrapper(
        climt.RRTMGLongwave(), radiation_interval)
    radiation_sw = UpdateFrequencyWrapper(
        climt.RRTMGShortwave(), radiation_interval)

    dycore = GFSDynamicalCore(
        [land, radiation_sw, radiation_lw, convection, condensation,
         boundary_layer],
        number_of_damped_levels=5,
    )
    grid = climt.get_grid(nx=nx, ny=ny)
    # Only the dycore seeds the state. LandMask/DataOcean/Instellation declare
    # latitude with dims ['*'] while GFSDynamicalCore declares it ['lat', 'lon'],
    # and get_default_state cannot combine the two. They are diagnostic-only and
    # their own inputs (latitude, longitude, time, area_type) are either already
    # in the dycore's state or stamped in by set_boundary_conditions.
    state = climt.get_default_state([dycore], grid_state=grid)
    return state, dycore, mask, ocean, insolation


def _bundled_sst():
    return os.path.join(os.path.dirname(climt.__file__), "_data",
                        "data_ocean", "earth_sst_climatology_1deg.nc")


def set_boundary_conditions(state, mask):
    """Stamp geography + topography once and set the static surface fields."""
    state.update(mask(state))                    # area_type, surface_geopotential,
                                                 # land_ice_thickness
    area = state["area_type"].values.astype(str)

    # static shortwave albedo per surface type (land_ice is bright)
    albedo = np.full(area.shape, _ALBEDO["sea"])
    albedo[area == "land"] = _ALBEDO["land"]
    albedo[area == "land_ice"] = _ALBEDO["land_ice"]
    for name in ("surface_albedo_for_direct_shortwave",
                 "surface_albedo_for_diffuse_shortwave"):
        state[name].values[:] = albedo

    # an initial surface temperature; DataOcean overwrites sea cells each step,
    # BucketHydrology evolves land cells, land_ice cells stay at this value.
    lat = np.radians(state["latitude"].values)
    state["surface_temperature"].values[:] = 290.0 - 40.0 * np.sin(lat) ** 2
    # a modest, non-zero surface humidity so latent exchange is active
    state["surface_specific_humidity"].values[:] = 0.006
    # small random seed on the wind to break symmetry
    state["eastward_wind"].values[:] = 1e-2 * np.random.randn(
        *state["eastward_wind"].shape)


def step(state, dycore, ocean, insolation, timestep):
    """Advance the model one dynamics step."""
    state.update(insolation(state))     # zenith angle for this time
    state.update(ocean(state))          # prescribe looping SST on sea cells
    diagnostics, state = dycore(state, timestep)
    state.update(diagnostics)
    state["time"] += timestep
    return state


def make_monitor():
    from sympl import PlotFunctionMonitor

    def plot_function(fig, state):
        ax = fig.add_subplot(2, 2, 1)
        state["surface_temperature"].plot.contourf(ax=ax, levels=16, robust=True)
        ax.set_title("Surface temperature")
        ax = fig.add_subplot(2, 2, 2)
        state["air_temperature"].mean(dim="lon").plot.contourf(
            ax=ax, levels=16, robust=True)
        ax.set_title("Zonal-mean temperature")
        ax = fig.add_subplot(2, 2, 3)
        state["eastward_wind"].mean(dim="lon").plot.contourf(
            ax=ax, levels=16, robust=True)
        ax.set_title("Zonal-mean zonal wind")
        ax = fig.add_subplot(2, 2, 4)
        state["surface_geopotential"].plot.contourf(ax=ax, levels=16)
        ax.set_title("Orography (surface geopotential)")
        fig.tight_layout()

    return PlotFunctionMonitor(plot_function)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=24,
                        help="number of dynamics steps (default: 24, a smoke test)")
    parser.add_argument("--long", action="store_true",
                        help="run a long climate integration (~200 model days)")
    parser.add_argument("--timestep", type=float, default=600.0,
                        help="dynamics timestep in seconds (default: 600)")
    parser.add_argument("--radiation-every", type=int, default=6,
                        help="call radiation every N steps (default: 6)")
    parser.add_argument("--nx", type=int, default=128)
    parser.add_argument("--ny", type=int, default=62)
    parser.add_argument("--plot", action="store_true",
                        help="show a live PlotFunctionMonitor")
    parser.add_argument("--plot-every", type=int, default=6)
    args = parser.parse_args()

    timestep = timedelta(seconds=args.timestep)
    radiation_interval = args.radiation_every * timestep
    steps = int(200 * 24 * 3600 / args.timestep) if args.long else args.steps

    state, dycore, mask, ocean, insolation = build_model(
        args.nx, args.ny, timestep, radiation_interval)
    set_boundary_conditions(state, mask)

    monitor = make_monitor() if args.plot else None
    print("Running %d steps (dt=%gs) on a %dx%d grid ..."
          % (steps, args.timestep, args.nx, args.ny))

    for i in range(steps):
        state = step(state, dycore, ocean, insolation, timestep)

        surf = state["surface_temperature"].values
        air = state["air_temperature"].values
        wind = state["eastward_wind"].values
        if not (np.isfinite(surf).all() and np.isfinite(air).all()
                and np.isfinite(wind).all()):
            raise RuntimeError("non-finite field at step %d (%s)"
                               % (i, state["time"]))

        if i % max(1, args.plot_every) == 0:
            print("step %4d  %s  T_surf[%.1f, %.1f] K  |u|max %.1f m/s"
                  % (i, state["time"], surf.min(), surf.max(),
                     np.abs(wind).max()))
            if monitor is not None:
                monitor.store(state)

    print("Done. Everything ran and stayed finite.")


if __name__ == "__main__":
    main()

"""Step a column forward in time, and draw what happened.

Every page from 07 on runs one of these two loops. They are the single copy:
the deleted ``docs/radiative-transfer/09-live-rce.qmd`` kept one version in a
Quarto include and a testable twin in a ``.py`` beside it, with a docstring
contract asking the two to be edited together and a test enforcing it. Pages
here import this module through ``pyodide: resources:`` instead, so there is
nothing to keep in sync.

Two component categories, and only two:

* **tendency components** (``CorkLongwaveRadiation``, ``SlabSurface``,
  ``EmanuelConvectionPython``) return rates of change, and are wrapped in a
  ``sympl.AdamsBashforth`` time integrator that turns them into a new state;
* **stepper components** (``SimpleBoundaryLayer``,
  ``DryConvectiveAdjustment``, ``GridScaleCondensation``) return the new state
  themselves, and are called directly.

``EmanuelConvectionPython`` is an ``ImplicitTendencyComponent``, a third sympl
protocol, and putting one inside a ``TendencyStepper`` makes sympl warn that
results "may lead to scientifically invalid results". It goes in the tendency
list anyway -- that is how climt has always run it, and page 12 explains the
warning rather than restructuring the loop around it. So: two categories.

**Update order is load-bearing.** The new prognostic state is applied *before*
the diagnostics. Reverse it and the integrator carries flux diagnostics forward
at their pre-step value, ``SlabSurface`` consumes stale longwave fluxes, and
the surface heats without bound. This broke the demo once already.

**Timestep type is load-bearing.** Pass ``climt.UnytTimeDelta``, not
``datetime.timedelta``. Under ``UnytBackend`` a bare ``timedelta`` makes
``total_seconds()`` return a unitless float, which will not cancel a tendency's
``/s``, and sympl raises ``UnitOperationError``.
"""
import numpy as np
import sympl
from sympl import AdamsBashforth

DAY_SECONDS = 86400.0


def integrate(tendency_components, stepper_components, state, timestep,
              n_steps):
    """Step ``state`` forward ``n_steps`` times, in place.

    Args:
        tendency_components: components returning tendencies; wrapped in one
            ``AdamsBashforth``.
        stepper_components: components returning a new state; called in the
            order given, after the tendency update, within each step.
        state: a climt state dict. Modified in place.
        timestep: a ``climt.UnytTimeDelta``.
        n_steps: number of steps.

    Returns:
        The same ``state`` object, for convenience.
    """
    model = AdamsBashforth(*tendency_components)
    for _ in range(n_steps):
        diagnostics, new_state = model(state, timestep)
        state.update(new_state)
        state.update(diagnostics)
        for stepper in stepper_components:
            stepper_diagnostics, stepper_state = stepper(state, timestep)
            state.update(stepper_state)
            state.update(stepper_diagnostics)
        state["time"] += timestep
    return state


def integrate_with_history(tendency_components, stepper_components, state,
                           timestep, n_steps, n_snapshots=8):
    """``integrate``, recording the evolution as it goes.

    Two scalars are recorded every step -- they cost nothing -- and the full
    profiles at ``n_snapshots`` evenly spaced steps, which do not. The first
    snapshot is taken after step 1 rather than before step 1, so the radiation
    diagnostics it reads are populated.

    Returns:
        ``(state, history)``. ``history`` is a plain dict of numpy arrays and
        a list of snapshot dicts -- no sympl objects, so it survives being
        passed around a page and pickles cleanly:

        * ``days`` -- (n_steps,) elapsed time in days
        * ``olr`` -- (n_steps,) upwelling longwave flux at the top interface
        * ``surface_temperature`` -- (n_steps,) K
        * ``snapshots`` -- list of dicts with ``day``, ``T`` (mid levels),
          ``H`` (longwave heating rate, K/day), ``U`` and ``D`` (up and
          downwelling longwave flux on interface levels)
    """
    model = AdamsBashforth(*tendency_components)
    day = timestep.total_seconds() / DAY_SECONDS
    # total_seconds() returns a unyt quantity under UnytBackend; the history is
    # plain numbers, so strip the units here rather than in six page cells.
    day = float(day)

    snap_at = set(np.unique(
        np.linspace(1, n_steps, n_snapshots).round().astype(int)).tolist())
    snapshots = []
    days, olr, surface_temperature = [], [], []

    for step in range(1, n_steps + 1):
        diagnostics, new_state = model(state, timestep)
        state.update(new_state)
        state.update(diagnostics)
        for stepper in stepper_components:
            stepper_diagnostics, stepper_state = stepper(state, timestep)
            state.update(stepper_state)
            state.update(stepper_diagnostics)
        state["time"] += timestep

        days.append(step * day)
        olr.append(float(
            state["upwelling_longwave_flux_in_air"].values[-1, 0, 0]))
        surface_temperature.append(float(
            state["surface_temperature"].values.ravel()[0]))
        if step in snap_at:
            snapshots.append(dict(
                day=step * day,
                T=state["air_temperature"].values[:, 0, 0].copy(),
                H=state["air_temperature_tendency_from_longwave"]
                .values[:, 0, 0].copy(),
                U=state["upwelling_longwave_flux_in_air"]
                .values[:, 0, 0].copy(),
                D=state["downwelling_longwave_flux_in_air"]
                .values[:, 0, 0].copy(),
            ))

    history = dict(
        days=np.array(days),
        olr=np.array(olr),
        surface_temperature=np.array(surface_temperature),
        snapshots=snapshots,
    )
    return state, history


def draw_evolution(history, state, solar, title=""):
    """Render one annotated 2x2 figure of the integration's evolution.

    quarto-live runs each Pyodide cell to completion in a Web Worker and paints
    only the cell's *final* figure, so a frame-by-frame animation is impossible
    here -- an earlier version tried and only ever showed the last frame.
    Instead the whole run happens quietly behind quarto-live's spinner and this
    draws the evolution at the end: temperature, longwave heating rate and
    up/downwelling flux profiles at several times, plus OLR and surface
    temperature as time series approaching energy balance.

    Args:
        history: the dict returned by :func:`integrate_with_history`.
        state: the state it was recorded from -- read only for the pressure
            coordinates.
        solar: prescribed absorbed shortwave at the surface, W m^-2, drawn as
            the equilibrium target on the OLR panel.
        title: figure suptitle.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    snaps = history["snapshots"]
    days = history["days"]
    olr = history["olr"]
    tsurf = history["surface_temperature"]
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    p_int = state["air_pressure_on_interface_levels"].values[:, 0, 0] / 100.0

    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=snaps[0]["day"], vmax=snaps[-1]["day"])

    fig = plt.figure(figsize=(11.5, 8.4))
    # right=0.88 reserves a margin for the colorbar. Letting colorbar(ax=[...])
    # steal the space instead does not work: its length comes from the default
    # aspect, so the bar is centred at half height and the panels it did not
    # shrink -- the OLR panel and its twinned temperature axis -- run their
    # labels straight through it.
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30, right=0.88)
    ax_T = fig.add_subplot(gs[0, 0])
    ax_H = fig.add_subplot(gs[0, 1], sharey=ax_T)
    ax_F = fig.add_subplot(gs[1, 0], sharey=ax_T)
    ax_O = fig.add_subplot(gs[1, 1])

    for k, snap in enumerate(snaps):
        colour = cmap(norm(snap["day"]))
        final = (k == len(snaps) - 1)
        width = 2.6 if final else 1.3
        order = 5 if final else 2
        ax_T.plot(snap["T"], p, color=colour, lw=width, zorder=order)
        ax_H.plot(snap["H"], p, color=colour, lw=width, zorder=order)
        ax_F.plot(snap["U"], p_int, color=colour, lw=width, zorder=order)
        ax_F.plot(snap["D"], p_int, color=colour, lw=width, ls="--",
                  zorder=order)

    ax_T.set_ylabel("pressure (hPa)")
    ax_F.set_ylabel("pressure (hPa)")
    ax_T.set_xlabel("temperature (K)")
    ax_T.set_title("Temperature profile")
    ax_H.set_xlabel("LW heating rate (K day$^{-1}$)")
    ax_H.set_title("Longwave heating rate")
    ax_F.set_xlabel("LW flux (W m$^{-2}$)")
    ax_F.set_title("LW flux: up (solid) / down (dashed)")
    # Log-pressure y axis, always. On a linear axis the entire upper atmosphere
    # is crushed into the top sliver of the panel: the gray column's isothermal
    # layer spans ~20 -> 4 hPa, which is 1.6% of a linear axis and ~19% of this
    # one. Limits run over the *interface* levels so the flux panel's topmost
    # point is not clipped; high -> low inverts, and the shared y applies that
    # to all three profile panels.
    ax_T.set_yscale("log")
    ax_T.set_ylim(p_int.max(), p_int.min())
    for ax in (ax_T, ax_H, ax_F):
        ax.grid(alpha=0.3, which="both")
        ax.yaxis.set_major_formatter(ScalarFormatter())   # 1000, not 10^3
    ax_H.axvline(0.0, color="0.4", lw=0.8)

    T_final = snaps[-1]["T"]
    ax_T.annotate("warm surface", xy=(T_final[0], p[0]),
                  xytext=(0.42, 0.14), textcoords="axes fraction", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="0.3"))
    # Only claim a skin temperature when the profile has one. The gray column
    # relaxes to a uniform top; the non-grey column keeps cooling to the model
    # lid, because its topmost layers absorb only in the strong bands, where
    # the upwelling flux was already emitted by cold levels below. Labelling
    # that "isothermal" would be a caption contradicting its own figure.
    isothermal_top = abs(T_final[-1] - T_final[-2]) < 2.0
    ax_T.annotate("isothermal top\n(skin temperature)" if isothermal_top
                  else "no isothermal top:\nstill cooling at the model lid",
                  xy=(T_final[-1], p[-1]),
                  xytext=(0.30, 0.72), textcoords="axes fraction", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="0.3"))
    ax_H.text(0.04, 0.06,
              "cooling (< 0) aloft\nrelaxes toward 0\nat equilibrium",
              transform=ax_H.transAxes, fontsize=8, va="bottom")
    ax_F.text(0.55, 0.9, "surface emits up;\nspace sends nothing down",
              transform=ax_F.transAxes, fontsize=8, va="top")

    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    cax = fig.add_axes([0.945, 0.30, 0.013, 0.42])
    fig.colorbar(mappable, cax=cax).set_label(
        "time (days)  —  early (dark) → late (yellow)")

    ax_O.plot(days, olr, color="#c92a2a", lw=1.9,
              label="OLR (top of atmosphere)")
    ax_O.axhline(solar, color="0.35", ls="--", lw=1.2,
                 label=f"absorbed solar = {solar:.0f} W m$^{{-2}}$")
    ax_O.set_xlabel("time (days)")
    ax_O.set_ylabel("longwave flux (W m$^{-2}$)", color="#c92a2a")
    ax_O.tick_params(axis="y", labelcolor="#c92a2a")
    ax_O.set_title("Approach to energy balance")
    ax_O.grid(alpha=0.3)

    ax_Ts = ax_O.twinx()
    ax_Ts.plot(days, tsurf, color="#1c7ed6", lw=1.6)
    ax_Ts.set_ylabel("surface temperature (K)", color="#1c7ed6")
    ax_Ts.tick_params(axis="y", labelcolor="#1c7ed6")

    gap = solar - olr[-1]
    message = (f"still {gap:+.0f} W m$^{{-2}}$ from balance\n"
               "(raise n_steps to close it)"
               if abs(gap) > 2 else
               "OLR ≈ absorbed solar:\nnear equilibrium")
    ax_O.annotate(message, xy=(days[-1], olr[-1]), xytext=(0.12, 0.30),
                  textcoords="axes fraction", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="0.3"))
    ax_O.legend(loc="lower right", fontsize=7)

    if title:
        fig.suptitle(title, fontsize=13, y=0.98)
    # Required in the browser: quarto-live paints the figure a cell shows, and
    # a cell that only creates one paints nothing.
    plt.show()


class UnytRelaxation(sympl.RelaxationTendencyComponent):
    """``sympl.RelaxationTendencyComponent`` with unyt-parseable tendency units.

    The stock component derives its tendency units with pint --
    ``str(ureg("m s^-1") / ureg("s"))`` is ``"1.0 meter / second ** 2"`` -- and
    unyt cannot parse that string, so constructing a state around it under
    ``climt.UnytBackend`` raises ``ValueError: Failed to parse unit``. The
    relaxation arithmetic itself is fine; only the units label is wrong for
    this backend. So: take the tendency units as an argument and hand unyt
    something it understands.
    """

    def __init__(self, quantity_name, units, tendency_units, **kwargs):
        self._tendency_units = tendency_units
        super(UnytRelaxation, self).__init__(quantity_name, units, **kwargs)

    @property
    def tendency_properties(self):
        return {self._quantity_name: {"dims": ["*"],
                                      "units": self._tendency_units}}


def wind_relaxation(state, speed, timescale_hours=24.0,
                    quantity_name="eastward_wind", initialise=True):
    """Give a column a wind, and keep it -- returns a tendency component.

    A single column has no dynamics, so nothing drives a wind, and
    ``SimpleBoundaryLayer``'s surface drag removes any wind you prescribe
    within days. Every configuration from page 8 on therefore needs a momentum
    source, or its surface fluxes are those of a dead-calm column.

    climt's examples do this by assigning the wind back at the end of every
    timestep (``examples/column_code_with_slab.py``). This is the same idea as
    a component: ``dx/dt = -(x - x_eq)/tau``, which composes with everything
    else, carries its units, and puts the assumption in the component list
    where a reader can see it.

    Writes ``equilibrium_<quantity>`` and ``<quantity>_relaxation_timescale``
    into ``state`` -- ``get_default_state`` does not know about them, because
    they belong to the relaxation rather than to any physics component.

    Args:
        state: the state to add the target fields to. Modified in place.
        speed: the wind the column relaxes toward, m/s.
        timescale_hours: relaxation timescale. Long compared with a timestep
            and short compared with the run, or it is either a hard reset or
            no forcing at all.
        quantity_name: the wind component to relax.
        initialise: also set the wind itself to ``speed``, which is what you
            want when building a fresh state. Pass ``False`` when attaching
            this to a state loaded from disk -- that state already carries a
            spun-up wind profile, sheared by the drag, and flattening it to a
            uniform ``speed`` would throw away part of the equilibrium.

    Returns:
        A tendency component to put in the tendency list.
    """
    backend = sympl.get_backend()
    template = state[quantity_name]

    def _field(name, value, units):
        return backend.create_quantity(
            np.full(template.values.shape, value), name, units, template.dims)

    state["equilibrium_" + quantity_name] = _field(
        "equilibrium_" + quantity_name, speed, "m s^-1")
    state[quantity_name + "_relaxation_timescale"] = _field(
        quantity_name + "_relaxation_timescale", timescale_hours * 3600.0, "s")
    if initialise:
        state[quantity_name].values[:] = speed

    return UnytRelaxation(quantity_name, "m s^-1", "m s^-2")

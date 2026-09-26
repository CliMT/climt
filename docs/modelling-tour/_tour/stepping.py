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
              n_steps, after_step=None):
    """Step ``state`` forward ``n_steps`` times, in place.

    Args:
        tendency_components: components returning tendencies; wrapped in one
            ``AdamsBashforth``.
        stepper_components: components returning a new state; called in the
            order given, after the tendency update, within each step.
        state: a climt state dict. Modified in place.
        timestep: a ``climt.UnytTimeDelta``.
        n_steps: number of steps.
        after_step: optional ``f(state)``, called at the very end of every
            step, after the clock has advanced -- to record something (a
            :class:`Recorder`), or to overwrite something. Page 8 uses the
            second to assign the wind back every step, the way climt's
            examples do.

            Use this rather than calling ``integrate(..., 1)`` in a loop of
            your own. Each call builds a fresh ``AdamsBashforth``, which has
            no memory of earlier steps, so a loop of one-step calls silently
            runs forward Euler instead of the third-order scheme.

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
        if after_step is not None:
            after_step(state)
    return state


class Recorder:
    """Record scalars off the state every step -- an ``after_step`` for
    :func:`integrate`.

    Give it the timestep, and one function per series, each taking a state and
    returning a number; it keeps the elapsed time in days alongside them::

        record = Recorder(timestep, jump=budgets.surface_air_jump,
                          flux=budgets.sensible_heat_flux)
        integrate(tendencies, steppers, state, timestep, n, after_step=record)
        record["days"], record["flux"]            # numpy arrays

    ``profiles`` names state quantities to snapshot as whole columns every
    ``profile_every`` steps, the first after step 1. They come back as a 2-D
    array, one row per snapshot, with the snapshot days under
    ``"profile_days"``.

    Days count from the first step this recorder saw, so a recorder attached
    to a state that has already run for a year still starts its axis near 0.
    """

    def __init__(self, timestep, profiles=(), profile_every=24, **series):
        self._dt_days = float(timestep.total_seconds()) / DAY_SECONDS
        self._series = series
        self._profile_every = profile_every
        self._values = {name: [] for name in series}
        self._values["days"] = []
        self._snapshots = {name: [] for name in profiles}
        self._snapshot_days = []
        self._step = 0

    def __call__(self, state):
        self._step += 1
        day = self._step * self._dt_days
        self._values["days"].append(day)
        for name, read in self._series.items():
            self._values[name].append(float(read(state)))
        if self._snapshots and (self._step - 1) % self._profile_every == 0:
            self._snapshot_days.append(day)
            for name in self._snapshots:
                values = np.asarray(state[name].values, dtype=float)
                self._snapshots[name].append(
                    values.reshape(values.shape[0], -1)[:, 0].copy())

    def __getitem__(self, name):
        if name == "profile_days":
            return np.array(self._snapshot_days)
        if name in self._snapshots:
            return np.array(self._snapshots[name])
        return np.array(self._values[name])

    def mean(self, name, days=30.0):
        """Mean of series ``name`` over its last ``days`` days."""
        elapsed = self["days"]
        return float(np.mean(self[name][elapsed > elapsed[-1] - days]))


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

    ``T`` is the state *after* the step, while ``H``, ``U`` and ``D`` come from
    the diagnostics the components returned, which were evaluated at the
    start-of-step state -- so the fluxes lag the temperature by one step. That
    is the design and not a bug: it is the radiation that produced this
    temperature. Do not "fix" it by reordering the loop, which would apply
    diagnostics before the new state and break the update order above.

    ``n_snapshots`` is a ceiling rather than a count: the snapshot steps are
    rounded onto integers and de-duplicated, so asking for more snapshots than
    there are steps yields one per step and no more.
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


def draw_discontinuity(radiative_state, turbulent_state, record,
                       depth_hpa=200.0, title=""):
    """Page 8's headline: the surface--air jump, before and after turbulence.

    Two profile panels side by side, zoomed to the lowest ``depth_hpa`` of the
    column, each with the surface temperature drawn as a point at the surface
    pressure and the jump to the lowest model layer labelled in K; beneath
    them, the sensible heat flux the boundary layer applied, step by step
    (faint) and as daily means.

    Args:
        radiative_state: the converged column with no boundary layer.
        turbulent_state: the converged column with one.
        record: the :class:`Recorder` from the turbulent run, carrying
            ``flux`` (W m^-2) and ``days``.
        depth_hpa: how much of the column above the surface to show.
        title: figure suptitle.
    """
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11.0, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.38,
                          wspace=0.12)
    ax_rad = fig.add_subplot(gs[0, 0])
    ax_turb = fig.add_subplot(gs[0, 1], sharex=ax_rad, sharey=ax_rad)
    ax_flux = fig.add_subplot(gs[1, :])

    panels = ((ax_rad, radiative_state, "Radiative equilibrium", "#c92a2a"),
              (ax_turb, turbulent_state, "+ boundary layer", "#1c7ed6"))
    for ax, state, label, colour in panels:
        p = state["air_pressure"].values[:, 0, 0] / 100.0
        T = state["air_temperature"].values[:, 0, 0]
        ps = float(state["surface_air_pressure"].values.ravel()[0]) / 100.0
        Ts = float(state["surface_temperature"].values.ravel()[0])
        shown = p > ps - depth_hpa
        ax.plot(T[shown], p[shown], "o-", color=colour, lw=2, ms=4,
                label="air (model layers)")
        ax.plot([Ts], [ps], "s", color="k", ms=8, label="surface")
        ax.plot([T[0], Ts], [p[0], ps], ":", color="0.35", lw=1.2)
        note = f"jump {Ts - T[0]:.2f} K"
        if state is turbulent_state:
            # The drawn profile is one step of a flickering equilibrium; say
            # what it averages to, which is the number the page quotes.
            note += f"\n(30-day mean {record.mean('jump'):.2f} K)"
        ax.annotate(note, xy=(0.5 * (T[0] + Ts), 0.5 * (p[0] + ps)),
                    xytext=(0.60, 0.30), textcoords="axes fraction",
                    fontsize=10, fontweight="bold", ha="center",
                    arrowprops=dict(arrowstyle="->", color="0.3"))
        ax.set_title(label)
        ax.set_xlabel("temperature (K)")
        ax.grid(alpha=0.3)
    ax_rad.set_ylabel("pressure (hPa)")
    ps = float(radiative_state["surface_air_pressure"].values.ravel()[0])
    ax_rad.set_ylim(ps / 100.0 + 8.0, ps / 100.0 - depth_hpa)
    ax_rad.legend(loc="upper right", fontsize=8)
    plt.setp(ax_turb.get_yticklabels(), visible=False)

    days, flux = record["days"], record["flux"]
    ax_flux.plot(days, flux, color="#1c7ed6", lw=0.4, alpha=0.35,
                 label="every step")
    per_day = int(round(1.0 / (days[1] - days[0]))) if len(days) > 1 else 1
    n_days = len(flux) // per_day
    if n_days:
        daily = flux[:n_days * per_day].reshape(n_days, per_day).mean(axis=1)
        ax_flux.plot((np.arange(n_days) + 0.5) * per_day * (days[1] - days[0]),
                     daily, color="#1c7ed6", lw=1.8, label="daily mean")
    ax_flux.set_xlabel("time (days)")
    ax_flux.set_ylabel("sensible heat flux (W m$^{-2}$)")
    ax_flux.set_title("Sensible heat flux the boundary layer applied, "
                      "surface to air")
    ax_flux.grid(alpha=0.3)
    ax_flux.legend(loc="upper right", fontsize=8)

    if title:
        fig.suptitle(title, fontsize=13, y=0.99)
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


class SurfaceHumidity(sympl.Stepper):
    """Keep the surface's specific humidity at a fixed relative humidity.

    ``SimpleBoundaryLayer`` computes the latent heat flux from
    ``surface_specific_humidity``, and treats it as an *input*: nothing in
    climt's browser-safe stack writes it, so whatever value it is given stays
    put while the surface warms and cools underneath it. A fixed number is
    then a surface whose relative humidity drifts with its temperature --
    0.015 kg/kg, for instance, is saturated at about 293 K and 130 % saturated
    at 289 K. (``SimplePhysics``, which does not run in a browser, sets it
    from the surface temperature internally, which is what this does.)

    This sets ``surface_specific_humidity = relative_humidity * q_sat(Ts, ps)``
    every step, from the surface temperature the slab has just produced. Put
    it in the stepper list *before* the boundary layer, so the boundary layer
    sees the humidity of the surface it is exchanging with.

    ``q_sat`` is ``soundings.saturation_specific_humidity``, the form
    ``GridScaleCondensation`` uses, so relative humidity 1 here and in the
    air mean the same thing. A page using this must list ``_tour/soundings.py``
    in its ``pyodide: resources:``.
    """

    input_properties = {
        "surface_temperature": {"dims": ["*"], "units": "degK"},
        "surface_air_pressure": {"dims": ["*"], "units": "Pa"},
    }
    diagnostic_properties = {}
    output_properties = {
        "surface_specific_humidity": {"dims": ["*"], "units": "kg/kg"},
    }

    def __init__(self, relative_humidity=1.0, **kwargs):
        self.relative_humidity = float(relative_humidity)
        super(SurfaceHumidity, self).__init__(**kwargs)

    def array_call(self, state, timestep):
        from soundings import saturation_specific_humidity

        q_surface = self.relative_humidity * saturation_specific_humidity(
            np.asarray(state["surface_temperature"], dtype=float),
            np.asarray(state["surface_air_pressure"], dtype=float))
        return {}, {"surface_specific_humidity": q_surface}


def draw_moisture(state, record, depth_hpa=250.0, title=""):
    """Page 9's headline: what the moisture does to the column's buoyancy.

    Three panels. Left, potential temperature and virtual potential
    temperature over the lowest ``depth_hpa`` of the column, with the gap
    between them shaded -- that gap *is* the moisture's contribution to
    buoyancy. Middle, specific humidity against its saturation value over
    the same depth. Right, the sensible and latent heat fluxes the boundary
    layer applied, step by step (faint) and as daily means.

    Args:
        state: the column at the end of the run.
        record: the :class:`Recorder` from the run, carrying ``sh`` and
            ``lh`` (W m^-2) and ``days``.
        depth_hpa: how much of the column above the surface to show.
        title: figure suptitle.
    """
    import matplotlib.pyplot as plt
    from soundings import saturation_specific_humidity

    Rd = float(sympl.get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Cp = float(sympl.get_constant(
        "heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"))
    p = np.asarray(state["air_pressure"].values, dtype=float)[:, 0, 0]
    T = np.asarray(state["air_temperature"].values, dtype=float)[:, 0, 0]
    q = np.asarray(state["specific_humidity"].values, dtype=float)[:, 0, 0]
    theta = T * (1.0e5 / p) ** (Rd / Cp)
    theta_v = theta * (1.0 + 0.61 * q)
    q_sat = saturation_specific_humidity(T, p)
    ps = float(np.asarray(state["surface_air_pressure"].values).ravel()[0])
    shown = p > ps - depth_hpa * 100.0
    hpa = p / 100.0

    fig = plt.figure(figsize=(12.0, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.6], wspace=0.32)
    ax_th = fig.add_subplot(gs[0, 0])
    ax_q = fig.add_subplot(gs[0, 1], sharey=ax_th)
    ax_f = fig.add_subplot(gs[0, 2])

    ax_th.fill_betweenx(hpa[shown], theta[shown], theta_v[shown],
                        color="#1c7ed6", alpha=0.25, lw=0)
    ax_th.plot(theta[shown], hpa[shown], "o-", color="#c92a2a", ms=3, lw=1.8,
               label="θ")
    ax_th.plot(theta_v[shown], hpa[shown], "o-", color="#1c7ed6", ms=3,
               lw=1.8, label="θ$_v$ = θ(1 + 0.61q)")
    gap = theta_v[0] - theta[0]
    ax_th.annotate(f"θ$_v$ − θ = {gap:.2f} K\nat the lowest level",
                   xy=(0.5 * (theta[0] + theta_v[0]), hpa[0]),
                   xytext=(0.40, 0.30), textcoords="axes fraction",
                   fontsize=9, arrowprops=dict(arrowstyle="->", color="0.3"))
    ax_th.set(xlabel="potential temperature (K)", ylabel="pressure (hPa)",
              title="Moisture adds buoyancy")
    ax_th.set_ylim(ps / 100.0 + 5.0, ps / 100.0 - depth_hpa)
    ax_th.legend(loc="upper left", fontsize=8)
    ax_th.grid(alpha=0.3)

    ax_q.plot(q[shown] * 1e3, hpa[shown], "o-", color="#2b8a3e", ms=3,
              lw=1.8, label="specific humidity q")
    ax_q.plot(q_sat[shown] * 1e3, hpa[shown], "--", color="0.4", lw=1.2,
              label="saturation q$_{sat}$")
    ax_q.set(xlabel="specific humidity (g kg$^{-1}$)",
             title="Water the surface supplied")
    ax_q.legend(loc="upper right", fontsize=8)
    ax_q.grid(alpha=0.3)
    plt.setp(ax_q.get_yticklabels(), visible=False)

    days = record["days"]
    per_day = int(round(1.0 / (days[1] - days[0]))) if len(days) > 1 else 1
    n_days = len(days) // per_day
    for name, label, colour in (("sh", "sensible", "#c92a2a"),
                                ("lh", "latent", "#1c7ed6")):
        series = record[name]
        ax_f.plot(days, series, color=colour, lw=0.4, alpha=0.3)
        if n_days:
            daily = series[:n_days * per_day].reshape(n_days, per_day).mean(1)
            ax_f.plot((np.arange(n_days) + 0.5) * per_day * (days[1] - days[0]),
                      daily, color=colour, lw=2.0,
                      label=f"{label} (daily mean)")
    ax_f.set(xlabel="time (days)", ylabel="upward flux at the surface "
             "(W m$^{-2}$)", title="How the surface sheds its heat")
    ax_f.set_ylim(bottom=0.0)
    ax_f.legend(loc="upper right", fontsize=8)
    ax_f.grid(alpha=0.3)

    if title:
        fig.suptitle(title, fontsize=13, y=1.02)
    plt.show()

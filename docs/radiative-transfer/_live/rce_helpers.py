"""Time-integration helper for the in-browser (and native-test) RCE demo.

This module is the plain, importable, natively-testable twin of the
`integrate_to_equilibrium` definition embedded in the Pyodide boot cell of
`docs/_includes/climt-live-setup.qmd`. Keep the body of the function
IDENTICAL in both places: `tests/test_live_rce_demo.py` unit-tests this
copy, and the flagship `docs/radiative-transfer/09-live-rce.qmd` page runs
the copy embedded in the `.qmd` include (Pyodide cannot import this file
directly since `docs/radiative-transfer` contains a hyphen and is not a
valid Python package path).
"""
from sympl import AdamsBashforth


def integrate_to_equilibrium(tendency_components, stepper_components, state,
                              timestep, n_steps):
    """Step tendency + stepper components forward toward radiative equilibrium.

    Keep this body in sync with the copy in
    docs/_includes/climt-live-setup.qmd.

    Update order matters: the new (prognostic) state is applied *before* the
    diagnostics, so flux diagnostics (e.g. the longwave fluxes SlabSurface
    consumes) survive into the next step instead of being clobbered by the
    stepper carrying those fields forward at their pre-step value.
    """
    model = AdamsBashforth(list(tendency_components))
    for _ in range(n_steps):
        diagnostics, new_state = model(state, timestep)
        state.update(new_state)
        state.update(diagnostics)
        for stepper in stepper_components:
            s_diag, s_new = stepper(state, timestep)
            state.update(s_new)
            state.update(s_diag)
        state["time"] += timestep
    return state

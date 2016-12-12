from .base_components import MonitorComposite
from .timestepping import AdamsBashforth
from collections import Iterable


class Federation(object):

    def __init__(self, initial_state, prognostic=None, diagnostic=None,
                 implicit=None, monitor=None):
        """

        Args:
            initial_state (dict): A model state dictionary to use as the
                initial model state.
            prognostic (Prognostic or iterable of Prognostic, optional):
                Prognostics to use when integrating forward in time. Default is
                no prognostics.
            diagnostic (Diagnostic or iterable of Diagnostic, optional):
                Diagnostics to compute before integrating forward in time.
                Default is no diagnostics.
            implicit (Implicit, optional): Implicit scheme to use to integrate
                forward in time. Default is no implicit scheme.
            monitor (Monitor or iterable of Monitor, optional): Monitors to
                receive model state as integration is performed. Default is no
                monitors.
        """
        self._state = initial_state
        self._implicit = implicit
        if isinstance(prognostic, Iterable):
            self._prognostic_list = prognostic
        elif prognostic is None:
            self._prognostic_list = []
        else:
            self._prognostic_list = [diagnostic]
        if isinstance(diagnostic, Iterable):
            self._diagnostic_list = diagnostic
        elif diagnostic is None:
            self._diagnostic_list = []
        else:
            self._diagnostic_list = [diagnostic]
        self._time_stepper = AdamsBashforth(
            prognostic_list=self._prognostic_list,
            diagnostic_list=self._diagnostic_list,
            order=3)
        if isinstance(monitor, Iterable):
            self._monitor = MonitorComposite(monitor)
        elif monitor is None:
            self._monitor = MonitorComposite([])
        else:
            self._monitor = monitor

    def step(self):
        pass

===============
Component Types
===============

Deriving from Sympl_, most components in climt are either :py:class:`TendencyComponent`,
:py:class:`DiagnosticComponent`, :py:class:`ImplicitTendencyComponent`
or :py:class:`Stepper`.

* :py:class:`TendencyComponent` takes the model state as input and returns tendencies and optional
  diagnostics. A radiation component is a good example of such components: any radiation component
  returns the heating rate of each layer of the atmosphere (in :math:`degK/s`)
* :py:class:`DiagnosticComponent` takes the model state as input and returns some other diagnostic
  quantities. An example would be a component that takes the model state and returns the optical
  depth.
* :py:class:`Stepper` takes the model state, and returns new values for some some quantities
  in the model state. A dynamical core is a good example, which returns new values of winds,
  temperature and pressure.
* :py:class:`ImplicitTendencyComponent`  takes the model state and outputs tendencies (like a :py:class:`TendencyComponent`)
  but require the model timestep (like a :py:class:`Stepper`) for various reasons, including to ensure
  that a vertical CFL_ criterion is met.


You can read more about this schema of model components in Sympl's
documentation_.

.. _Sympl: http://sympl.readthedocs.io
.. _documentation: http://sympl.readthedocs.io/en/latest/computation.html
.. _CFL: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

=============== 
Component Types
===============

Deriving from Sympl_, most components in CliMT are either :py:class:`ClimtPrognostic`, :py:class:`ClimtDiagnostic`
or :py:class:`ClimtImplicit`. To summarise their behaviours:

* :py:class:`ClimtPrognostic` takes the model state as input and returns tendencies and optional
  diagnostics. A radiation component is a good example of such components: any radiation component
  returns the heating rate of each layer of the atmosphere (in :math:`degK/s`)
* :py:class:`ClimtDiagnostic` takes the model state as input and returns some other diagnostic
  quantities. An example would be a component that takes the model state and returns the optical
  depth.
* :py:class:`ClimtImplicit` takes the model state, and returns new values for some some quantities
  in the model state. A dynamical core is a good example, which returns new values of winds,
  temperature and pressure.

You can read more about this schema of model components in Sympl's
documentation_. 

Unfortunately, not all components used in a typical climate model fit into this
schema. Convection schemes output tendencies (like a :py:class:`ClimtPrognostic`) but require the model timestep
(like a :py:class:`ClimtImplicit`) for various reasons, including to ensure that a vertical CFL_ criterion is met.
Spectral dynamical cores require model state and tendencies in spectral space, which means they
don't play well with :py:class:`~sympl.Implicit` components which modify model state in grid space.

To account for this diversity of model components, CliMT introduces two additional entities: :py:class:`ClimtImplicitPrognostic`
which subclasses :py:class:`~sympl.Prognostic` and :py:class:`~ClimtSpectralDynamicalCore`, which subclasses :py:class:`~sympl.TimeStepper`. The awkward
name :py:class:`ImplicitPrognostic` mirrors the awkwardness of the way in which convection schemes are constructed.
Ideally, those parts of a convection scheme which require the model time step and those that do not
should live in separate components; until someone writes such a scheme, :py:class:`ClimtImplicitPrognostic` is here
to stay!

:py:class:`ClimtSpectralDynamicalCore` is a generalisation of the standard
:py:class:`~sympl.TimeStepper` in that it allows you to assign Prognostics after object
initialisation and it can function without any Prognostics at all (just adiabatic stepping
forward of the primitive equations).

.. _Sympl: http://sympl.readthedocs.io
.. _documentation: http://sympl.readthedocs.io/en/latest/computation.html
.. _CFL: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

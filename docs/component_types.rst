=============== 
Component Types
===============

Deriving from Sympl_, most components in CliMT are either `Prognostic`, `Diagnostic`
or `Implicit`. You can read more about this schema of model components in Sympl's
documentation_. Unfortunately, not all components used in a typical climate model fit into this
schema. Convection schemes output tendencies (like a `Prognostic`) but require the model timestep
(like an `Implicit`) for various reasons, including to ensure that a vertical CFL_ criterion is met.
Spectral dynamical cores require model state and tendencies in spectral space, which means they
don't play well with `Implicit` components which modify model state in grid space.

To account for this diversity of model components, CliMT introduces two additional entities: `ImplicitPrognostic`
which subclasses `Prognostic` and `SpectralDynamicalCore`, which subclasses `TimeStepper`. The awkward
name `ImplicitPrognostic` mirrors the awkwardness of the way in which convection schemes are constructed.
Ideally, those parts of a convection scheme which require the model time step and those that do not
should live in separate components; until someone writes such a scheme, `ImplicitPrognostic` is here
to stay! `SpectralDynamicalCore` is expected to make it to version 1.0 of CliMT, though it is
currently working in an older version of CliMT.

.. _Sympl: http://sympl.readthedocs.io
.. _documentation: http://sympl.readthedocs.io/en/latest/computation.html
.. _CFL: https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

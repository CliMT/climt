.. highlight:: python

=========================
Configuring climt
=========================

A typical climate model allows the user the following
kinds of configuration:

* Algorithmic Configuration: Changing the "tuning" parameters
  of the algorithms underlying various components. For example,
  the kind of advection used in a dynamical core could be configured.
* Memory/Array Configuration: Deciding the grid shapes, and distribution
  over multiple processors, if using MPI.
* "Physical" Configuration: Choosing the physical constants that are used
  during the simulation. For example, gravitational acceleration or planetary
  rotation rate could be varied.
* Behavioural Configuration: Allows modification of component behaviour. For example,
  radiative transfer components are called only once every N (>> 1) time steps and the
  output is kept constant for the remaining N-1 time steps.
* Compositional Configuration: Describing the components that make up the model, the order
  in which components need to be called, among other things.

Climate models can be configured in many ways, including hard coded configurations, namelists,
shell environment variables. These diverse ways of configuring a climate model make it difficult
to keep track of all configurations and changes made.

climt aims to keep all configuration in the main run script, but separated logically to ensure
the script is still intuitive to read.

Algorithmic Configuration
--------------------------

Configuring the algorithm used by each component is done by various keyword arguments passed
to the component while creating it. See, for example, the documentation for
:py:mod:`climt.RRTMGShortwave`.

Memory/Array Configuration
--------------------------

climt does not yet support MPI, so there is no API yet to handle distributed arrays.
However, the shape of arrays used by a model can be set while calling
:py:func:`climt.get_default_state`. See, for example, the configuration of arrays in a
`GCM`_.

Physical Configuration
----------------------

climt provides an interface to set and reset constants
required by various components. The constants are put into different categories (:py:data:`boltzmann_constant`
is a 'physical constant' whereas :py:data:`planetary_rotation_rate` is a 'planetary constant', for example).

The constants can be reset to their default values so that climt is in a known state at the end of
a simulation. In the future, climt will provide a context manager to clean up modified constants
at the end of a run.

You can read more about this functionality in :ref:`utility_functions`.

Interfacial Configuration
--------------------------

Wrappers are the preferred way of changing the inputs or outputs of a component to make
it apparently work in a different way.

* Piecewise constant output: Computationally expensive modules like radiative transfer
  are sometimes called only once every few timesteps, and the same values is used for
  the intermediate timesteps of a model. For example a GCM with a time step of 10 minutes
  might only call radiation after 1 hour of model time has elapsed. To allow for such
  behaviour, :py:class:`sympl.UpdateFrequencyWrapper` can be used.
  See how this can be used practically in this `example`_.

* TendencyComponent version: Spectral dynamical cores step the model forward in spectral space,
  and therefore, they do not play well with :py:mod:`Stepper`
  components that step forward the model in grid space. Typically, this is handled by
  finite differencing the output of Stepper components and providing them as time tendencies.
  :py:mod:`Stepper` components can be wrapped with :py:class:`sympl.TimeDifferencingWrapper` which
  returns a component which provides the
  time differenced tendencies. The time differencing is done using a first order scheme:

  :math:`\frac{dX}{dt} = (X_{out} - X_{in})/\delta t`.

  See how this is used in the `Grey GCM`_.

* Scaled version: Very often, we perform experiments where we want to study the sensitivity of the simulation
  to a particular quantity or the effect of a certain quantity on the output (mechanism denial).
  This is in some instances done by scaling the quantity or setting it to zero (which
  is also a scaling). To allow for this kind of modification, :py:class:`sympl.ScalingWrapper` can be used. This is a method
  available to all kinds of components (Stepper, TendencyComponent, etc.,). See the documentation for this
  method in the description of the base components in :ref:`component_list`.

Compositional Configuration
----------------------------

This kind of configuration will allow the automatic building of models given certain
components selected by the user.
Currently, the user has to write the script to build the model and run it. It is clear that
a lot of this code is repetitive and can be replaced by an entity (Which will be called
:py:mod:`Federation`).

.. note::
    This functionality is currently unavailble, and will be present in a future version of climt.



.. _GCM: https://github.com/CliMT/climt/blob/a69a23fc2470cc516a41c057976bb3d31ac6f0d7/examples/grey_gcm.py#L56-L59
.. _example: https://github.com/CliMT/climt/blob/a69a23fc2470cc516a41c057976bb3d31ac6f0d7/examples/full_radiation_gcm_energy_balanced.py#L61
.. _Grey GCM: https://github.com/CliMT/climt/blob/a69a23fc2470cc516a41c057976bb3d31ac6f0d7/examples/grey_gcm.py#L48

.. highlight:: python

=========================
Configuring CliMT
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

CliMT aims to keep all configuration in the main run script, but separated logically to ensure
the script is still intuitive to read.

Algorithmic Configuration
--------------------------

Configuring the algorithm used by each component is done by various keyword arguments passed
to the component while creating it. See, for example, the documentation for
:py:mod:`climt.RRTMGShortwave`.

Memory/Array Configuration
--------------------------

CliMT does not yet support MPI, so there is no API yet to handle distributed arrays.
However, the shape of arrays used by a model can be set while calling 
:py:func:`climt.get_default_state`. See, for example, the configuration of arrays in a
`GCM`_.

Physical Configuration
----------------------

CliMT provides an interface to set and reset constants
required by various components. The constants are put into different categories (:py:data:`boltzmann_constant`
is a 'physical constant' whereas :py:data:`planetary_rotation_rate` is a 'planetary constant', for example).

The constants can be reset to their default values so that CliMT is in a known state at the end of
a simulation. In the future, CliMT will provide a context manager to clean up modified constants
at the end of a run.

You can read more about this functionality in :ref:`utility_functions`.

Behavioural Configuration
--------------------------

.. warning::
        This API is currently unstable. Expect the names of methods to change
        in future versions.

Currently, CliMT allows for two kinds of behavioural modification of components.

* Piecewise constant output: Computationally expensive modules like radiative transfer
  are sometimes called only once every few timesteps, and the same values is used for
  the intermediate timesteps of a model. For example a GCM with a time step of 10 minutes
  might only call radiation after 1 hour of model time has elapsed. To allow for such
  behaviour, :py:meth:`climt.ClimtPrognostic.piecewise_constant_version` can be used.
  See how this can be used practically in this `example`_.

* Prognostic version: Spectral dynamical cores step the model forward in spectral space,
  and therefore, they do not play well with :py:mod:`climt.ClimtImplicit`
  components that step forward the model in grid space. Typically, this is handled by
  finite differencing the output of Implicit components and providing them as time tendencies.
  :py:mod:`ClimtImplicit` components have a method :py:mod:`ClimtImplicit.prognostic_version` which 
  returns a component which provides the
  time differenced tendencies. The time differencing is done using a first order scheme:
 
  :math:`\frac{dX}{dt} = (X_{out} - X_{in})/\delta t`.
  
  See how this is used in the `Grey GCM`_.

* Scaled version: Very often, we perform experiments where we want to study the sensitivity of the simulation
  to a particular quantity or the effect of a certain quantity on the output (mechanism denial).
  This is in some instances done by scaling the quantity or setting it to zero (which
  is also a scaling). To allow for this kind of modification, :py:meth:`scaled_version` can be used. This is a method
  available to all kinds of components (Implicit, Prognostic, etc.,).



.. _GCM: https://github.com/CliMT/climt/blob/e171ebef945535f9f82df716da01b4a7c3b1221a/examples/grey_gcm_energy_balanced.py#L51
.. _example: https://github.com/CliMT/climt/blob/e171ebef945535f9f82df716da01b4a7c3b1221a/examples/full_radiation_gcm_energy_balanced.py#L70
.. _Grey GCM: https://github.com/CliMT/climt/blob/5bdac431413f122ae5f46ed4e6610f6a314593c6/examples/grey_gcm_energy_balanced.py#L44 

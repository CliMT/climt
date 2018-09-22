.. highlight:: python

=================
A Realistic Model
=================

As mentioned before, climt includes some components
which returns the a new version of the model state,
and some which return just tendencies.

Since tendencies by themselves are not useful for much
other than plotting, we need to couple them with numerical
integration components to march the model forward in time.
Again, we will use the grey radiation scheme as an example.

The following script is used to obtain the temperature profile
of the atmosphere if no physical process other than radiation
(specifically, grey gas radiation in this example) are present.
The temperature profile obtained is called the **radiative equilibrium** 
profile.

As before, we will create the radiation component and the model state:

.. ipython:: python

   import climt
   import matplotlib.pyplot as plt
   import numpy as np
   # Two new imports
   from sympl import AdamsBashforth
   from datetime import timedelta

   # Create some components
   radiation = climt.GrayLongwaveRadiation()

   # Get a state dictionary filled with required quantities
   # for the components to run
   state = climt.get_default_state([radiation])

We have two new imports, :code:`AdamsBashforth` and :code:`timedelta`.
The former is a numerical `integrator`_ which will step the model forward
in time, and the latter is a standard python module which will be used to represent
the time step of the model.

Now, to create the integrator and the timestep:

.. ipython:: python

    model_time_step = timedelta(hours=1)
    model = AdamsBashforth([radiation])

We now have a model ready to run! The integrator will return the new model
state and any diagnostics that :code:`radiation` has generated. We can then
update the current model state with the new model state and continue to step
the model forward in time:

.. ipython:: python

    for step in range(10):
        diagnostics, new_state = model(state, model_time_step)
        ''' Update state with diagnostics.
        This updated state can be saved if necessary '''
        state.update(diagnostics)
        '''Update state quantities'''
        state.update(new_state)
        '''Update model time'''
        state['time'] += model_time_step
        '''See if the maximum temperature is changing'''
        print(state['time'], ': ', state['air_temperature'].max().values)

And voila, we have a model that actually evolves over time! Many example
scripts that illustrate standard model configurations used in climate
modelling are available in the github `repository`_. These scripts include
examples which setup graphics to view the evolution of the model over time.

.. note::
    A more user friendly API called :code:`Federation` will be available in
    a later version of climt. However, setting up models is easy enough even
    without :code:`Federation` once you get used to the workflow.

.. _integrator: https://en.wikipedia.org/wiki/Linear_multistep_method
.. _repository: https://github.com/CliMT/climt/tree/master/examples

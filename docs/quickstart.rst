.. highlight:: python

===========
Quickstart
===========

Let us start with a simple climt model which is not very useful,
but helps illustrate how to use climt:

.. ipython:: python

   import climt

   # Create some components
   radiation = climt.GrayLongwaveRadiation()

   surface = climt.SlabSurface()

   # Get a state dictionary filled with required quantities
   # for the components to run
   state = climt.get_default_state([radiation, surface])

   # Run components
   tendencies, diagnostics = radiation(state)

   # See output
   tendencies.keys()
   tendencies['air_temperature']

Here, all the essential aspects of creating and running a model
in climt are present:

* Import the :code:`climt` package

* Create one or many components

* Create a state dictionary using :code:`get_default_state`

* Run the components

* Do something with the output

Variables :code:`radiation` and :code:`surface` are two components that we
create. All climt components take a lot of optional arguments: However, by
design, the default options (which are used if you don't specify any arguments)
are meant to be scientifically meaningful.

The variables :code:`state`, :code:`tendencies` and
:code:`diagnostics` are dictionaries which contain quantities
which act either as inputs **to** components or outputs **from** components.

The function :code:`get_default_state()`, if called only with a list of components,
will provide a set of quantities which represent
a single column of the atmosphere. These default values may or may not be
meaningful in your context, so it is best to see what they are and change them
according to your needs.

.. note::
    The square brackets are required in the call to :code:`get_default_state`, even
    if it is one component: :code:`climt.get_default_state([radiation])` is the
    correct syntax.

Building more sophisticated models and running them is merely an extended version
of the above simple example. climt makes heavy use of `Sympl`_, and knowledge of
Sympl is necessary to use climt to its full capabilities. So, do go through Sympl's
docs!

.. _Sympl: http://sympl.readthedocs.io

.. highlight:: python

======================
Interacting with CliMT
======================

As we saw in the Quickstart section, CliMT has two
kinds of entities for the user to interact with:
model components and model state. Here, we will take
a closer look at both these elements.

Model State
------------

The model state is a dictionary whose keys are names of
quantities and values are `xarray`_ DataArrays. To ensure
model scripts are readable not just by specialists, names
of quantities use the descriptive `CF Convention`_. Only
in the case where the CF Convention names are really
unwieldy, like :code:`air_temperature_at_effective_cloud_top_defined_by_infrared_radiation` for
example, we use more convenient names.

DataArrays
are a more human-friendly way of handling numerical arrays.
DataArrays label the dimensions of an array provide
various mathematical functions which can be directly
applied to arrays.

Let's create a 3-d model state to see how useful DataArrays are:

.. ipython:: python

   import climt
   import matplotlib.pyplot as plt
   import numpy as np

   # Create some components
   radiation = climt.GrayLongwaveRadiation()
   condensation = climt.GridScaleCondensation()

We need to tell CliMT what the model dimensions are. Dimension
names, values and units can be arbitrary. Dimension values can be
two-dimensional as well (useful for irregular grids).

.. ipython:: python

   # x coordinate is called 'some_x_coord'
   x = dict(label='some_x_coord',
            values=np.linspace(0, 20, 10),
            dims='some_x_coord',
            units='kilometer')


   # y coordinate is called 'some_y_coord'
   y = dict(label='some_y_coord',
            values=np.linspace(0, 20, 10),
            dims='some_y_coord',
            units='degrees_north')

   # Get a state dictionary filled with required quantities
   # for the components to run
   state = climt.get_default_state([radiation], x=x, y=y)

   state['air_temperature']

CliMT **does not** interpret any of the dimension attributes in
state quantities other than :code:`units`. The values and labels of coordinates
are mainly for users and components. For instance, :py:class:`SimplePhysics`
requires that the y dimension be called :code:`latitude`. So, any
model that uses :py:class:`SimplePhysics` has to label one of the
dimensions as latitude.

As you can see, :code:`air_temperature` has

* a uniform value of 290
* coordinates of some_x_coord, some_y_coord and mid_levels
* units of *degK*, which is the notation used in CliMT (and Sympl) for
  *degrees Kelvin*.

.. note::

    CliMT always names the vertical coordinate as :code:`mid_levels` or :code:`interface_levels`,
    however, the state dictionary will contain a key corresponding to the name
    of the vertical coordinate specified by the user.

You can perform any of the functions `supported`_ by xarray on
the model state quantities.

.. ipython:: python

    state['air_temperature'].sum()

You can access data within DataArrays using their labels:

.. ipython:: python

    state['air_temperature'].loc[dict(mid_levels=[10, 11],
                                     some_x_coord=slice(0,5))]

You can also directly plot DataArrays:

.. ipython:: python

    state['air_temperature'].plot()
    plt.show()

DataArrays are a very powerful way of dealing with array-oriented data, and
you should read more about `xarray`_, and not just for using CliMT!

Model Components
-----------------

Components are representations of physical processes. You can see
all available components in CliMT in the section :ref:`component_list`.

All components take some inputs from the model state, and return **outputs** or
**tendencies** along with diagnostics (if any).

Diagnostics are quantities computed while calculating **outputs** or **tendencies**.
For example, a radiation component calculates heating rates. However, in the process
of calculating these heating rates, it also calculates the radiative flux at each
interface level.

.. ipython:: python

    # These are the tendencies returned by radiation
    radiation.tendencies

    # These are the diagnostics returned by radiation
    radiation.diagnostics

    # These are the outputs returned by condensation
    condensation.outputs

    # These are the diagnostics returned by condensation
    condensation.diagnostics

No component will return **both** outputs and tendencies. The
tendency of a quantity :math:`X` is given by :math:`\frac{dX}{dt}`, and so
the units of a quantity returned as a tendency will always have per second
as as suffix: i.e, if a component is returning :code:`air_temperature` as
a tendency, then its units will be :code:`degK/s`.



.. _xarray: http://xarray.pydata.org

.. _supported: http://xarray.pydata.org/en/stable/computation.html

.. _CF Convention: http://cfconventions.org/Data/cf-standard-names/41/build/cf-standard-name-table.html

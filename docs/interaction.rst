.. highlight:: python

======================
Interacting with climt
======================

As we saw in the Quickstart section, climt has two
kinds of entities for the user to interact with:
model components and model state. Here, we will take
a closer look at both these elements.

Model State
------------

The model state is a dictionary whose keys are names of
quantities and values are `sympl`_ DataArrays. The `sympl`_ DataArray is
a thin wrapper over the `xarray`_ DataArray that makes it units aware. To ensure
model scripts are readable not just by specialists, names
of quantities use the descriptive `CF Convention`_. Only
in the case where the CF Convention names are really
unwieldy, like :code:`air_temperature_at_effective_cloud_top_defined_by_infrared_radiation` for
example, we use more convenient names.

DataArrays are a more human-friendly way of handling numerical arrays.
DataArrays label the dimensions of an array and provide
various mathematical functions which can be directly
applied to arrays. `sympl`_ DataArrays in addition allow conversion
between units, a feature required to allow interoperability between
components which expect inputs in different units.


Let's create a 3-d model state to see how useful DataArrays are:

.. ipython:: python

   import climt
   import matplotlib.pyplot as plt
   import numpy as np

   # Create some components
   radiation = climt.GrayLongwaveRadiation()
   convection = climt.DryConvectiveAdjustment()

We need to tell climt what the model dimensions are. This is done
by the :code:`get_grid` function. This function takes three arguments
which are the number of grid points in the three directions, and
provides a state dictionary containing the definition of a grid.

Passing this grid state dictionary onto :code:`get_default_state` makes
climt aware of the dimensions required by the model:

.. ipython:: python

   grid = climt.get_grid(ny=3, nz=5)

   # Get a state dictionary filled with required quantities
   # for the components to run
   state = climt.get_default_state([radiation, convection], grid_state=grid)

   state['air_temperature']

climt **does not** interpret any of the dimension attributes in
state quantities other than :code:`units`. The values and labels of coordinates
are mainly for users and components. For instance, :py:class:`SimplePhysics`
requires that the y dimension be called :code:`latitude`. So, any
model that uses :py:class:`SimplePhysics` has to label one of the
dimensions as latitude.

As you can see, :code:`air_temperature` has

* a uniform value of 290
* coordinates of latitude and mid_levels
* units of *degK*, which is the notation used in climt (and Sympl) for
  *degrees Kelvin*.

It is also fairly easy to change units. The :py:func:`.to_units()` method can
be used as below to return a DataArray with the equivalent temperature in degrees Farenheit:

.. ipython:: python

    state['air_temperature'].to_units('degF')


.. note::

    climt always names the vertical coordinate as :code:`mid_levels` or :code:`interface_levels`,
    however, the state dictionary will contain a key corresponding to the name
    of the vertical coordinate specified by :code:`get_grid`.

As mentioned previously, DataArrays are a user-friendly way of handling numerical or numpy
arrays. The numpy array underlying any DataArray is easily accessed using the :code:`values`
attribute:

.. ipython:: python

    type(state['air_temperature'].values)

and can also be modified easily:

.. ipython:: python

    state['air_temperature'].values[:] = 291

The right hand side can also be any numpy array, as long as it has the same dimensions (or can
be broadcasted to the same dimensions) as the
current numpy array.

.. note::

    It is recommended to use the syntax :code:`...values[:] = ...` rather than :code:`...values =
    ...`, as the former modifies the numpy array in-place. In either case, DataArrays check to
    ensure the dimensions (or shape) of the new data matches with the current dimensions. 

You can perform any of the functions `supported`_ by xarray on
the model state quantities.

.. ipython:: python

    state['air_temperature'].sum()


You can also directly plot DataArrays:

.. ipython:: python

    state['air_temperature'].plot()

DataArrays are a very powerful way of dealing with array-oriented data, and
you should read more about `xarray`_, and not just for using climt!

Model Components
-----------------

Components are representations of physical processes. You can see
all available components in climt in the section :ref:`component_list`.

All components take some inputs from the model state, and return **outputs** or
**tendencies** along with diagnostics (if any).

Diagnostics are quantities computed while calculating **outputs** or **tendencies**.
For example, a radiation component calculates heating rates. However, in the process
of calculating these heating rates, it also calculates the radiative flux at each
interface level.

.. ipython:: python

    # These are the tendencies returned by radiation
    radiation.tendency_properties

    # These are the diagnostics returned by radiation
    radiation.diagnostic_properties

    # These are the outputs returned by convection
    convection.output_properties

    # convection returns no diagnostics
    convection.diagnostic_properties

No component will return **both** outputs and tendencies. The
tendency of a quantity :math:`X` is given by :math:`\frac{dX}{dt}`, and so
the units of a quantity returned as a tendency will always have per second
as as suffix: i.e, if a component is returning :code:`air_temperature` as
a tendency, then its units will be :code:`degK/s`.



.. _xarray: http://xarray.pydata.org

.. _sympl: http://sympl.readthedocs.io

.. _supported: http://xarray.pydata.org/en/stable/computation.html

.. _CF Convention: http://cfconventions.org/Data/cf-standard-names/41/build/cf-standard-name-table.html

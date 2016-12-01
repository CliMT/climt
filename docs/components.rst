==========
Components
==========

In CliMT, "components" is a loose term to refer to :py:class:`Prognostic`,
:py:class:`Diagnostic`, :py:class:`Implicit`, and :py:class:`Monitor` objects.
These are not objects of their own (you can't run ``prog = Prognostic()``). They
define interfaces for subclasses to follow.

Each of them, once initialized, can be passed in a current model state.
``Prognostic`` objects use the state to return tendencies and diagnostics at
the current time. ``Diagnostic`` objects return only diagnostics from the
current time. ``Implicit`` objects will take in a timestep along with the
state, and then return the next state as well as modifying the current state
to include more diagnostics (it is similar to a :py:class:`TimeStepper` in how
it is called). And ``Monitor`` objects will store the current state in some
way, whether it is by displaying the new state on a plot that is shown to the
user, updating information on a web server, or saving the state to a file.

These classes themselves (listed in the previous paragraph) are not ones you
can initialize (e.g. there is no one 'prognostic' scheme), but they each have
subclasses defined in CliMT, which can all be used in the ways described in
this section.

Ensuring consistency
--------------------

To make sure that a given state dictionary is valid to pass in to a
component, you can write (for example, with a ``Prognostic``):

.. code-block:: python

    prognostic = RRTMRadiation()
    prognostic.ensure_state_is_valid_input(state)

where ``state`` is your state dictionary. This call will raise an
:py:class:`InvalidStateException` if the input state is not valid. This is good,
because if you call a ``Prognostic`` with an invalid input, it may give results
which are not scientifically valid. For instance, if the object requires the
horizontal grid to have constant spacing, it may assume this to be true
without checking when it is called (for performance reasons).

The ``ensure_state_is_valid_input`` method exists on all component classes, and
is called in the same way as shown above.

Initializing a Component
~~~~~~~~~~~~~~~~~~~~~~~~

You can create a component object, such as ``RRTMRadiation`` like so:

.. code-block:: python

    radiation = RRTMRadiation()

All component objects (``Prognostic``, ``Diagnostic``, ``Implicit``,
and ``Monitor``) should have sensible defaults that are used when
no options are passed in on initialization. Each object can have its own
options, however. For example, to compute only longwave radiation you might
write:

.. code-block:: python

    radiation = RRTMRadiation(do_shortwave=False, do_longwave=True)

The options available ('keyword arguments' or 'kwargs') will depend on the
particular object being initialized.

Prognostic
----------

As stated above, :py:class:`Prognostic` objects use the state to return
tendencies and diagnostics at the current time. In a full model, the tendencies
are used by a time stepping scheme (in CliMT, a :py:class:`TimeStepper`) to
determine the values of quantities at the next time.

You can call a ``Prognostic`` directly to get diagnostics and tendencies like
so:

.. code-block:: python

    radiation = RRTMRadiation()
    diagnostics, tendencies = radiation(state)

``diagnostics`` and ``tendencies`` in this case will both be dictionaries,
similar to ``state``. Even if the ``Prognostic`` being called does not compute
any diagnostics, it will still return an empty diagnostics dictionary.

Usually, you will call a Prognostic object through a :py:class:`TimeStepper`
that uses it to determine values at the next timestep.

Diagnostic
----------

:py:class:`Diagnostic` objects use the state to return quantities
('diagnostics') from the same timestep as the input state. You can call a
``Diagnostic`` directly to get diagnostic quantities like so:

.. code-block:: python

    diagnostic_component = MyDiagnostic()
    diagnostics = diagnostic_component(state)

Instead of returning a new dictionary with the additional diagnostic quantities,
a ``Diagnostic`` can update the state dictionary in-place with the new
quantities. You do this like so:

.. code-block:: python

    diagnostic_component = MyDiagnostic()
    diagnostic_component.update_state(state)

The ``update_state`` call has the advantage that it will automatically check to
see if it is overwriting any quantities already present in state, and will
raise a :py:class:`SharedKeyException` before doing so. This ensures you
don't have multiple pieces of code trying to output the same diagnostic, with
one overwriting the other.


Implicit
--------

:py:class:`Implicit` objects use a state and a timestep to return the next
state, and update the input state with any relevant diagnostic quantities. You
can call an Implicit object like so:

.. code-block:: python

    from datetime import timedelta
    implicit = MyImplicit()
    timestep = timedelta(minutes=10)
    next_state = implicit(state, timestep)

Following the ``implicit`` call, ``state`` will have been modified in-place to
include any diagnostics produced by the ``Implicit`` component for the timestep
of the input state.

Monitor
-------

:py:class:`Monitor` objects can store states in some way, whether it is by
displaying the new state on a plot that is shown to the user, updating
information on a web server, or saving the state to a file. They are called
like so:

.. code-block:: python

    monitor = MyMonitor()
    monitor.store(state)

The ``Monitor`` will take advantage of the 'time' key in the ``state``
dictionary in order to determine the model time of the state. This is
particularly important for a ``Monitor`` which outputs a series of states to
disk.

Timestepping
============

:py:class:`climt.TimeStepper` objects use :py:class:`climt.Diagnostic` and
:py:class:`climt.Prognostic` objects to step a model state forward in time.

Initialization
--------------

:py:class:``climt.TimeStepper`` objects are initialized using lists of
:py:class:``climt.Prognostic`` and :py:class:``climt.Diagnostic`` objects.

.. code-block:: python

    prognostic_list = [MyPrognostic(), MyOtherPrognostic()]
    diagnostic_list = [MyDiagnostic()]
    time_stepper = AdamsBashforth(prognostic_list, diagnostic_list)

If no :py:class:``climt.Diagnostic`` objects are being used, its list does not
need to be given.

.. code-block:: python

    time_stepper = AdamsBashforth([MyPrognostic()])

It is possible to pass in an empty list for either the
:py:class:``climt.Prognostic`` or :py:class:``climt.Diagnostic`` objects, but
this is likely not a useful thing to do.

Usage
-----

:py:class:``climt.TimeStepper`` objects support the same consistency checks
as components (the ``ensure_state_is_valid_input`` method) outlined in
:ref:`consistency`.

Once initialized, a :py:class:``climt.TimeStepper`` object has a very similar
interface to the :py:class:`Implicit` object.

.. code-block:: python

    from datetime import timedelta
    time_stepper = AdamsBashforth([MyPrognostic()])
    timestep = timedelta(minutes=10)
    next_state = time_stepper(state, timestep)

Following the ``time_stepper`` call, ``state`` will have been modified
in-place to include any diagnostics produced by the
:py:class:``climt.Implicit`` component for the timestep of the input state.

This is important, so we'll repeat it:
**the input state will be modified by the call to the
:py:class:``climt.TimeStepper``**. In fact, it is only after calling the
:py:class:``climt.TimeStepper`` and getting the new state that the
previous/current state will have all the diagnostic quantities
produced by the :py:class:``climt.Diagnostic`` and
:py:class:``climt.Prognostic`` objects. This means you will sometimes want to
pass ``state`` to your :py:class:``climt.Monitor`` objects *after* calling
the :py:class:``climt.TimeStepper`` and getting ``next_state``.

.. autoclass:: climt.TimeStepper
    :members:
    :special-members:
    :exclude-members: __weakref__,__metaclass__

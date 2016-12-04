==========
Composites
==========

There are a set of objects in CliMT that wrap multiple components into a single
object so they can be called as if they were one component. There is one each
for :py:class:`climt.Prognostic`, :py:class:`climt.Diagnostic`, and
:py:class:`climt.Monitor`. These can be used to simplify code, so that
the way you call a list of components is the same as the way you would
call a single component. For example, *instead* of writing:

.. code-block:: python

    prognostic_list = [
        MyPrognostic(),
        MyOtherPrognostic(),
        YetAnotherPrognostic(),
    ]
    all_diagnostics = {}
    total_tendencies = {}
    for prognostic_component in prognostic_list:
        tendencies, diagnostics = prognostic_component(state)
        # this should actually check to make sure nothing is overwritten,
        # but this code does not
        all_tendencies.update(tendencies)
        for key in tendencies.keys():
            if key not in total_tendencies:
                total_tendencies[key] = tendencies[key]
            else:
                total_tendencies[key] += tendencies[key]

You could write:

.. code-block:: python

    prognostic_composite = PrognosticComposite([
        MyPrognostic(),
        MyOtherPrognostic(),
        YetAnotherPrognostic(),
    ])
    tendencies, diagnostics = prognostic_composite(state)

This second call is much cleaner. It will also automatically detect whether
multiple components are trying to write out the same diagnostic, and raise
an exception if that is the case (so no results are being silently
overwritten). You can get similar simplifications for
:py:class:`climt.Diagnostic` and :py:class:`climt.Monitor`.

API Reference
-------------

.. autoclass:: climt.PrognosticComposite
    :members:
    :special-members:
    :exclude-members: __weakref__,__metaclass__

.. autoclass:: climt.DiagnosticComposite
    :members:
    :special-members:
    :exclude-members: __weakref__,__metaclass__


.. autoclass:: climt.MonitorComposite
    :members:
    :special-members:
    :exclude-members: __weakref__,__metaclass__

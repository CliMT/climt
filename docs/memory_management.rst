=================
Memory Management
=================

Arrays
------

If possible, you should try to be aware of when there are two code references
to the same in-memory array. This can help avoid some common bugs. Let's start
with an example. Say you create a ConstantTendencyComponent object like so::

    >>> import numpy as np
    >>> from climt import ConstantTendencyComponent, DataArray
    >>> array = DataArray(
            np.ones((5, 5, 10)),
            dims=('lon', 'lat', 'lev'), attrs={'units': 'K/s'})
    >>> tendencies = {'air_temperature': array}
    >>> prognostic = ConstantTendencyComponent(tendencies)

This is all fine so far. But it's important to know that now ``array`` is the
same array stored inside ``prognostic``::

    >>> out_tendencies, out_diagnostics = prognostic({})
    >>> out_tendencies['air_temperature'] is array  # same place in memory
    True

So if you were to modify ``array``, it would *change the output given by
prognostic*::

    >>> array[:] = array[:] * 5.
    >>> out_tendencies, out_diagnostics = prognostic({})
    >>> out_tendencies['air_temperature'] is array
    True
    >>> np.all(out_tendencies['air_temperature'].values == array.values)
    True

When in doubt, assume that any array you put into a component when it is
initialized should not be modified any more, unless changing the values in the
component is intentional. Below is some less (but potentially) useful
information for those interested.

If instead of modifying ``array``, you make a new array for the python variable
``array`` to refer to, it doesn't modify the array in ``prognostic``::

    >>> array = array * 5.
    >>> out_tendencies, out_diagnostics = prognostic({})
    >>> out_tendencies['air_temperature'] is array
    False
    >>> np.all(out_tendencies['air_temperature'].values == array.values)
    False

This is because having the ``[:]`` on the left hand side of the assignment
operator ``\=`` tells python that you want to modify the existing memory of the
array on the left hand side. More precisely, having ``array =`` tells python
that you want to change what the variable ``array`` refers to, and set it to
be the thing on the right hand side, while ``array[:] =`` tells python to
call the ``__setitem__(key, value)`` method of ``array`` with the contents
of the square parentheses as the key and the right hand side as the value.

Interestingly, ``array = array * 5.`` has different behavior from
``array *= 5.``. The first one will change what ``array`` refers to, as before,
while the second one will modify ``array`` in-place without changing the
reference. All similarly written operations (``-=``, ``+=``, ``/=``, etc.) are
in-place operations. When you want to avoid copying data, ``array *= 5.`` is
better since the values of the array will be modified where they already are
in memory, instead of allocating an entirely new array.

Dictionaries
------------

Unlike arrays, the dictionary containers are copied when passed in. Copying
dictionaries is fairly cheap, since the new dictionary will still refer to the
same values (arrays) as before, and all that has to be copied is the key-value
pairs::

    >>> tendencies['new_quantity'] = array
    >>> out_tendencies, out_diagnostics = prognostic({})
    >>> 'new_quantity' in out_tendencies.keys()
    False


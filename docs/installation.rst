.. highlight:: shell

============
Installation
============

Dependencies
------------

CliMT depends on some easily installable libraries. For
an installation from source, it also requires that C and fortran
compilers be installed.

On Ubuntu Linux, for example, these dependencies can be
installed by:

.. code-block:: console

    $ sudo apt-get install gcc
    $ sudo apt-get install gfortran

On Mac OSX, you will need to do the following:

.. code-block:: console

    $ brew install gcc
    $ export CC=gcc-x

where :code:`gcc-x` is the name of the C compiler that `Homebrew`_
installs. Exporting the name of the compiler is essential on Mac since the
default compiler that ships with Mac (called :code:`gcc`, but is actually a
different compiler) cannot
compile OpenMP programs, like the dynamical core in CliMT.


.. _Homebrew: https://brew.sh/


Stable release
--------------

Once you have the dependencies installed,
you can install CliMT by running this command in your terminal:

.. code-block:: console

    $ pip install climt

This is the preferred method to install CliMT, as it will always install the most recent stable release.
On Ubuntu Linux, you might need to prefix the above command with :code:`sudo`.

.. NOTE::
    for the moment, a pip install compiles from source. In the future, we plan to have
    binary distributions of CliMT.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for CliMT can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/CliMT/climt

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/CliMT/climt/tarball/master


Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ pip install -r requirements_dev.txt
    $ python setup.py install

Both commands may require the use of :code:`sudo`.

.. _Github repo: https://github.com/CliMT/climt
.. _tarball: https://github.com/CliMT/climt/tarball/master

=====
CliMT
=====


.. image:: https://img.shields.io/pypi/v/climt.svg
    :target: https://pypi.python.org/pypi/climt
    :alt: PyPI

.. image:: https://img.shields.io/travis/CliMT/climt.svg
    :target: https://travis-ci.org/CliMT/climt
    :alt: Continuous Integration

.. image:: https://img.shields.io/codecov/c/github/CliMT/climt.svg
    :target: https://travis-ci.org/CliMT/climt
    :alt: Coverage

.. image:: https://readthedocs.org/projects/climt/badge/?version=latest
    :target: https://climt.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


CliMT is a Toolkit for building Earth system models in Python. CliMT stands for *Climate Modelling
and diagnostics Toolkit* -- it is meant both for creating models and for generating diagnostics
(radiative fluxes for an atmospheric column, for example). However, since it might eventually
include model components for purposes other than climate modelling (local area models, large-eddy
simulation), we prefer to keep the abbreviation un-expanded!

CliMT hopes to enable researchers to easily perform online analysis and make
modifications to existing models by increasing the ease with which models
can be understood and modified. It also enables educators to write
accessible models that serve as an entry point for students into Earth
system modeling, while also containing state-of-the-art components.

Initially CliMT contains only components for the atmosphere, and does not yet
include a coupler. But there are plans to extend CliMT to a fully coupled Earth
system model in the future. The toolkit is also written in such a way that it
could enable the development of non-climate models (e.g. weather prediction,
large-eddy simulation). To do so requires only that the prognostic and
diagnostic schemes are wrapped into the correct Python-accessible interface.

CliMT builds on Sympl_, which provides the base classes and  array and constants handling
functionality. Thanks to Sympl_ and Pint_, CliMT is also a fully units aware model. It is
useful to know how Sympl_ works to use CliMT better. Read more about Sympl_ at
https://sympl.readthedocs.io.

* Free software: BSD license
* Documentation: https://climt.readthedocs.io.

Installation
-------------

The `documentation`_ has installations to install CliMT.


Features
--------

* CliMT is fully units-aware!
* Uses the xarray_ `DataArray` abstraction to build self describing model arrays. 
* Provides different levels of abstraction towards building a climate model.
* Like Sympl_, CliMT consciously uses descriptive names in the user API to ensure
  model scripts are self-documenting.
* Allows for quick prototyping of earth system model components.
* Provides a clean and convenient interface to add new components.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _Sympl: https://github.com/mcgibbon/sympl
.. _Pint: https://pint.readthedocs.io
.. _xarray: http://xarray.pydata.org
.. _documentation: http://climt.readthedocs.io/en/latest/installation.html

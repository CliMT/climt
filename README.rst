=====
climt
=====


.. image:: https://img.shields.io/pypi/v/climt.svg
    :target: https://pypi.python.org/pypi/climt
    :alt: PyPI

.. image:: https://img.shields.io/travis/climt/climt.svg
    :target: https://travis-ci.org/climt/climt
    :alt: Continuous Integration

.. image:: https://ci.appveyor.com/api/projects/status/h9ayx22cxyfwh5rh?svg=true
    :target: https://ci.appveyor.com/project/JoyMonteiro/climt
    :alt: Continuous Integration

.. image:: https://img.shields.io/codecov/c/github/climt/climt.svg
    :target: https://travis-ci.org/climt/climt
    :alt: Coverage

.. image:: https://readthedocs.org/projects/climt/badge/
    :target: https://climt.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://zenodo.org/badge/74854230.svg
    :target: https://zenodo.org/badge/latestdoi/74854230
    :alt: Zenodo DOI


.. image:: ./docs/climt_logo.jpg
    :height: 512px
    :width: 512px
    :align: center

**climt** is a Toolkit for building Earth system models in Python. climt stands for *Climate Modelling
and Diagnostics Toolkit* -- it is meant both for creating models and for generating diagnostics
(radiative fluxes for an atmospheric column, for example). However, since it might eventually
include model components for purposes other than climate modelling (local area models, large-eddy
simulation), we prefer to keep the abbreviation un-expanded!

climt hopes to enable researchers to easily perform online analysis and make
modifications to existing models by increasing the ease with which models
can be understood and modified. It also enables educators to write
accessible models that serve as an entry point for students into Earth
system modeling, while also containing state-of-the-art components.

Initially climt contains only components for the atmosphere, and does not yet
include a coupler. But there are plans to extend climt to a fully coupled Earth
system model in the future. The toolkit is also written in such a way that it
could enable the development of non-climate models (e.g. weather prediction,
large-eddy simulation). To do so requires only that the prognostic and
diagnostic schemes are wrapped into the correct Python-accessible interface.

climt builds on sympl_, which provides the base classes and  array and constants handling
functionality. Thanks to sympl_ and Pint_, climt is also a fully units aware model. It is
useful to know how sympl_ works to use climt better. Read more about sympl_ at
https://sympl.readthedocs.io.

* Free software: BSD license
* Documentation: https://climt.readthedocs.io.

Installation
-------------

climt can be installed directly from the python package index using pip.

    pip install climt

should work on most systems. From version 0.9.2 onwards, this command will
install binary wheels, eliminating the requirement of a compiler on your
system.

Detailed instructions for Mac and Linux systems are available in the `documentation`_.

Features
--------

* climt is fully units-aware!
* Uses the xarray_ `DataArray` abstraction to build self describing model arrays. 
* Provides different levels of abstraction towards building a climate model.
* Like sympl_, climt consciously uses descriptive names in the user API to ensure
  model scripts are self-documenting.
* Allows for quick prototyping of earth system model components.
* Provides a clean and convenient interface to add new components.

Citing climt
------------

If you use climt in your research, please cite the following paper documenting sympl_ and climt

    https://www.geosci-model-dev.net/11/3781/2018/

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _sympl: https://github.com/mcgibbon/sympl
.. _Pint: https://pint.readthedocs.io
.. _xarray: http://xarray.pydata.org
.. _documentation: http://climt.readthedocs.io/en/latest/installation.html

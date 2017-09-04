.. _component_list:
.. currentmodule:: climt

==========
Components
==========

This page documents the different components available through CliMT.

Base Components
----------------

.. autosummary::
    :toctree: generated/

    ClimtPrognostic

    ClimtDiagnostic

    ClimtImplicit

    ClimtImplicitPrognostic

    ClimtSpectralDynamicalCore

Dynamics
---------

.. autosummary::
    :toctree: generated/

    GfsDynamicalCore
    GfsDynamicalCore.__call__

Radiation
---------

.. autosummary::
    :toctree: generated/

    RRTMGLongwave
    RRTMGLongwave.__call__

    RRTMGShortwave
    RRTMGShortwave.__call__

    GrayLongwaveRadiation
    GrayLongwaveRadiation.__call__

    Frierson06LongwaveOpticalDepth
    Frierson06LongwaveOpticalDepth.__call__

Convection
----------

.. autosummary::
    :toctree: generated/

    EmanuelConvection
    EmanuelConvection.__call__

Surface Processes
-----------------

.. autosummary::
    :toctree: generated/

    SimplePhysics
    SimplePhysics.__call__

    SlabSurface
    SlabSurface.__call__

Ice and Snow
------------

.. autosummary::
    :toctree: generated/

    IceSheet
    IceSheet.__call__

Test Cases
-----------

.. autosummary::
    :toctree: generated/

    HeldSuarez
    HeldSuarez.__call__

    DcmipInitialConditions
    DcmipInitialConditions.__call__

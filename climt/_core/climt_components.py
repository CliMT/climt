from sympl import Prognostic
from datetime import timedelta


class ImplicitPrognostic(Prognostic):
    """
    The base class used mainly for convection schemes.

    Convection schemes return tendencies, but require the current model time step
    to work. This can be worked around by simply adding an attribute storing the
    time step in :code:`Prognostic`, but then we lose the distinction between these
    kind of schemes and pure Prognostics like radiative transfer codes. Depending on
    availability of time and the possibility of resolving this situation (by
    decomposing convection schemes into multiple components?), **this class might eventually
    be deprecated.**

    """

    current_time_step = timedelta(seconds=100.)

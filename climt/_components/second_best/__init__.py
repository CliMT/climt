try:
    from .component import SecondBEST

    __all__ = ["SecondBEST"]
except ImportError:
    # component.py lands in a later task; keep the package importable so
    # subpackages (e.g. processes/) can be developed/tested independently.
    __all__ = []

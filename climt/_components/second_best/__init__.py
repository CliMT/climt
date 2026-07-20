try:
    from .component import SecondBEST

    __all__ = ["SecondBEST"]
except ModuleNotFoundError as _e:
    # TODO(task-9): remove this guard once component.py exists — the real
    # import above should then run unconditionally.
    # component.py lands in a later task; keep the package importable so
    # subpackages (e.g. processes/) can be developed/tested independently.
    # Narrowly catch only the missing-component-module case, so a genuine
    # import error inside a future component.py still surfaces.
    if _e.name != __name__ + ".component":
        raise
    __all__ = []

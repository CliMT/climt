import inspect
import sys

import climt
from climt._core import fortran, registry


def by_name(components):
    return {c.name: c for c in components}


def test_every_exported_component_is_listed():
    listed = by_name(climt.list_components())
    for cls in climt._components.__all__:
        assert listed[cls.__name__].source == "climt"
    assert listed["GrayLongwaveRadiation"] == registry.ComponentInfo(
        "GrayLongwaveRadiation", "TendencyComponent", True, "", "climt")


def test_fortran_backed_table_matches_the_code():
    """A component is in FORTRAN_BACKED exactly when its module asks the loader for an extension."""
    for cls in climt._components.__all__:
        uses_loader = "load_extension(" in inspect.getsource(sys.modules[cls.__module__])
        assert uses_loader == (cls.__name__ in registry.FORTRAN_BACKED), cls.__name__
    assert set(registry.FORTRAN_BACKED.values()) == set(fortran.MODULES)


def test_unavailable_components_say_what_they_need(monkeypatch):
    status = dict(fortran.fortran_status(), state="stub", modules={m: False for m in fortran.MODULES})
    monkeypatch.setattr(registry, "fortran_status", lambda: status)
    listed = by_name(climt.list_components())
    for name in registry.FORTRAN_BACKED:
        assert listed[name].available is False and listed[name].requires == "climt-fortran"
    assert listed["CorkLongwaveRadiation"].available is True


def test_entry_points_are_listed_without_being_loaded(monkeypatch):
    class FakeDist:
        name = "some-package"

    class FakeEntryPoint:
        name = "SomeComponent"
        value = "some_package.components:SomeComponent"
        dist = FakeDist()

        def load(self):
            raise AssertionError("entry points must not be loaded")

    monkeypatch.setattr(registry, "_entry_points", lambda: [FakeEntryPoint()])
    record = by_name(climt.list_components())["SomeComponent"]
    assert record == registry.ComponentInfo(
        "SomeComponent", "external", True, "", "some-package (some_package.components:SomeComponent)")


def test_repr_is_a_table_with_a_numba_line():
    text = repr(climt.list_components())
    lines = text.splitlines()
    assert lines[0].split()[:4] == ["name", "kind", "available", "requires"]
    assert any(line.startswith("RRTMGLongwave") for line in lines)
    assert lines[-1].startswith("numba: ")
